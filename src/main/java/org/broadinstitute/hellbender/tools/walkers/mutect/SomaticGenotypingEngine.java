package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.*;
import org.apache.commons.collections4.ListUtils;
import org.apache.commons.lang.mutable.MutableDouble;
import org.apache.commons.lang.mutable.MutableInt;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.analysis.integration.gauss.GaussIntegrator;
import org.apache.commons.math3.analysis.integration.gauss.GaussIntegratorFactory;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.log4j.Logger;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.AFCalculator;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.AFCalculatorProvider;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.*;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.broadinstitute.hellbender.utils.MathUtils;

import static org.broadinstitute.hellbender.utils.OptimizationUtils.argmax;

public class SomaticGenotypingEngine extends AssemblyBasedCallerGenotypingEngine {

    public static final String IN_COSMIC_VCF_ATTRIBUTE = "IN_COSMIC";
    public static final String IN_DBSNP_VCF_ATTRIBUTE = "IN_DBSNP";
    public static final String IN_PON_VCF_ATTRIBUTE = "IN_PON";
    public static final String NORMAL_ARTIFACT_LOD_ATTRIBUTE = "N_ART_LOD";
    public static final String STRAND_ARTIFACT_POSTERIOR_PROBABILITIES_KEY = "SA_POST_PROB";
    public static final String STRAND_ARTIFACT_MAP_ALLELE_FRACTIONS_KEY = "SA_MAP_AF";

    private final M2ArgumentCollection MTAC;

    private final String tumorSampleName;
    private final String matchedNormalSampleName;
    final boolean hasNormal;

    //Mutect2 does not run in GGA mode
    private static final List<VariantContext> NO_GIVEN_ALLELES = Collections.emptyList();

    private static final OptionalDouble NO_FIXED_TUMOR_ALT_FRACTION = OptionalDouble.empty();
    private static final OptionalDouble GERMLINE_HET_ALT_FRACTION = OptionalDouble.of(0.5);

    // {@link GenotypingEngine} requires a non-null {@link AFCalculatorProvider} but this class doesn't need it.  Thus we make a dummy
    private static AFCalculatorProvider DUMMY_AF_CALCULATOR_PROVIDER = new AFCalculatorProvider() {
        @Override
        public AFCalculator getInstance(final int ploidy, final int maximumAltAlleles) { return null; }
    };

    private final static Logger logger = Logger.getLogger(SomaticGenotypingEngine.class);

    @Override
    protected String callSourceString() {
        return "M2_call";
    }

    public SomaticGenotypingEngine(final SampleList samples,
                                   final M2ArgumentCollection MTAC,
                                   final String tumorSampleName,
                                   final String matchedNormalSampleName) {
        super(MTAC, samples, DUMMY_AF_CALCULATOR_PROVIDER, !MTAC.doNotRunPhysicalPhasing);
        this.MTAC = MTAC;
        this.tumorSampleName = tumorSampleName;
        this.matchedNormalSampleName = matchedNormalSampleName;
        hasNormal = matchedNormalSampleName != null;

        final double errorProbability = QualityUtils.qualToErrorProb(MTAC.POWER_CONSTANT_QSCORE);
    }

    /**
     * Main entry point of class - given a particular set of haplotypes, samples and reference context, compute
     * genotype likelihoods and assemble into a list of variant contexts and genomic events ready for calling
     *
     * The list of samples we're working with is obtained from the readLikelihoods
     * @param readLikelihoods                       Map from reads->(haplotypes,likelihoods)
     * @param perSampleFilteredReadList              Map from sample to reads that were filtered after assembly and before calculating per-read likelihoods.
     * @param activeRegionWindow                     Active window
     *
     * @return                                       A CalledHaplotypes object containing a list of VC's with genotyped events and called haplotypes
     */
    public CalledHaplotypes callMutations(
            final ReadLikelihoods<Haplotype> readLikelihoods,
            final Map<String, List<GATKRead>> perSampleFilteredReadList,
            final AssemblyResultSet assemblyResultSet,
            final ReferenceContext referenceContext,
            final SimpleInterval activeRegionWindow,
            final FeatureContext featureContext,
            final SAMFileHeader header) {
        Utils.nonNull(readLikelihoods, "likelihoods are null");
        Utils.validateArg(readLikelihoods.numberOfSamples() > 0, "likelihoods have no samples");
        Utils.nonNull(activeRegionWindow, "activeRegionWindow is null");
        Utils.validateArg(readLikelihoods.samples().contains(tumorSampleName), "readLikelihoods does not contain the tumor sample ");

        final List<Haplotype> haplotypes = readLikelihoods.alleles();

        // update the haplotypes so we're ready to call, getting the ordered list of positions on the reference
        // that carry events among the haplotypes
        final List<Integer> startPosKeySet = decomposeHaplotypesIntoVariantContexts(haplotypes, assemblyResultSet.getFullReferenceWithPadding(), assemblyResultSet.getPaddedReferenceLoc(), NO_GIVEN_ALLELES).stream()
                .filter(loc -> activeRegionWindow.getStart() <= loc && loc <= activeRegionWindow.getEnd())
                .collect(Collectors.toList());

        final Set<Haplotype> calledHaplotypes = new HashSet<>();
        final List<VariantContext> returnCalls = new ArrayList<>();

        for( final int loc : startPosKeySet ) {
            final List<VariantContext> eventsAtThisLoc = getVCsAtThisLocation(haplotypes, loc, NO_GIVEN_ALLELES);
            final VariantContext mergedVC = AssemblyBasedCallerUtils.makeMergedVariantContext(eventsAtThisLoc);
            if( mergedVC == null ) {
                continue;
            }
            final Map<Allele, List<Haplotype>> alleleMapper = createAlleleMapper(eventsAtThisLoc, mergedVC, loc, haplotypes);

            // converting ReadLikelihoods<Haplotype> to ReadLikeliHoods<Allele>
            ReadLikelihoods<Allele> readAlleleLikelihoods = readLikelihoods.marginalize(alleleMapper,
                    new SimpleInterval(mergedVC).expandWithinContig(ALLELE_EXTENSION, header.getSequenceDictionary()));

            //TODO: downsampling goes here but gatk does not have somatic downsampling as of 12/6/2016

            filterOverlappingReads(readAlleleLikelihoods, tumorSampleName, mergedVC.getReference(), loc, false);
            if (hasNormal) {
                filterOverlappingReads(readAlleleLikelihoods, matchedNormalSampleName, mergedVC.getReference(), loc, true);
            }

            final PerAlleleCollection<Double> tumorLods = getHetGenotypeLogOdds(readAlleleLikelihoods, false, tumorSampleName);
            final Optional<PerAlleleCollection<Double>> normalLods = !hasNormal ? Optional.empty() :
                    Optional.of(getHetGenotypeLogOdds(readAlleleLikelihoods, true, matchedNormalSampleName));
            final Optional<PerAlleleCollection<Double>> normalArtifactLods = !hasNormal ? Optional.empty() :
                    Optional.of(getHetGenotypeLogOdds(readAlleleLikelihoods, false, matchedNormalSampleName));

            final List<Allele> somaticAltAlleles = mergedVC.getAlternateAlleles().stream()
                    .filter(allele -> tumorLods.getAlt(allele) > MTAC.TUMOR_LOD_THRESHOLD)
                    .filter(allele -> hasNormal ? normalLods.get().getAlt(allele) > MTAC.NORMAL_LOD_THRESHOLD : true)
                    .collect(Collectors.toList());

            if (somaticAltAlleles.isEmpty()) {
                continue;
            }

            final VariantContextBuilder callVcb = new VariantContextBuilder(mergedVC);
            callVcb.attribute(GATKVCFConstants.TUMOR_LOD_KEY, somaticAltAlleles.stream().mapToDouble(tumorLods::getAlt).toArray());

            if (hasNormal) {
                callVcb.attribute(GATKVCFConstants.NORMAL_LOD_KEY, somaticAltAlleles.stream().mapToDouble(a -> normalLods.get().getAlt(a)).toArray());
                callVcb.attribute(NORMAL_ARTIFACT_LOD_ATTRIBUTE, somaticAltAlleles.stream().mapToDouble(a -> normalArtifactLods.get().getAlt(a)).toArray());
            }

            final PerAlleleCollection<Double> tumorAlleleFractions = getAlleleFractions(readAlleleLikelihoods, tumorSampleName);

            //TODO: multiple alt alleles -- per-allele strand bias
            final List<Allele> allSomaticAlleles = ListUtils.union(Arrays.asList(mergedVC.getReference()), somaticAltAlleles);
            addStrandArtifactAnnotations(readAlleleLikelihoods, callVcb, allSomaticAlleles);

            if (!featureContext.getValues(MTAC.cosmicFeatureInput, loc).isEmpty()) {
                callVcb.attribute(IN_COSMIC_VCF_ATTRIBUTE, true);
            }

            if (!featureContext.getValues(MTAC.dbsnp.dbsnp, loc).isEmpty()) {
                callVcb.attribute(IN_DBSNP_VCF_ATTRIBUTE, true);
            }

            if (!featureContext.getValues(MTAC.normalPanelFeatureInput, mergedVC.getStart()).isEmpty()) {
                callVcb.attribute(IN_PON_VCF_ATTRIBUTE, true);
            }

            final VariantContext call = addGenotypes(hasNormal, allSomaticAlleles, readAlleleLikelihoods, tumorAlleleFractions, callVcb);
            // how should we be making use of _perSampleFilteredReadList_?
            readAlleleLikelihoods = prepareReadAlleleLikelihoodsForAnnotation(readLikelihoods, perSampleFilteredReadList,
                    false, alleleMapper, readAlleleLikelihoods, call);

            final VariantContext annotatedCall =  annotationEngine.annotateContext(call, featureContext, referenceContext, readAlleleLikelihoods, a -> true);

            call.getAlleles().stream().map(alleleMapper::get).filter(Objects::nonNull).forEach(calledHaplotypes::addAll);
            returnCalls.add( annotatedCall );
        }

        final List<VariantContext> outputCalls = doPhysicalPhasing ? phaseCalls(returnCalls, calledHaplotypes) : returnCalls;
        final int eventCount = outputCalls.size();
        final List<VariantContext> outputCallsWithEventCountAnnotation = outputCalls.stream()
                .map(vc -> new VariantContextBuilder(vc).attribute(GATKVCFConstants.EVENT_COUNT_IN_HAPLOTYPE_KEY, eventCount).make())
                .collect(Collectors.toList());
        return new CalledHaplotypes(outputCallsWithEventCountAnnotation, calledHaplotypes);
    }

    // compute the likelihoods of: AA, AB, AC. . . where A is ref and B, C. . . are somatic alts
    private PerAlleleCollection<Double> getHetGenotypeLogOdds(ReadLikelihoods<Allele> likelihoods,
                                                              final boolean isGermline,
                                                              final String sample) {
        final OptionalDouble givenAltAlleleFraction = isGermline ? GERMLINE_HET_ALT_FRACTION : NO_FIXED_TUMOR_ALT_FRACTION;
        final String sampleForAlleleFractions = (sample.equals(matchedNormalSampleName) && isGermline) ? matchedNormalSampleName : tumorSampleName;
        final PerAlleleCollection<Double> hetGenotypeLogLks = getHetGenotypeLogLikelihoods(likelihoods, sampleForAlleleFractions, sample, givenAltAlleleFraction);
        final PerAlleleCollection<Double> lods = new PerAlleleCollection<>(PerAlleleCollection.Type.ALT_ONLY);
        final int flipLodFactor = isGermline ? -1 : 1;
        final List<Allele> altAlleles = likelihoods.alleles().stream().filter(Allele::isNonReference).collect(Collectors.toList());
        lods.set(altAlleles, a -> flipLodFactor * (hetGenotypeLogLks.get(a) - hetGenotypeLogLks.getRef()));
        return lods;
    }

    private VariantContext addGenotypes(final boolean hasNormal, final List<Allele> allSomaticAlleles, final ReadLikelihoods<Allele> likelihoods,
                                        final PerAlleleCollection<Double> altAlleleFractions, final VariantContextBuilder callVcb) {
        final PerAlleleCollection<Integer> tumorAlleleDepths = getAlleleCounts(likelihoods, tumorSampleName);
        final Genotype tumorGenotype = new GenotypeBuilder(tumorSampleName, allSomaticAlleles)
                .AD(allSomaticAlleles.stream().mapToInt(tumorAlleleDepths::get).toArray())
                .attribute(GATKVCFConstants.ALLELE_FRACTION_KEY, allSomaticAlleles.stream().filter(Allele::isNonReference).mapToDouble(altAlleleFractions::get).toArray())
                .make();

        final List<Genotype> genotypes = new ArrayList<>(Arrays.asList(tumorGenotype));

        // TODO: We shouldn't always assume that the genotype in the normal is hom ref
        final Allele ref = allSomaticAlleles.stream().filter(Allele::isReference).findFirst().get();
        final List<Allele> homRefAllelesforNormalGenotype = Collections.nCopies(2, ref);

        // if we are calling with a normal, build the genotype for the sample to appear in vcf
        if (hasNormal) {
            final PerAlleleCollection<Integer> normalAlleleDepths = getAlleleCounts(likelihoods, matchedNormalSampleName);
            final int normalTotalDepth = likelihoods.sampleReadCount(likelihoods.indexOfSample(matchedNormalSampleName));
            final double[] normalAlleleFractions = allSomaticAlleles.stream().filter(Allele::isNonReference)
                    .mapToDouble(a ->  normalAlleleDepths.get(a) / (double) normalTotalDepth).toArray();

            final Genotype normalGenotype = new GenotypeBuilder(matchedNormalSampleName, homRefAllelesforNormalGenotype)
                    .AD(allSomaticAlleles.stream().mapToInt(normalAlleleDepths::get).toArray())
                    .attribute(GATKVCFConstants.ALLELE_FRACTION_KEY, normalAlleleFractions)
                    .make();
            genotypes.add(normalGenotype);
        }

        return new VariantContextBuilder(callVcb).alleles(allSomaticAlleles).genotypes(genotypes).make();
    }

    // the latent variable z in the strand artifact filter model
    public enum StrandArtifactZ {
        ART_FWD, ART_REV, NO_ARTIFACT;
    }

    // (* model priors *)
    // pseudocounts for the beta distribution over epsilon
    // alpha > 0 and beta > 0. alpha = beta = 1 gives us the flat prior
    final int alpha = 2;
    final int beta = 6; // give more prior weight to beta (i.e. pseudocount for tails) so that the peak shifts towards 0

    // prior probabilities for the latent variable z = {ART_FWD, ART_REV, NO_ART}
    final double[] pi = new double[]{0.01, 0.01, 0.98};

    private void addStrandArtifactAnnotations(final ReadLikelihoods<Allele> likelihoods,
                                              final VariantContextBuilder callVcb,
                                              final List<Allele> allSomaticAlleles) {
        // We simply take the first alternate allele for now
        final Allele altAlelle = allSomaticAlleles.get(1);

        final Collection<ReadLikelihoods<Allele>.BestAllele> tumorBestAlleles = likelihoods.bestAlleles(tumorSampleName);
        final int numAltReadsForward = (int) tumorBestAlleles.stream().filter(ba -> !ba.read.isReverseStrand() && ba.isInformative() && ba.allele.equals(altAlelle)).count();
        final int numAltReadsReverse = (int) tumorBestAlleles.stream().filter(ba -> ba.read.isReverseStrand() && ba.isInformative() && ba.allele.equals(altAlelle)).count();
        final int numReadsForward = (int) tumorBestAlleles.stream().filter(ba -> !ba.read.isReverseStrand() && ba.isInformative()).count();
        final int numReadsReverse = (int) tumorBestAlleles.stream().filter(ba -> ba.read.isReverseStrand() && ba.isInformative()).count();
        final int numReads = numReadsForward + numReadsReverse;

        final double[] unnormalized_posterior_probabilities = new double[StrandArtifactZ.values().length];
        final double[] maximum_a_posteriori_allele_fraction_estimates = new double[StrandArtifactZ.values().length];

        // Compute the posterior probability of ARTIFACT_FWD and ARTIFACT_REV; it's a double integral over the allele fraction f and epsilon
        // Note we must calculate the unnormalized probabilities for NO_ART, ARTIFACT_FWD and ARTIFACT_REV in order to get the normalized probabilities

        // the integrand is a polynomial of degree n, where n is the number of reads at the locus
        // thus to integrate exactly we need (n/2)+1 points
        final int numIntegPointsForAlleleFraction = numReads / 2 + 1;
        final int numIntegPointsForEpsilon = (numReads + alpha + beta - 2) / 2 + 1;

        final GaussIntegratorFactory integratorFactory = new GaussIntegratorFactory();

        final GaussIntegrator gaussIntegratorForAlleleFraction = integratorFactory.legendre(numIntegPointsForAlleleFraction, 0.0, 1.0);
        final List<Double> gaussIntegrationWeightsForAlleleFraction = IntStream.range(0, numIntegPointsForAlleleFraction)
                .mapToDouble(gaussIntegratorForAlleleFraction::getWeight).boxed().collect(Collectors.toList());
        final List<Double> gaussIntegrationAbscissasForAlleleFraction = IntStream.range(0, numIntegPointsForAlleleFraction)
                .mapToDouble(gaussIntegratorForAlleleFraction::getPoint).boxed().collect(Collectors.toList());

        final GaussIntegrator gaussIntegratorForEpsilon = integratorFactory.legendre(numIntegPointsForEpsilon, 0.0, 1.0);
        final List<Double> gaussIntegrationWeightsForEpsilon = IntStream.range(0, numIntegPointsForEpsilon)
                .mapToDouble(gaussIntegratorForEpsilon::getWeight).boxed().collect(Collectors.toList());
        final List<Double> gaussIntegrationAbscissasForEpsilon = IntStream.range(0, numIntegPointsForEpsilon)
                .mapToDouble(gaussIntegratorForEpsilon::getPoint).boxed().collect(Collectors.toList());

        double likelihoodForArtifactFwd = 0.0;
        double likelihoodForArtifactRev = 0.0;
        double[] forwardLikelihoodsOfAlleleFractions = new double[numIntegPointsForAlleleFraction];
        double[] reverseLikelihoodsOfAlleleFractions = new double[numIntegPointsForAlleleFraction];


        for (int i = 0; i < numIntegPointsForAlleleFraction; i++) {
            final double f = gaussIntegrationAbscissasForAlleleFraction.get(i);
            for (int j = 0; j < numIntegPointsForEpsilon; j++) {
                final double epsilon = gaussIntegrationAbscissasForEpsilon.get(j);
                final double integrandArtifactFwd = getIntegrandGivenArtifact(f, epsilon, numReadsForward, numReadsReverse, numAltReadsForward, numAltReadsReverse);
                final double integrandArtifactRev = getIntegrandGivenArtifact(f, epsilon, numReadsReverse, numReadsForward, numAltReadsReverse, numAltReadsForward);
                likelihoodForArtifactFwd += gaussIntegrationWeightsForAlleleFraction.get(i) * gaussIntegrationWeightsForEpsilon.get(j) * integrandArtifactFwd;
                likelihoodForArtifactRev += gaussIntegrationWeightsForAlleleFraction.get(i) * gaussIntegrationWeightsForEpsilon.get(j) * integrandArtifactRev;

                // for a fixed f, integrate over epsilons. this gives us the likelihood p(x^+, x^- | f, z) for a fixed f, which is proportional to
                // the posterior probability of p(f | x^+, x^-, z)
                forwardLikelihoodsOfAlleleFractions[i] += gaussIntegrationWeightsForEpsilon.get(j) * integrandArtifactFwd;
                reverseLikelihoodsOfAlleleFractions[i] += gaussIntegrationWeightsForEpsilon.get(j) * integrandArtifactRev;
            }
        }

        unnormalized_posterior_probabilities[StrandArtifactZ.ART_FWD.ordinal()] = pi[StrandArtifactZ.ART_FWD.ordinal()] * likelihoodForArtifactFwd;
        unnormalized_posterior_probabilities[StrandArtifactZ.ART_REV.ordinal()] = pi[StrandArtifactZ.ART_REV.ordinal()] * likelihoodForArtifactRev;

        maximum_a_posteriori_allele_fraction_estimates[StrandArtifactZ.ART_FWD.ordinal()] = MathUtils.arrayMax(forwardLikelihoodsOfAlleleFractions);
        maximum_a_posteriori_allele_fraction_estimates[StrandArtifactZ.ART_REV.ordinal()] = MathUtils.arrayMax(reverseLikelihoodsOfAlleleFractions);

        // Compute the posterior probability of NO_ARTIFACT; evaluate a single integral over the allele fraction
        // TODO: write a dot product function and add to MathUtils
        final List<Double> integrandsForNoArtifact = gaussIntegrationAbscissasForAlleleFraction.stream()
                .map(f -> getIntegrandGivenNoArtifact(f, numReadsForward, numReadsReverse, numAltReadsForward, numAltReadsReverse))
                .collect(Collectors.toList());
        final double likelihoodForNoArtifact = IntStream.range(0, numIntegPointsForAlleleFraction).
                mapToDouble(i -> gaussIntegrationWeightsForAlleleFraction.get(i) * integrandsForNoArtifact.get(i)).sum();

        unnormalized_posterior_probabilities[StrandArtifactZ.NO_ARTIFACT.ordinal()] = pi[StrandArtifactZ.NO_ARTIFACT.ordinal()] * likelihoodForNoArtifact;
        // if z = NO_ARTIFACT, MAP estimate for f is the sample alt allele fraction
        maximum_a_posteriori_allele_fraction_estimates[StrandArtifactZ.NO_ARTIFACT.ordinal()] = -1.0; // TODO: get the allele fraction from the variant context

        final double[] posterior_probabilities = MathUtils.normalizeFromRealSpace(unnormalized_posterior_probabilities);
        callVcb.attribute(STRAND_ARTIFACT_POSTERIOR_PROBABILITIES_KEY, posterior_probabilities);
        callVcb.attribute(STRAND_ARTIFACT_MAP_ALLELE_FRACTIONS_KEY, maximum_a_posteriori_allele_fraction_estimates);



    }

    private double getIntegrandGivenNoArtifact(final double f, final int n_plus, final int n_minus, final int x_plus, final int x_minus){
        final BinomialDistribution binomFwd = new BinomialDistribution(n_plus, f);
        final BinomialDistribution binomRev = new BinomialDistribution(n_minus, f);

        return binomFwd.probability(x_plus)*binomRev.probability(x_minus);
    }

    private double getIntegrandGivenArtifact(final double f, final double epsilon, final int nWithArtifact, final int nNoArtifact, final int xWithArtifact, final int xNoArtifact){
        final BinomialDistribution binomWithArtifact = new BinomialDistribution(nWithArtifact, f + epsilon * (1-f));
        final BinomialDistribution binomNoArtifact = new BinomialDistribution(nNoArtifact, f);
        final BetaDistribution betaPrior = new BetaDistribution(alpha, beta);

        return betaPrior.density(epsilon) * binomWithArtifact.probability(xWithArtifact) * binomNoArtifact.probability(xNoArtifact);
    }


    /** Calculate the likelihoods of hom ref and each het genotype of the form ref/alt

     * @param givenAltAlleleFraction                     constant fixed alt allele fraction, e.g. f=1/2 for evaluating
     *                                                   germline het likelihoods.  If not given, estimate from the reads
     *
     * @return                                      genotype likelihoods for homRef and het for each alternate allele
     */
    private PerAlleleCollection<Double> getHetGenotypeLogLikelihoods(final ReadLikelihoods<Allele> likelihoods,
                                                                     final String sampleNameForAlleleFractions,
                                                                     final String sampleNameForLikelihoods,
                                                                     final OptionalDouble givenAltAlleleFraction) {
        final Optional<PerAlleleCollection<Double>> alleleFractions = givenAltAlleleFraction.isPresent() ?
                Optional.empty() : Optional.of(getAlleleFractions(likelihoods, sampleNameForAlleleFractions));
        final PerAlleleCollection<MutableDouble> genotypeLogLikelihoods = new PerAlleleCollection<>(PerAlleleCollection.Type.REF_AND_ALT);
        genotypeLogLikelihoods.set(likelihoods.alleles(), a -> new MutableDouble(0));

        final Allele refAllele = genotypeLogLikelihoods.getRefAllele();
        final int refAlleleIndex = likelihoods.indexOfAllele(refAllele);

        // note: indexed by allele (row), read (column)
        final LikelihoodMatrix<Allele> matrix = likelihoods.sampleMatrix(likelihoods.indexOfSample(sampleNameForLikelihoods));

        final int numReads = matrix.numberOfReads();
        for (int alleleIndex = 0; alleleIndex < matrix.numberOfAlleles(); alleleIndex++) {
            if (alleleIndex == refAlleleIndex) {
                continue;
            }
            final Allele altAllele = matrix.getAllele(alleleIndex);
            final double altAlleleFraction = givenAltAlleleFraction.orElseGet(() -> alleleFractions.get().getAlt(altAllele));

            for (int readIndex = 0; readIndex < numReads; readIndex++) {
                final GATKRead read = matrix.getRead(readIndex);
                //TODO: do we need to check for MQ = 0? Haven't such reads already been filtered?
                if (read.getMappingQuality() != 0) {
                    final double readRefLogLikelihood = matrix.get(refAlleleIndex, readIndex);
                    final double readAltLogLikelihood = matrix.get(alleleIndex, readIndex);
                    genotypeLogLikelihoods.get(altAllele).add(hetLog10Likelihood(readRefLogLikelihood, readAltLogLikelihood, altAlleleFraction));
                }
            }
        }

        final PerAlleleCollection<Double> result = new PerAlleleCollection<>(PerAlleleCollection.Type.REF_AND_ALT);
        final double refLogLikelihood = IntStream.range(0, numReads)
                .mapToDouble(i -> matrix.get(refAlleleIndex, i))
                .sum();
        result.setRef(refAllele, refLogLikelihood);
        likelihoods.alleles().stream().filter(Allele::isNonReference).forEach(a -> result.setAlt(a, genotypeLogLikelihoods.get(a).toDouble()));

        return result;
    }

    // get log likelihood for a read to have come from a ref/alt genotype with given alt allele fraction
    private static double hetLog10Likelihood(final double refLogLikelihood, final double altLogLikelihood, final double altAlleleFraction) {
        return MathUtils.log10SumLog10(refLogLikelihood + Math.log10(1 - altAlleleFraction), altLogLikelihood + Math.log10(altAlleleFraction));
    }

    // TODO: calculate using the uncertainty rather than this cheap approach
    private PerAlleleCollection<Double> getAlleleFractions(final ReadLikelihoods<Allele> likelihoods,
                                                           final String sampleName) {
        final List<Allele> alleles = likelihoods.alleles();
        final PerAlleleCollection<Integer> alleleCounts = getAlleleCounts(likelihoods, sampleName);
        final PerAlleleCollection<Double> alleleFractions = new PerAlleleCollection<>(PerAlleleCollection.Type.REF_AND_ALT);
        final int totalCount = alleleCounts.getRef() + alleleCounts.getAltAlleles().stream().mapToInt(a -> alleleCounts.getAlt(a)).sum();
        alleleFractions.set(alleles, a -> totalCount == 0 ? 0 : alleleCounts.get(a) / (double) totalCount);

        return alleleFractions;
    }


    private PerAlleleCollection<Integer> getAlleleCounts(final ReadLikelihoods<Allele> likelihoods,
                                                         final String sampleName) {
        final List<Allele> alleles = likelihoods.alleles();

        final PerAlleleCollection<MutableInt> alleleCounts = new PerAlleleCollection<>(PerAlleleCollection.Type.REF_AND_ALT);
        likelihoods.alleles().stream().forEach(a -> alleleCounts.set(a, new MutableInt(0)));

        for (final ReadLikelihoods<Allele>.BestAllele bestAllele : likelihoods.bestAlleles(sampleName)) {
            final GATKRead read = bestAllele.read;
            if (read.getMappingQuality() > 0 && bestAllele.isInformative()) {
                alleleCounts.get(bestAllele.allele).increment();
            }
        }

        final PerAlleleCollection<Integer> result = new PerAlleleCollection<>(PerAlleleCollection.Type.REF_AND_ALT);
        alleles.stream().forEach(a -> result.set(a, alleleCounts.get(a).toInteger()));

        return result;
    }

    private void filterOverlappingReads(final ReadLikelihoods<Allele> likelihoods, final String sample, final Allele ref, final int location, final boolean retainMismatches) {
        Utils.validateArg(likelihoods.indexOfSample(sample) >= 0, "Sample is missing from likelihoods");

        // Get the best alleles of each read and group them by the read name.
        // This puts paired reads from the same fragment together
        final Map<String, List<ReadLikelihoods<Allele>.BestAllele>> fragments = likelihoods.bestAlleles(sample).stream()
                .collect(Collectors.groupingBy(ba -> ba.read.getName()));

        // We only potentially filter read pairs that overlap at this position
        final List<Pair<ReadLikelihoods<Allele>.BestAllele, ReadLikelihoods<Allele>.BestAllele>> overlappingReadPairs =
                fragments.values().stream()
                        .filter(l -> l.size() == 2)
                        .map(l -> new ImmutablePair<>(l.get(0), l.get(1)))
                        .filter(p -> ReadUtils.isInsideRead(p.getLeft().read, location) && ReadUtils.isInsideRead(p.getRight().read, location))
                        .collect(Collectors.toList());

        final Set<GATKRead> readsToDiscard = new HashSet<>();

        for (final Pair<ReadLikelihoods<Allele>.BestAllele, ReadLikelihoods<Allele>.BestAllele> pair : overlappingReadPairs) {
            final ReadLikelihoods<Allele>.BestAllele read = pair.getLeft();
            final ReadLikelihoods<Allele>.BestAllele mate = pair.getRight();

            if (read.allele.equals(mate.allele)) {
                // keep the higher-quality read
                readsToDiscard.add(read.likelihood < mate.likelihood ? read.read : mate.read);
            } else if (retainMismatches) {
                // keep the alt read
                readsToDiscard.add(read.allele.equals(ref) ? read.read : mate.read);
            } else {
                // throw out both
                readsToDiscard.add(read.read);
                readsToDiscard.add(mate.read);
            }
        }

        likelihoods.removeSampleReads(likelihoods.indexOfSample(sample), readsToDiscard, likelihoods.numberOfAlleles());
    }

}
