package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.samtools.CigarOperator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodsCalculationModel;
import org.broadinstitute.hellbender.tools.walkers.genotyper.VariantCallContext;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.activityprofile.ActivityProfileState;
import org.broadinstitute.hellbender.utils.downsampling.DownsamplingMethod;
import org.broadinstitute.hellbender.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.iterators.ReadFilteringIterator;
import org.broadinstitute.hellbender.utils.locusiterator.AlignmentStateMachine;
import org.broadinstitute.hellbender.utils.locusiterator.LocusIteratorByState;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.utils.samples.Sample;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class HaplotypeCallerEngineUnitTest extends BaseTest {
    private static final Allele FAKE_REF_ALLELE = Allele.create("N", true); // used in isActive function to call into UG Engine. Should never appear anywhere in a VCF file
    private static final Allele FAKE_ALT_ALLELE = Allele.create("<FAKE_ALT>", false); // used in isActive function to call into UG Engine. Should never appear anywhere in a VCF file
    private static final double AVERAGE_HQ_SOFTCLIPS_HQ_BASES_THRESHOLD = 6.0;

    private static final byte REF_MODEL_DELETION_QUAL = 30;

    @Test
    public void testIsActive() throws IOException {
        final File testBam = new File(NA12878_20_21_WGS_bam);
        final File reference = new File(b37_reference_20_21);
        final SimpleInterval shardInterval = new SimpleInterval("20", 10000000, 10001000);
        final SimpleInterval paddedShardInterval = new SimpleInterval(shardInterval.getContig(), shardInterval.getStart() - 100, shardInterval.getEnd() + 100);
        final HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();

        // We expect isActive() to return 1.0 for the sites below, and 0.0 for all other sites
        final List<SimpleInterval> expectedActiveSites = Arrays.asList(
                new SimpleInterval("20", 9999996, 9999996),
                new SimpleInterval("20", 9999997, 9999997),
                new SimpleInterval("20", 10000117, 10000117),
                new SimpleInterval("20", 10000211, 10000211),
                new SimpleInterval("20", 10000439, 10000439),
                new SimpleInterval("20", 10000598, 10000598),
                new SimpleInterval("20", 10000694, 10000694),
                new SimpleInterval("20", 10000758, 10000758),
                new SimpleInterval("20", 10001019, 10001019)
        );

        try ( final ReadsDataSource reads = new ReadsDataSource(testBam.toPath());
              final ReferenceDataSource ref = new ReferenceFileSource(reference);
              final CachingIndexedFastaSequenceFile referenceReader = new CachingIndexedFastaSequenceFile(reference);) {

            final HaplotypeCallerEngine hcEngine = new HaplotypeCallerEngine(hcArgs, reads.getHeader(), referenceReader);

            List<ReadFilter> hcFilters = HaplotypeCallerEngine.makeStandardHCReadFilters();
            hcFilters.forEach(filter -> filter.setHeader(reads.getHeader()));
            ReadFilter hcCombinedFilter = hcFilters.get(0);
            for ( int i = 1; i < hcFilters.size(); ++i ) {
                hcCombinedFilter = hcCombinedFilter.and(hcFilters.get(i));
            }
            final Iterator<GATKRead> readIter = new ReadFilteringIterator(reads.query(paddedShardInterval), hcCombinedFilter);

            final LocusIteratorByState libs = new LocusIteratorByState(readIter, DownsamplingMethod.NONE, false, ReadUtils.getSamplesFromHeader(reads.getHeader()), reads.getHeader(), false);

            for ( final AlignmentContext pileup : libs ) {
                final SimpleInterval pileupInterval = new SimpleInterval(pileup.getLocation());
                final ReferenceContext pileupRefContext = new ReferenceContext(ref, pileupInterval);

                final ActivityProfileState isActiveResult = hcEngine.isActive(pileup, pileupRefContext, new FeatureContext(null, pileupInterval));

                final double expectedIsActiveValue = expectedActiveSites.contains(pileupInterval) ? 1.0 : 0.0;
                Assert.assertEquals(isActiveResult.isActiveProb(), expectedIsActiveValue, "Wrong isActive probability for site " + pileupInterval);
            }
        }
    }
    @Test
    public void testBlabla() throws IOException {
        final File testBam = new File(NA12878_20_21_WGS_bam);
        final File reference = new File(b37_reference_20_21);
        final SimpleInterval shardInterval = new SimpleInterval("20", 10000000, 10001000);
        final SimpleInterval paddedShardInterval = new SimpleInterval(shardInterval.getContig(), shardInterval.getStart() - 100, shardInterval.getEnd() + 100);
        final HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();

        // We expect isActive() to return 1.0 for the sites below, and 0.0 for all other sites
        final List<SimpleInterval> expectedActiveSites = Arrays.asList(
                new SimpleInterval("20", 9999996, 9999996),
                new SimpleInterval("20", 9999997, 9999997),
                new SimpleInterval("20", 10000117, 10000117),
                new SimpleInterval("20", 10000211, 10000211),
                new SimpleInterval("20", 10000439, 10000439),
                new SimpleInterval("20", 10000598, 10000598),
                new SimpleInterval("20", 10000694, 10000694),
                new SimpleInterval("20", 10000758, 10000758),
                new SimpleInterval("20", 10001019, 10001019)
        );

        Map<SimpleInterval, ActivityProfileState> activityStates = new HashMap();
        Map<SimpleInterval, RefVsAnyResult> refVsAnyResults = new HashMap();


        try ( final ReadsDataSource reads = new ReadsDataSource(testBam.toPath());
              final ReferenceDataSource ref = new ReferenceFileSource(reference);
              final CachingIndexedFastaSequenceFile referenceReader = new CachingIndexedFastaSequenceFile(reference);) {

            final HaplotypeCallerEngine hcEngine = new HaplotypeCallerEngine(hcArgs, reads.getHeader(), referenceReader);

            List<ReadFilter> hcFilters = HaplotypeCallerEngine.makeStandardHCReadFilters();
            hcFilters.forEach(filter -> filter.setHeader(reads.getHeader()));
            ReadFilter hcCombinedFilter = hcFilters.get(0);
            for ( int i = 1; i < hcFilters.size(); ++i ) {
                hcCombinedFilter = hcCombinedFilter.and(hcFilters.get(i));
            }
            final Iterator<GATKRead> readIter = new ReadFilteringIterator(reads.query(paddedShardInterval), hcCombinedFilter);

            List<GATKRead> readList = new ArrayList<>();
            readIter.forEachRemaining(readList::add);
            Collections.shuffle(readList, new Random(1));
            System.out.println(readList.size());

            SampleList samples = hcEngine.getSamplesList();
            String thisSample = samples.getSample(0);

            final int ploidy = hcEngine.getActiveRegionEvaluationGenotyperEngine().getConfiguration().genotypeArgs.samplePloidy;
            final List<Allele> noCall = GATKVariantContextUtils.noCallAlleles(ploidy); // used to noCall all genotypes until the exact model is applied

            for(GATKRead read : readList){
                AlignmentStateMachine readASM = new AlignmentStateMachine(read);
                CigarOperator currentCigarOp = readASM.stepForwardOnGenome();
                System.out.println("Processing read " + read);
                while(currentCigarOp != null){
                    if(!readASM.isLeftEdge()) {
                        SimpleInterval currentLocation = readASM.getLocation();
                        ReferenceContext rc = new ReferenceContext(ref, currentLocation);

                        final GenotypesContext genotypes = GenotypesContext.create(1);
                        final MathUtils.RunningAverage averageHQSoftClips = new MathUtils.RunningAverage();

                        List<PileupElement> els = new ArrayList<PileupElement>();
                        PileupElement p = readASM.makePileupElement();
                        els.add(p);
                        ReadPileup readPileup = new ReadPileup(currentLocation, els);

                        RefVsAnyResult currentResult = refVsAnyResults.get(currentLocation);

                        if(currentResult != null){
                            final int likelihoodCount = ploidy + 1;
                            final double log10Ploidy = MathUtils.log10(ploidy);

                            final byte qual = p.isDeletion() ? REF_MODEL_DELETION_QUAL : p.getQual();
                            if (!p.isDeletion() && qual <= hcArgs.minBaseQualityScore){
                                currentCigarOp = readASM.stepForwardOnGenome();
                                continue;
                            }


                            hcEngine.
                                    getReferenceConfidenceModel().
                                    applyPileupElementRefVsNonRefLikelihoodAndCount(rc.getBase(),
                                        likelihoodCount,
                                        log10Ploidy,
                                        currentResult,
                                        p,
                                        qual,
                                        averageHQSoftClips);

                            for (int i = 0; i < likelihoodCount; i++) {
                                currentResult.addGenotypeLikelihood(i, -log10Ploidy);
                            }

                        }else{
                            currentResult = hcEngine.
                                    getReferenceConfidenceModel().
                                    calcGenotypeLikelihoodsOfRefVsAny(ploidy,
                                            readPileup,
                                            rc.getBase(),
                                            hcArgs.minBaseQualityScore,
                                            averageHQSoftClips);
                            refVsAnyResults.put(currentLocation, currentResult);
                        }

                        final double[] genotypeLikelihoods = currentResult.getGenotypeLikelihoods();

                        genotypes.add(new GenotypeBuilder(thisSample).alleles(noCall).PL(genotypeLikelihoods).make());


                        final List<Allele> alleles = Arrays.asList(FAKE_REF_ALLELE, FAKE_ALT_ALLELE);
                        final double isActiveProb;

                        isActiveProb = hcEngine.getActiveRegionEvaluationGenotyperEngine().calculateSingleSampleRefVsAnyActiveStateProfileValue(genotypes.get(0).getLikelihoods().getAsVector());

                        ActivityProfileState newState = new ActivityProfileState(currentLocation, isActiveProb, averageHQSoftClips.mean() > AVERAGE_HQ_SOFTCLIPS_HQ_BASES_THRESHOLD ? ActivityProfileState.Type.HIGH_QUALITY_SOFT_CLIPS : ActivityProfileState.Type.NONE, averageHQSoftClips.mean());
                        activityStates.put(currentLocation, newState);
                    }

                    currentCigarOp = readASM.stepForwardOnGenome();

                }

            }

        }

        List<ActivityProfileState> results = activityStates.values().stream().filter(activityProfileState -> activityProfileState.isActiveProb() >0 ).collect(Collectors.toList());
    }


}
