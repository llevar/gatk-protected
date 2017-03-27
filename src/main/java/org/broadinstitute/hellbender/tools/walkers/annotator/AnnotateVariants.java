package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.DbsnpArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyResultSet;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.LikelihoodEngineArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReadLikelihoodCalculationEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.File;
import java.util.*;

/**
 * Created by David Benjamin on 3/24/17.
 */
@CommandLineProgramProperties(
        summary = "Annotate vcf",
        oneLineSummary = "Annotate vcf",
        programGroup = VariantProgramGroup.class
)
public class AnnotateVariants extends VariantWalker {

    @Argument(fullName= StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The output filtered VCF file",
            optional=false)
    private final String outputVcf = null;

    @ArgumentCollection
    public LikelihoodEngineArgumentCollection likelihoodArgs = new LikelihoodEngineArgumentCollection();

    @Argument(fullName = "group", shortName = "G", doc = "One or more classes/groups of annotations to apply to variant calls", optional = true)
    public List<String> annotationGroupsToUse = new ArrayList<>(Arrays.asList(new String[]{ }));

    @Argument(fullName="annotation", shortName="A", doc="One or more specific annotations to apply to variant calls", optional = true)
    protected List<String> annotationsToUse = new ArrayList<>(Arrays.asList(new String[]{}));

    @Argument(fullName = "excludeAnnotation", shortName = "XA", doc = "One or more specific annotations to exclude", optional = true)
    public List<String> annotationsToExclude = new ArrayList<>();

    @ArgumentCollection
    public DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    @Argument(fullName = "comp", shortName = "comp", doc = "Comparison VCF file(s)", optional = true)
    public List<FeatureInput<VariantContext>> comps = Collections.emptyList();

    private VariantContextWriter vcfWriter;
    private  SAMFileHeader bamHeader;
    private org.broadinstitute.hellbender.utils.genotyper.SampleList sampleList;
    private ReadLikelihoodCalculationEngine likelihoodCalculationEngine;
    private VariantAnnotatorEngine annotationEngine;

    @Override
    public void onTraversalStart() {
        final VCFHeader inputHeader = getHeaderForVariants();
        final Set<VCFHeaderLine> headerLines = new HashSet<>(inputHeader.getMetaDataInSortedOrder());
        final VCFHeader originalVcfHEader = new VCFHeader(headerLines, inputHeader.getGenotypeSamples());
        vcfWriter = createVCFWriter(new File(outputVcf));

        bamHeader = getHeaderForReads();
        sampleList = new IndexedSampleList(new ArrayList<>(ReadUtils.getSamplesFromHeader(bamHeader)));
        likelihoodCalculationEngine = AssemblyBasedCallerUtils.createLikelihoodCalculationEngine(likelihoodArgs);

        annotationEngine = VariantAnnotatorEngine.ofSelectedMinusExcluded(annotationGroupsToUse, annotationsToUse, annotationsToExclude, dbsnp.dbsnp, comps);

        final Set<VCFHeaderLine> headerInfo = new HashSet<>(originalVcfHEader.getMetaDataInInputOrder());
        headerInfo.addAll(annotationEngine.getVCFAnnotationDescriptions());

        final VCFHeader newVcfHeader = new VCFHeader(headerInfo, sampleList.asListOfSamples());
        newVcfHeader.setSequenceDictionary(getBestAvailableSequenceDictionary());
        vcfWriter.writeHeader(newVcfHeader);
    }

    @Override
    public void apply(final VariantContext vc, final ReadsContext readsContext, final ReferenceContext refContext, final FeatureContext fc) {
        final SimpleInterval refWindow = refContext.getWindow();

        final Haplotype refHaplotype = createReferenceHaplotype(refContext);
        final List<Haplotype> variantHaplotypes = ReadThreadingAssembler.composeGivenHaplotypes(refHaplotype, Arrays.asList(vc), refWindow);

        final AssemblyRegion assemblyRegion = new AssemblyRegion(refWindow, 0, bamHeader);

        // NOTE: we don't actually perform assembly here!
        final AssemblyResultSet assemblyResultSet = new AssemblyResultSet();
        assemblyResultSet.add(refHaplotype);
        variantHaplotypes.forEach(assemblyResultSet::add);
        assemblyResultSet.setRegionForGenotyping(assemblyRegion);

        final Map<String,List<GATKRead>> reads = AssemblyBasedCallerUtils.splitReadsBySample(sampleList, bamHeader, assemblyRegion.getReads());

        final ReadLikelihoods<Haplotype> readLikelihoods = likelihoodCalculationEngine.computeReadLikelihoods(assemblyResultSet, sampleList, reads);
        final Map<GATKRead,GATKRead> readRealignments = AssemblyBasedCallerUtils.realignReadsToTheirBestHaplotype(readLikelihoods, assemblyResultSet.getReferenceHaplotype(), assemblyResultSet.getPaddedReferenceLoc());
        readLikelihoods.changeReads(readRealignments);

        final Map<Allele, List<Haplotype>> alleleMapper = new HashMap<>();
        alleleMapper.put(vc.getReference(), Arrays.asList(refHaplotype));
        for (int n = 0; n < vc.getAlternateAlleles().size(); n++) {
            alleleMapper.put(vc.getAlternateAllele(n), Arrays.asList(variantHaplotypes.get(n)));
        }

        ReadLikelihoods<Allele> readAlleleLikelihoods = readLikelihoods.marginalize(alleleMapper, vc);
        final VariantContext annotatedCall =  annotationEngine.annotateContext(vc, fc, refContext, readAlleleLikelihoods, a -> true);
        vcfWriter.add(annotatedCall);

    }

    private static Haplotype createReferenceHaplotype(ReferenceContext refContext) {
        final byte[] refBases = refContext.getBases();
        final SimpleInterval refWindow = refContext.getWindow();
        final Haplotype refHaplotype = new Haplotype(refBases, true);
        refHaplotype.setAlignmentStartHapwrtRef(refWindow.getStart());
        final Cigar c = new Cigar();
        c.add(new CigarElement(refHaplotype.getBases().length, CigarOperator.M));
        refHaplotype.setCigar(c);
        return refHaplotype;
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }
}
