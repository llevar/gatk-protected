package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

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


    private VariantContextWriter vcfWriter;
    private  SAMFileHeader bamHeader;
    private org.broadinstitute.hellbender.utils.genotyper.SampleList sampleList;

    @Override
    public void onTraversalStart() {
        final VCFHeader inputHeader = getHeaderForVariants();
        final Set<VCFHeaderLine> headerLines = new HashSet<>(inputHeader.getMetaDataInSortedOrder());
        final VCFHeader originalVcfHEader = new VCFHeader(headerLines, inputHeader.getGenotypeSamples());
        vcfWriter = createVCFWriter(new File(outputVcf));

        bamHeader = getHeaderForReads();
        sampleList = new IndexedSampleList(new ArrayList<>(ReadUtils.getSamplesFromHeader(bamHeader)));


        final Set<VCFHeaderLine> headerInfo = new HashSet<>(originalVcfHEader.getMetaDataInInputOrder());
        headerInfo.addAll(annotationEngine.getVCFAnnotationDescriptions());

        final VCFHeader newVcfHeader = new VCFHeader(headerInfo, sampleList.asListOfSamples());
        newVcfHeader.setSequenceDictionary(getBestAvailableSequenceDictionary());
        vcfWriter.writeHeader(newVcfHeader);
    }

    @Override
    public void apply(final VariantContext vc, final ReadsContext readsContext, final ReferenceContext refContext, final FeatureContext fc) {
        final ReadPileup pileup = GATKProtectedVariantContextUtils.getPileup(vc, readsContext);
        final byte refBase = vc.getReference().getBases()[0];

        final double averageAltPosition = Utils.stream(pileup)
                .filter(pe -> pe.getBase() != refBase)
                .mapToDouble(pe -> Math.min(pe.getRead().getLength() - pe.getOffset(), pe.getOffset()))
                .average().orElseGet(() -> 0.0);

        final double averageAltClipping = Utils.stream(pileup)
                .filter(pe -> pe.getBase() != refBase)
                .mapToDouble(pe -> AlignmentUtils.getNumHardClippedBases(pe.getRead()))
                .average().orElseGet(() -> 0.0);

        final double averageAltCigar = Utils.stream(pileup)
                .filter(pe -> pe.getBase() != refBase)
                .mapToDouble(pe -> pe.getRead().numCigarElements())
                .average().orElseGet(() -> 0.0);


    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }
}
