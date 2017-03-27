package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.apache.commons.math3.exception.MathIllegalArgumentException;
import org.apache.commons.math3.stat.descriptive.rank.Median;
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

    public static final String MEDIAN_ALT_POSITION_NAME = "MED_ALT_POS";
    public static final String MEDIAN_ALT_CLIPPING_NAME = "MED_ALT_CLIP";
    public static final String MEDIAN_ALT_COMPLEXITY_NAME = "MED_ALT_CIGAR";


    @Override
    public void onTraversalStart() {
        final VCFHeader inputHeader = getHeaderForVariants();
        final Set<VCFHeaderLine> headerLines = new HashSet<>(inputHeader.getMetaDataInSortedOrder());
        headerLines.add(new VCFInfoHeaderLine(MEDIAN_ALT_POSITION_NAME, 1, VCFHeaderLineType.Float, "median alt distance from end of read"));
        headerLines.add(new VCFInfoHeaderLine(MEDIAN_ALT_CLIPPING_NAME, 1, VCFHeaderLineType.Float, "median clipped bases in alt reads"));
        headerLines.add(new VCFInfoHeaderLine(MEDIAN_ALT_COMPLEXITY_NAME, 1, VCFHeaderLineType.Float, "median number of cigar elements in alt reads"));

        final VCFHeader vcfHeader = new VCFHeader(headerLines, inputHeader.getGenotypeSamples());
        vcfWriter = createVCFWriter(new File(outputVcf));
        vcfWriter.writeHeader(vcfHeader);
    }

    @Override
    public void apply(final VariantContext vc, final ReadsContext readsContext, final ReferenceContext refContext, final FeatureContext fc) {
        final ReadPileup pileup = GATKProtectedVariantContextUtils.getPileup(vc, readsContext);
        final byte refBase = vc.getReference().getBases()[0];
        if (vc.isSymbolicOrSV()) {
            return;
        }

        try {
            final double medianAltPosition = new Median().evaluate(Utils.stream(pileup)
                    .filter(pe -> pe.getBase() != refBase)
                    .mapToDouble(pe -> Math.min(pe.getRead().getLength() - pe.getOffset(), pe.getOffset())).toArray());

            final double medianAltClipping = new Median().evaluate(Utils.stream(pileup)
                    .filter(pe -> pe.getBase() != refBase)
                    .mapToDouble(pe -> AlignmentUtils.getNumHardClippedBases(pe.getRead())).toArray());

            final double medianAltCigar = new Median().evaluate(Utils.stream(pileup)
                    .filter(pe -> pe.getBase() != refBase)
                    .mapToDouble(pe -> pe.getRead().numCigarElements())
                    .toArray());

            final VariantContext annotatedVariant = new VariantContextBuilder(vc).attribute(MEDIAN_ALT_POSITION_NAME, medianAltPosition)
                    .attribute(MEDIAN_ALT_CLIPPING_NAME, medianAltClipping)
                    .attribute(MEDIAN_ALT_COMPLEXITY_NAME, medianAltCigar).make();

            vcfWriter.add(annotatedVariant);
        } catch (MathIllegalArgumentException ex) {
            vcfWriter.add(vc);
        }

    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }
}
