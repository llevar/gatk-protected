package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.contamination.ContaminationRecord;

import java.io.File;

public class M2FiltersArgumentCollection extends AssemblyBasedCallerArgumentCollection {
    private static final long serialVersionUID = 9345L;

    /**
     * This is a measure of the minimum evidence to support that a variant observed in the tumor is not also present in the normal
     * as an artifact i.e. not as a germline event.
     */
    @Argument(fullName = "normal_artifact_lod", optional = true, doc = "LOD threshold for calling normal artifacts")
    public double NORMAL_ARTIFACT_LOD_THRESHOLD = 0.0;

    /**
     * The LOD threshold for the normal is typically made more strict if the variant has been seen in dbSNP (i.e. another
     * normal sample). We thus require MORE evidence that a variant is NOT seen in this tumor's normal if it has been observed as a germline variant before.
     */
    @Argument(fullName = "dbsnp_normal_lod", optional = true, doc = "LOD threshold for calling normal non-variant at dbsnp sites")
    public double NORMAL_DBSNP_LOD_THRESHOLD = 5.5;

    @Argument(fullName = "maxAltAllelesThreshold", optional = true, doc="filter variants with too many alt alleles")
    public int numAltAllelesThreshold = 1;

    @Argument(fullName = "maxEventsInHaplotype", optional = true, doc="Variants coming from a haplotype with more than this many events are filtered")
    public int maxEventsInHaplotype = 2;

    // TODO: add short names, make them constants
    @Argument(fullName = "strand_artifact_threshold", optional = true, doc = "Filter a variant if the probability of strand artifact exceeds this number")
    public double STRAND_ARTIFACT_POSTERIOR_PROB_THRESHOLD = 1.0;

    @Argument(shortName = "contaminationTable", fullName = "contaminationTable", optional = true, doc="Table containing contamination information.")
    public File contaminationTable = null;

}
