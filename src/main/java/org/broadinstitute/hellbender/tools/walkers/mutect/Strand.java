package org.broadinstitute.hellbender.tools.walkers.mutect;

import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Created by davidben on 3/14/17.
 */
public enum Strand {
    FORWARD, REVERSE, BOTH;

    public static boolean readComesFromStrand(final GATKRead read, final Strand strand) {
        return strand == BOTH || (read.isReverseStrand() ? strand == REVERSE : strand == FORWARD);
    }
}
