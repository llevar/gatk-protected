package org.broadinstitute.hellbender.tools.coveragemodel;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import javax.annotation.Nonnull;
import java.io.Serializable;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * This class represents a key for a contiguous block of elements in an indexed linear space.
 *
 * The begin index is inclusive
 * The end index is exclusive
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
final class LinearlySpacedIndexBlock implements Serializable {

    private static final long serialVersionUID = 4571138983754570342L;

    /* begin index inclusive, end index exclusive */
    private final int begIndex, endIndex, numElements;

    public LinearlySpacedIndexBlock(final int begIndex, final int endIndex) {
        this.begIndex = ParamUtils.isPositiveOrZero(begIndex, "The begin index of a block must be non-negative.");
        this.endIndex = ParamUtils.inRange(endIndex, begIndex + 1, Integer.MAX_VALUE, "The block must at least" +
                " contain one element.");
        numElements = endIndex - begIndex;
    }

    public int getBegIndex() { return begIndex; }

    public int getEndIndex() { return endIndex; }

    public int getNumElements() { return numElements; }

    /**
     * Asserts that a collection of linearly spaced index blocks are non-overlapping and fully-covering (i.e. there is no
     * gap between them)
     *
     * @param blocks a collection of linear space blocks
     * @throws IllegalArgumentException if the collection is null
     * @throws AssertionError if blocks overlap or there is a gap between them
     */
    public static void assertNonOverlappingFullyCovering(@Nonnull final Collection<LinearlySpacedIndexBlock> blocks) {
        Utils.nonNull(blocks, "The collection of linear space blocks must be non-null");
        final List<LinearlySpacedIndexBlock> sortedBlocks = blocks.stream()
                .sorted(Comparator.comparingInt(LinearlySpacedIndexBlock::getBegIndex))
                .collect(Collectors.toList());
        if (!IntStream.range(0, sortedBlocks.size() - 1)
                .allMatch(idx -> sortedBlocks.get(idx).getEndIndex() == sortedBlocks.get(idx + 1).getBegIndex())) {
            throw new AssertionError("Some of the blocks in the collection either overlap or have gaps between them");
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (!(o instanceof LinearlySpacedIndexBlock)) {
            return false;
        }

        final LinearlySpacedIndexBlock block = (LinearlySpacedIndexBlock) o;
        return (begIndex == block.begIndex) && (endIndex == block.endIndex);
    }

    /**
     * The best hash code is the {@link LinearlySpacedIndexBlock#begIndex} for non-overlapping and fully-covering
     * blocks
     *
     * @return hash code
     */
    @Override
    public int hashCode() {
        return begIndex;
    }

    @Override
    public String toString() {
        return "[" + begIndex + ", " + endIndex + "]";
    }
}
