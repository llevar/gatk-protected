package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Created by David Benjamin 7/15/15
 */
public final class SegmentUtils {
    /**
     * read a list of segments without calls from a file
     */
    public static List<Segment> readUncalledSegments(final File segmentsFile) throws IOException {
        try (final TableReader<Segment> reader = TableUtils.reader(segmentsFile,
                (columns, formatExceptionFactory) -> {
                    if (!columns.matchesAll(0, "Sample", "Chromosome", "Start", "End")) {//ignore last two columns: NTARGETS and AVERAGE
                        throw formatExceptionFactory.apply("Bad header");
                    }

                    // return the lambda to translate dataLines into uncalled segments.
                    return (dataLine) -> new Segment(dataLine.get(0),
                            new SimpleInterval(dataLine.get(1), dataLine.getInt(2), dataLine.getInt(3)));
                })) {
            return reader.stream().collect(Collectors.toList());
        } catch (final UncheckedIOException e) {
            throw e.getCause();
        }
    }

    /**
     * read a list of segments with calls from a file
     */
    public static List<Segment> readCalledSegments(final File segmentsFile) throws IOException {
        try (final TableReader<Segment> reader = TableUtils.reader(segmentsFile,
                (columns, formatExceptionFactory) -> {
                    if (!columns.matchesExactly("Sample", "Chromosome", "Start", "End", "Call")) {
                        throw formatExceptionFactory.apply("Bad header");
                    }
                    // return the lambda to translate dataLines into called segments.
                    return (dataLine) -> new Segment(dataLine.get(0),
                            new SimpleInterval(dataLine.get(1), dataLine.getInt(2), dataLine.getInt(3)), dataLine.get(4));
                })) {
            return reader.stream().collect(Collectors.toList());
        } catch (final UncheckedIOException e) {
            throw e.getCause();
        }
    }

    /**
     * write a list of segments with calls to file
     */
    public static void writeCalledSegments(final File outFile, List<Segment> segments) throws IOException {
        try (final TableWriter<Segment> writer = TableUtils.writer(outFile,
                    new TableColumnCollection("Sample", "Chromosome", "Start", "End", "Call"),

                    //lambda for filling an initially empty DataLine
                    (segment, dataLine) -> {
                        SimpleInterval interval = segment.getInterval();
                        dataLine.append(segment.getSample(), interval.getContig()).append(interval.getStart(), interval.getEnd())
                                .append(segment.getCall());
                    })) {
            for (final Segment segment : segments) {
                if (segment.getCall() == null) {
                    throw new IOException("Segment at " + segment.getInterval().toString() + " has no call.");
                }
                writer.writeRecord(segment);
            }
        }
    }

    /**
     * Returns the overlapping targets.
     *
     * Delegates to ExonCollection's binary search.
     */
    public static List<TargetCoverage> overlappingTargets(Segment segment, final ExonCollection<TargetCoverage> targets) {
        Utils.nonNull(segment, "Can't get overlapping targets of null segment.");
        Utils.nonNull(targets, "Overlapping targets requires non-null targets collection.");
        return targets.exons(segment.getInterval());
    }

    /**
     * the mean of all overlapping targets' coverages
     *
     * @throws IllegalStateException if overlapping targets have not been assigned or if no overlapping targets were found.
     */
    public static double meanTargetCoverage(final Segment segment, final ExonCollection<TargetCoverage> targets) {
        Utils.nonNull(segment, "Can't get mean coverage of null segment.");
        Utils.nonNull(targets, "Mean target coverage requires non-null targets collection.");

        final List<TargetCoverage> myTargets = overlappingTargets(segment, targets);

        if (myTargets.size() == 0) {
            throw new IllegalStateException("Empty segment -- no overlapping targets.");
        }
        return myTargets.stream().mapToDouble(TargetCoverage::getCoverage).average().getAsDouble();
    }

    public static int numTargets(Segment segment, final ExonCollection<TargetCoverage> targets) {
        return overlappingTargets(segment, targets).size();
    }
}