package picard.sam.markduplicates;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.util.SortingLongCollection;
import picard.sam.markduplicates.util.ReadEndsForMarkDuplicates;
import picard.sam.markduplicates.util.RepresentativeReadIndexerCodec;
import picard.sam.util.RepresentativeReadIndexer;

import java.io.File;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/**
 * MarkDuplicates calculation helper class for flow based mode
 */
public class MarkDuplicatesForFlowHelper implements MarkDuplicatesHelper {

    private final Log log = Log.getInstance(MarkDuplicatesForFlowHelper.class);

    private static  final int END_INSIGNIFICANT_VALUE = 0;
    private static final String ATTR_DUPLICATE_SCORE = "ForFlowDuplicateScore";

    // constants for clippingTagContains
    public static String        CLIPPING_TAG_NAME = "tm";
    public static final char[]  CLIPPING_TAG_CONTAINS_A = {'A'};
    public static final char[]  CLIPPING_TAG_CONTAINS_AQ = {'A', 'Q'};
    public static final char[]  CLIPPING_TAG_CONTAINS_QZ = {'Q', 'Z'};

    // instance of hosting MarkDuplicates
    private MarkDuplicates  md;

    public MarkDuplicatesForFlowHelper(MarkDuplicates md) {
        this.md = md;
    }

    /**
     * This method is identical in function to generateDuplicateIndexes except that it accomodates for
     * the possible significance of the end side of the reads (w/ or wo/ uncertainty). This is only
     * applicable for flow mode invocation.
     */
    public void generateDuplicateIndexes(final boolean useBarcodes, final boolean indexOpticalDuplicates) {
        final int entryOverhead;
        if (md.TAG_DUPLICATE_SET_MEMBERS) {
            // Memory requirements for RepresentativeReadIndexer:
            // three int entries + overhead: (3 * 4) + 4 = 16 bytes
            entryOverhead = 16;
        } else {
            entryOverhead = SortingLongCollection.SIZEOF;
        }
        // Keep this number from getting too large even if there is a huge heap.
        int maxInMemory = (int) Math.min((Runtime.getRuntime().maxMemory() * 0.25) / entryOverhead, (double) (Integer.MAX_VALUE - 5));
        // If we're also tracking optical duplicates, reduce maxInMemory, since we'll need two sorting collections
        if (indexOpticalDuplicates) {
            maxInMemory /= ((entryOverhead + SortingLongCollection.SIZEOF) / entryOverhead);
            md.opticalDuplicateIndexes = new SortingLongCollection(maxInMemory, md.TMP_DIR.toArray(new File[md.TMP_DIR.size()]));
        }
        log.info("Will retain up to " + maxInMemory + " duplicate indices before spilling to disk.");
        md.duplicateIndexes = new SortingLongCollection(maxInMemory, md.TMP_DIR.toArray(new File[md.TMP_DIR.size()]));
        if (md.TAG_DUPLICATE_SET_MEMBERS) {
            final RepresentativeReadIndexerCodec representativeIndexCodec = new RepresentativeReadIndexerCodec();
            md.representativeReadIndicesForDuplicates = SortingCollection.newInstance(RepresentativeReadIndexer.class,
                    representativeIndexCodec,
                    Comparator.comparing(read -> read.readIndexInFile),
                    maxInMemory,
                    md.TMP_DIR);
        }

        ReadEndsForMarkDuplicates firstOfNextChunk = null;
        int nextChunkRead1Coordinate2Min = Integer.MAX_VALUE;
        int nextChunkRead1Coordinate2Max = Integer.MIN_VALUE;
        final List<ReadEndsForMarkDuplicates> nextChunk = new ArrayList<>(200);

        // First just do the pairs
        log.info("Traversing read pair information and detecting duplicates.");
        for (final ReadEndsForMarkDuplicates next : md.pairSort) {
            if (firstOfNextChunk != null && areComparableForDuplicatesWithEndSignificance(firstOfNextChunk, next, true, useBarcodes,
                    nextChunkRead1Coordinate2Min, nextChunkRead1Coordinate2Max)) {
                nextChunk.add(next);
                if ( next.read1Coordinate2 != END_INSIGNIFICANT_VALUE) {
                    nextChunkRead1Coordinate2Min = Math.min(nextChunkRead1Coordinate2Min, next.read1Coordinate2);
                    nextChunkRead1Coordinate2Max = Math.max(nextChunkRead1Coordinate2Max, next.read1Coordinate2);
                }
            } else {
                md.handleChunk(nextChunk);
                nextChunk.clear();
                nextChunk.add(next);
                firstOfNextChunk = next;
                if ( next.read1Coordinate2 != END_INSIGNIFICANT_VALUE)
                    nextChunkRead1Coordinate2Min = nextChunkRead1Coordinate2Max = next.read1Coordinate2;
                else {
                    nextChunkRead1Coordinate2Min = Integer.MAX_VALUE;
                    nextChunkRead1Coordinate2Max = Integer.MIN_VALUE;
                }
            }
        }
        md.handleChunk(nextChunk);

        md.pairSort.cleanup();
        md.pairSort = null;

        // Now deal with the fragments
        log.info("Traversing fragment information and detecting duplicates.");
        boolean containsPairs = false;
        boolean containsFrags = false;

        firstOfNextChunk = null;
        nextChunkRead1Coordinate2Min = Integer.MAX_VALUE;
        nextChunkRead1Coordinate2Max = Integer.MIN_VALUE;

        for (final ReadEndsForMarkDuplicates next : md.fragSort) {
            if (firstOfNextChunk != null && areComparableForDuplicatesWithEndSignificance(firstOfNextChunk, next, false, useBarcodes,
                    nextChunkRead1Coordinate2Min, nextChunkRead1Coordinate2Max)) {
                nextChunk.add(next);
                containsPairs = containsPairs || next.isPaired();
                containsFrags = containsFrags || !next.isPaired();
                if ( next.read1Coordinate2 != END_INSIGNIFICANT_VALUE) {
                    nextChunkRead1Coordinate2Min = Math.min(nextChunkRead1Coordinate2Min, next.read1Coordinate2);
                    nextChunkRead1Coordinate2Max = Math.max(nextChunkRead1Coordinate2Max, next.read1Coordinate2);

                    if ( firstOfNextChunk.read1Coordinate2 == END_INSIGNIFICANT_VALUE)
                        firstOfNextChunk = next;
                }
            } else {
                if (nextChunk.size() > 1 && containsFrags) {
                    md.markDuplicateFragments(nextChunk, containsPairs);
                }
                nextChunk.clear();
                nextChunk.add(next);
                firstOfNextChunk = next;
                if ( next.read1Coordinate2 != END_INSIGNIFICANT_VALUE)
                    nextChunkRead1Coordinate2Min = nextChunkRead1Coordinate2Max = next.read1Coordinate2;
                else {
                    nextChunkRead1Coordinate2Min = Integer.MAX_VALUE;
                    nextChunkRead1Coordinate2Max = Integer.MIN_VALUE;
                }
                containsPairs = next.isPaired();
                containsFrags = !next.isPaired();
            }
        }
        md.markDuplicateFragments(nextChunk, containsPairs);
        md.fragSort.cleanup();
        md.fragSort = null;

        log.info("Sorting list of duplicate records.");
        md.duplicateIndexes.doneAddingStartIteration();
        if (md.opticalDuplicateIndexes != null) {
            md.opticalDuplicateIndexes.doneAddingStartIteration();
        }
        if (md.TAG_DUPLICATE_SET_MEMBERS) {
            md.representativeReadIndicesForDuplicates.doneAdding();
        }
    }

    /**
     * Builds a read ends object that represents a single read - for flow based read
     */
    @Override
    public ReadEndsForMarkDuplicates buildReadEnds(final SAMFileHeader header, final long index, final SAMRecord rec, final boolean useBarcodes) {
        final ReadEndsForMarkDuplicates ends = md.buildReadEnds(header, index, rec, useBarcodes);

        // this code only supported unpaired reads
        if (rec.getReadPairedFlag() && !rec.getMateUnmappedFlag()) {
            throw new IllegalArgumentException("FLOW_MODE does not support paired reads. offending read: " + rec);
        }

        // adjust start/end coordinates
        ends.read1Coordinate = getReadEndCoordinate(rec, !rec.getReadNegativeStrandFlag(), true);
        if (md.USE_END_IN_UNPAIRED_READS) {
            ends.read1Coordinate2 = getReadEndCoordinate(rec, rec.getReadNegativeStrandFlag(), false);
        }

        // adjust score
        if ( md.FLOW_QUALITY_SUM_STRATEGY ) {
            ends.score = computeFlowDuplicateScore(rec, ends.read1Coordinate, ends.read1Coordinate2);
        }

        return ends;
    }

    /**
     * update score for pairedEnds
     */
    @Override
    public void updatePairedEndsScore(final SAMRecord rec, final ReadEndsForMarkDuplicates pairedEnds) {
        if (md. FLOW_QUALITY_SUM_STRATEGY ) {
            pairedEnds.score += computeFlowDuplicateScore(rec, pairedEnds.read1Coordinate, pairedEnds.read1Coordinate2);
        } else {
            md.updatePairedEndsScore(rec, pairedEnds);
        }
    }

    /**
     * This method is identical in function to areComparableForDuplicates except that it accomodates for
     * the possible significance of the end side of the reads (w/ or wo/ uncertainty). This is only
     * applicable for flow mode invocation.
     */
    private boolean areComparableForDuplicatesWithEndSignificance(final ReadEndsForMarkDuplicates lhs, final ReadEndsForMarkDuplicates rhs, final boolean compareRead2, final boolean useBarcodes,
                                                                  final int lhsRead1Coordinate2Min, final int lhsRead1Coordinate2Max) {
        boolean areComparable = md.areComparableForDuplicates(lhs, rhs, compareRead2, useBarcodes);

        if (areComparable) {
            areComparable = (!endCoorSignificant(lhs.read1Coordinate2, rhs.read1Coordinate2) ||
                    endCoorInRangeWithUncertainty(lhsRead1Coordinate2Min, lhsRead1Coordinate2Max,
                            rhs.read1Coordinate2, md.UNPAIRED_END_UNCERTAINTY));
        }

        return areComparable;
    }

    private boolean endCoorSignificant(final int lhsCoor, final int rhsCoor) {
        return lhsCoor != END_INSIGNIFICANT_VALUE && rhsCoor != END_INSIGNIFICANT_VALUE;
    }

    private boolean endCoorInRangeWithUncertainty(int lhsCoorMin, int lhsCoorMax, int rhsCoor, int uncertainty) {
        return (rhsCoor >= (lhsCoorMin - uncertainty)) && (rhsCoor <= (lhsCoorMax + uncertainty));
    }

    /**
     * A quality summing scoring strategy used for flow based reads.
     *
     * The method walks on the bases of the read, in the synthesis direction. For each base, the effective
     * quality value is defined as the value on the first base on the hmer to which the base belongs to. The score
     * is defined to be the sum of all effective values above a given threshold.
     *
     * @param rec - SAMRecord to get a score for
     * @param threshold - threshold above which effective quality is included
     * @return - calculated score (see method description)
     */
    static protected int getFlowSumOfBaseQualities(final SAMRecord rec, int threshold) {
        int score = 0;

        // access qualities and bases
        final byte[] quals = rec.getBaseQualities();
        final byte[]  bases = rec.getReadBases();

        // create iteration range and direction
        final int startingOffset = !rec.getReadNegativeStrandFlag() ? 0 : bases.length;
        final int endOffset = !rec.getReadNegativeStrandFlag() ? bases.length : 0;
        final int  iterIncr = !rec.getReadNegativeStrandFlag() ? 1 : -1;

        // loop on bases, extract qual related to homopolymer from start of homopolymer
        byte lastBase = 0;
        byte effectiveQual = 0;
        for ( int i = startingOffset ; i != endOffset ; i += iterIncr ) {
            final byte base = bases[i];
            if ( base != lastBase ) {
                effectiveQual = quals[i];
            }
            if ( effectiveQual >= threshold) {
                score += effectiveQual;
            }
            lastBase = base;
        }

        return score;
    }

    private short computeFlowDuplicateScore(SAMRecord rec, int start, int end) {

        if ( end == END_INSIGNIFICANT_VALUE)
            return -1;

        Short storedScore = (Short)rec.getTransientAttribute(ATTR_DUPLICATE_SCORE);
        if ( storedScore == null ) {
            short score = 0;

            score += (short) Math.min(getFlowSumOfBaseQualities(rec, md.FLOW_EFFECTIVE_QUALITY_THRESHOLD), Short.MAX_VALUE / 2);

            score += rec.getReadFailsVendorQualityCheckFlag() ? (short) (Short.MIN_VALUE / 2) : 0;
            storedScore = score;
            rec.setTransientAttribute(ATTR_DUPLICATE_SCORE, storedScore);
        }

        return storedScore;
    }
    private int getReadEndCoordinate(final SAMRecord rec, final boolean start, final boolean certain) {
        final FlowOrder flowOrder = new FlowOrder(rec);
        final int unclippedCoor = start ? rec.getUnclippedStart() : rec.getUnclippedEnd();
        final int alignmentCoor = start ? rec.getAlignmentStart() : rec.getAlignmentEnd();

        if ( !flowOrder.isValid() ) {
            return unclippedCoor;
        } else if ( certain && md.FLOW_SKIP_FIRST_N_FLOWS != 0 ) {
            final byte[] bases = rec.getReadBases();
            byte hmerBase = start ? bases[0] : bases[bases.length - 1];
            int  hmersLeft = md.FLOW_SKIP_FIRST_N_FLOWS;      // number of hmer left to trim

            // advance flow order to base
            while ( flowOrder.current() != hmerBase ) {
                flowOrder.advance();
                hmersLeft--;
            }

            int hmerSize;
            for ( hmerSize = 1; hmerSize < bases.length ; hmerSize++ ) {
                if ((start ? bases[hmerSize] : bases[bases.length - 1 - hmerSize]) != hmerBase) {
                    if (--hmersLeft <= 0) {
                        break;
                    } else {
                        hmerBase = start ? bases[hmerSize] : bases[bases.length - 1 - hmerSize];
                        flowOrder.advance();
                        while (flowOrder.current() != hmerBase) {
                            flowOrder.advance();
                            hmersLeft--;
                        }
                        if (hmersLeft <= 0) {
                            break;
                        }
                    }
                }
            }
            final int  coor = unclippedCoor + (start ? hmerSize : -hmerSize);
            return md.USE_UNPAIRED_CLIPPED_END
                    ? (start ? Math.max(coor, alignmentCoor) : Math.min(coor, alignmentCoor))
                    : coor;
        } else if (md. FLOW_Q_IS_KNOWN_END ? isAdapterClipped(rec) : isAdapterClippedWithQ(rec) ) {
            return unclippedCoor;
        } else if ( !certain && isQualityClipped(rec) ) {
            return END_INSIGNIFICANT_VALUE;
        } else if (md.USE_UNPAIRED_CLIPPED_END) {
            return alignmentCoor;
        } else {
            return unclippedCoor;
        }
    }

    public static boolean isAdapterClipped(final SAMRecord rec) {
        return clippingTagContains(rec, CLIPPING_TAG_CONTAINS_A);
    }

    public static boolean isAdapterClippedWithQ(final SAMRecord rec) {
        return clippingTagContains(rec, CLIPPING_TAG_CONTAINS_AQ);
    }

    public static boolean isQualityClipped(final SAMRecord rec) {
        return clippingTagContains(rec, CLIPPING_TAG_CONTAINS_QZ);
    }

    private static boolean clippingTagContains(final SAMRecord rec, final char[] chars) {
        final String clippingTagValue = rec.getStringAttribute(CLIPPING_TAG_NAME);

        if ( clippingTagValue == null ) {
            return false;
        } else {
            for ( final char ch : chars ) {
                if ( clippingTagValue.indexOf(ch) >= 0 ) {
                    return true;
                }
            }
            return false;
        }
    }

    /**
     * private class used to represent use a SAMRecord's flow order, if such is present
     */
    static private class FlowOrder {

        final byte[] flowOrder; // the flow order byte string
        int flowIndex = 0; // the current position on the flow order

        private FlowOrder(final SAMRecord rec) {

            // find flow order
            final SAMFileHeader header = rec.getHeader();
            for ( final SAMReadGroupRecord rg : header.getReadGroups() ) {
                if (rg.getFlowOrder() != null) {
                    flowOrder = rg.getFlowOrder().getBytes();
                    return;
                }
            }
            flowOrder = null;
        }

        private boolean isValid() {
            return flowOrder != null;
        }

        private void advance() {
            if (++flowIndex >= flowOrder.length) {
                flowIndex = 0;
            }
        }

        private byte current() {
            return flowOrder[flowIndex];
        }
    }
}
