package picard.sam.markduplicates;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.util.SortingLongCollection;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import picard.sam.markduplicates.util.ReadEndsForMarkDuplicates;
import picard.sam.markduplicates.util.RepresentativeReadIndexerCodec;
import picard.sam.util.RepresentativeReadIndexer;

import java.io.File;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/**
 * This class/tool extends the base MarkDuplicates with additional parameters and strategies.
 *
 * The tool should only be applied to flow base reads
 */
@CommandLineProgramProperties(
        summary = MarkDuplicates.USAGE_SUMMARY + MarkDuplicates.USAGE_DETAILS + MarkDuplicatesForFlow.FLOW_ONLY_SUFFIX,
        oneLineSummary = MarkDuplicates.USAGE_SUMMARY + MarkDuplicatesForFlow.FLOW_ONLY_SUFFIX,
        programGroup = ReadDataManipulationProgramGroup.class)
@DocumentedFeature
public class MarkDuplicatesForFlow extends MarkDuplicates {

    static final String FLOW_ONLY_SUFFIX = " (flow based read version)";
    private final Log log = Log.getInstance(MarkDuplicatesForFlow.class);


    @Argument(doc = "Use specific quality summing strategy for flow based reads. The strategy ensures that the same " +
            "(and correct) quality value is used for all bases of the same homopolymer. Default false.")
    public boolean FLOW_QUALITY_SUM_STRATEGY = false;

    @Argument(doc = "Make the end location of single end read be significant when considering duplicates, " +
            "in addition to the start location, which is always significant (i.e. require single end reads to start and" +
            "end on the same position to be considered duplicate). Default false.")
    public boolean USE_END_IN_UNPAIRED_READS = false;

    @Argument(doc = "Use position of the clipping as the end position, when considering duplicates (or use the unclipped end position). Default false.")
    public boolean USE_UNPAIRED_CLIPPED_END = false;

    @Argument(doc = "Maximal difference of the read end position that counted as equal. Useful for flow based " +
            "reads where the end position might vary due to sequencing errors. Default 0.")
    public int UNPAIRED_END_UNCERTAINTY = 0;

    @Argument(doc = "Skip first N flows, when considering duplicates. Useful for flow based reads where sometimes there " +
            "is noise in the first flows. Default 0.")
    public int FLOW_SKIP_FIRST_N_FLOWS = 0;

    @Argument(doc = "Treat position of read trimming based on quality as the known end (relevant for flow based reads). Default false - if the read " +
            "is trimmed on quality its end is not defined and the read is duplicate with any read starting at the same place. ")
    public boolean FLOW_Q_IS_KNOWN_END = false;

    // constants for clippingTagContains
    public static String        CLIPPING_TAG_NAME = "tm";
    public static final char[]  CLIPPING_TAG_CONTAINS_A = {'A'};
    public static final char[]  CLIPPING_TAG_CONTAINS_AQ = {'A', 'Q'};
    public static final char[]  CLIPPING_TAG_CONTAINS_QZ = {'Q', 'Z'};

    /**
     * This method is identical in function to generateDuplicateIndexes except that it accomodates for
     * the possible significance of the end side of the reads (w/ or wo/ uncertainty). This is only
     * applicable for flow mode invocation.
     */
    @Override
    protected void generateDuplicateIndexes(final boolean useBarcodes, final boolean indexOpticalDuplicates) {
        final int entryOverhead;
        if (TAG_DUPLICATE_SET_MEMBERS) {
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
            this.opticalDuplicateIndexes = new SortingLongCollection(maxInMemory, TMP_DIR.toArray(new File[TMP_DIR.size()]));
        }
        log.info("Will retain up to " + maxInMemory + " duplicate indices before spilling to disk.");
        this.duplicateIndexes = new SortingLongCollection(maxInMemory, TMP_DIR.toArray(new File[TMP_DIR.size()]));
        if (TAG_DUPLICATE_SET_MEMBERS) {
            final RepresentativeReadIndexerCodec representativeIndexCodec = new RepresentativeReadIndexerCodec();
            this.representativeReadIndicesForDuplicates = SortingCollection.newInstance(RepresentativeReadIndexer.class,
                    representativeIndexCodec,
                    Comparator.comparing(read -> read.readIndexInFile),
                    maxInMemory,
                    TMP_DIR);
        }

        ReadEndsForMarkDuplicates firstOfNextChunk = null;
        int nextChunkRead1Coordinate2Min = Integer.MAX_VALUE;
        int nextChunkRead1Coordinate2Max = Integer.MIN_VALUE;
        final List<ReadEndsForMarkDuplicates> nextChunk = new ArrayList<>(200);

        // First just do the pairs
        log.info("Traversing read pair information and detecting duplicates.");
        for (final ReadEndsForMarkDuplicates next : this.pairSort) {
            if (firstOfNextChunk != null && areComparableForDuplicatesWithEndSignificance(firstOfNextChunk, next, true, useBarcodes,
                    nextChunkRead1Coordinate2Min, nextChunkRead1Coordinate2Max)) {
                nextChunk.add(next);
                if ( next.read1Coordinate2 != END_INSIGNIFICANT ) {
                    nextChunkRead1Coordinate2Min = Math.min(nextChunkRead1Coordinate2Min, next.read1Coordinate2);
                    nextChunkRead1Coordinate2Max = Math.max(nextChunkRead1Coordinate2Max, next.read1Coordinate2);
                }
            } else {
                handleChunk(nextChunk);
                nextChunk.clear();
                nextChunk.add(next);
                firstOfNextChunk = next;
                if ( next.read1Coordinate2 != END_INSIGNIFICANT )
                    nextChunkRead1Coordinate2Min = nextChunkRead1Coordinate2Max = next.read1Coordinate2;
                else {
                    nextChunkRead1Coordinate2Min = Integer.MAX_VALUE;
                    nextChunkRead1Coordinate2Max = Integer.MIN_VALUE;
                }
            }
        }
        handleChunk(nextChunk);

        this.pairSort.cleanup();
        this.pairSort = null;

        // Now deal with the fragments
        log.info("Traversing fragment information and detecting duplicates.");
        boolean containsPairs = false;
        boolean containsFrags = false;

        firstOfNextChunk = null;
        nextChunkRead1Coordinate2Min = Integer.MAX_VALUE;
        nextChunkRead1Coordinate2Max = Integer.MIN_VALUE;

        for (final ReadEndsForMarkDuplicates next : this.fragSort) {
            if (firstOfNextChunk != null && areComparableForDuplicatesWithEndSignificance(firstOfNextChunk, next, false, useBarcodes,
                    nextChunkRead1Coordinate2Min, nextChunkRead1Coordinate2Max)) {
                nextChunk.add(next);
                containsPairs = containsPairs || next.isPaired();
                containsFrags = containsFrags || !next.isPaired();
                if ( next.read1Coordinate2 != END_INSIGNIFICANT ) {
                    nextChunkRead1Coordinate2Min = Math.min(nextChunkRead1Coordinate2Min, next.read1Coordinate2);
                    nextChunkRead1Coordinate2Max = Math.max(nextChunkRead1Coordinate2Max, next.read1Coordinate2);

                    if ( firstOfNextChunk.read1Coordinate2 == END_INSIGNIFICANT )
                        firstOfNextChunk = next;
                }
            } else {
                if (nextChunk.size() > 1 && containsFrags) {
                    markDuplicateFragments(nextChunk, containsPairs);
                }
                nextChunk.clear();
                nextChunk.add(next);
                firstOfNextChunk = next;
                if ( next.read1Coordinate2 != END_INSIGNIFICANT )
                    nextChunkRead1Coordinate2Min = nextChunkRead1Coordinate2Max = next.read1Coordinate2;
                else {
                    nextChunkRead1Coordinate2Min = Integer.MAX_VALUE;
                    nextChunkRead1Coordinate2Max = Integer.MIN_VALUE;
                }
                containsPairs = next.isPaired();
                containsFrags = !next.isPaired();
            }
        }
        markDuplicateFragments(nextChunk, containsPairs);
        this.fragSort.cleanup();
        this.fragSort = null;

        log.info("Sorting list of duplicate records.");
        this.duplicateIndexes.doneAddingStartIteration();
        if (this.opticalDuplicateIndexes != null) {
            this.opticalDuplicateIndexes.doneAddingStartIteration();
        }
        if (TAG_DUPLICATE_SET_MEMBERS) {
            this.representativeReadIndicesForDuplicates.doneAdding();
        }
    }

    /**
     * Builds a read ends object that represents a single read - for flow based read
     */
    @Override
    protected ReadEndsForMarkDuplicates buildReadEnds(final SAMFileHeader header, final long index, final SAMRecord rec, final boolean useBarcodes) {
        final ReadEndsForMarkDuplicates ends = super.buildReadEnds(header, index, rec, useBarcodes);

        // this code only supported unpaired reads
        if (rec.getReadPairedFlag() && !rec.getMateUnmappedFlag()) {
            throw new IllegalArgumentException("FLOW_MODE does not support paired reads. offending read: " + rec);
        }

        // adjust start/end coordinates
        ends.read1Coordinate = getReadEndCoordinate(rec, !rec.getReadNegativeStrandFlag(), true);
        if (USE_END_IN_UNPAIRED_READS) {
            ends.read1Coordinate2 = getReadEndCoordinate(rec, rec.getReadNegativeStrandFlag(), false);
        }

        // adjust score
        if ( FLOW_QUALITY_SUM_STRATEGY ) {
            ends.score = computeFlowDuplicateScore(rec, ends.read1Coordinate, ends.read1Coordinate2);
        }

        return ends;
    }

    /**
     * This method is identical in function to areComparableForDuplicates except that it accomodates for
     * the possible significance of the end side of the reads (w/ or wo/ uncertainty). This is only
     * applicable for flow mode invocation.
     */
    private boolean areComparableForDuplicatesWithEndSignificance(final ReadEndsForMarkDuplicates lhs, final ReadEndsForMarkDuplicates rhs, final boolean compareRead2, final boolean useBarcodes,
                                                                  final int lhsRead1Coordinate2Min, final int lhsRead1Coordinate2Max) {
        boolean areComparable = areComparableForDuplicates(lhs, rhs, compareRead2, useBarcodes);

        if (areComparable) {
            areComparable = (!endCoorSignificant(lhs.read1Coordinate2, rhs.read1Coordinate2) ||
                    endCoorInRangeWithUncertainty(lhsRead1Coordinate2Min, lhsRead1Coordinate2Max,
                            rhs.read1Coordinate2, UNPAIRED_END_UNCERTAINTY));
        }

        return areComparable;
    }

    private boolean endCoorSignificant(final int lhsCoor, final int rhsCoor) {
        return lhsCoor != END_INSIGNIFICANT && rhsCoor != END_INSIGNIFICANT;
    }

    private boolean endCoorInRangeWithUncertainty(int lhsCoorMin, int lhsCoorMax, int rhsCoor, int uncertainty) {
        return (rhsCoor >= (lhsCoorMin - uncertainty)) && (rhsCoor <= (lhsCoorMax + uncertainty));
    }

    private int getFlowSumOfBaseQualities(final SAMRecord rec) {
        int score = 0;

        // access qualities and bases
        final byte[]      quals = rec.getBaseQualities();
        final byte[]      bases = rec.getReadBases();

        // create iteration range and direction
        final int         startingOffset = !rec.getReadNegativeStrandFlag() ? 0 : bases.length;
        final int         endOffset = !rec.getReadNegativeStrandFlag() ? bases.length : 0;
        final int         iterIncr = !rec.getReadNegativeStrandFlag() ? 1 : -1;

        // loop on bases, extract qual related to homopolymer from start of homopolymer
        byte        lastBase = 0;
        byte        effectiveQual = 0;
        for ( int i = startingOffset ; i != endOffset ; i += iterIncr ) {
            final byte        base = bases[i];
            if ( base != lastBase ) {
                effectiveQual = quals[i];
            }
            if ( effectiveQual >= FLOW_EFFECTIVE_QUALITY_THRESHOLD) {
                score += effectiveQual;
            }
            lastBase = base;
        }

        return score;
    }

    /**
     * update score for pairedEnds
     */
    @Override
    protected void updatePairedEndsScore(SAMRecord rec, ReadEndsForMarkDuplicates pairedEnds) {
        if ( FLOW_QUALITY_SUM_STRATEGY ) {
            pairedEnds.score += computeFlowDuplicateScore(rec, pairedEnds.read1Coordinate, pairedEnds.read1Coordinate2);
        } else {
            super.updatePairedEndsScore(rec, pairedEnds);
        }
    }

    private short computeFlowDuplicateScore(SAMRecord rec, int start, int end) {

        if ( end == END_INSIGNIFICANT )
            return -1;

        Short storedScore = (Short)rec.getTransientAttribute(ATTR_DUPLICATE_SCORE);
        if ( storedScore == null ) {
            short score = 0;

            score += (short) Math.min(getFlowSumOfBaseQualities(rec), Short.MAX_VALUE / 2);

            score += rec.getReadFailsVendorQualityCheckFlag() ? (short) (Short.MIN_VALUE / 2) : 0;
            storedScore = score;
            rec.setTransientAttribute(ATTR_DUPLICATE_SCORE, storedScore);
        }

        return storedScore;
    }
    private int getReadEndCoordinate(final SAMRecord rec, final boolean start, final boolean certain) {
        final FlowOrder     flowOrder = new FlowOrder(rec);
        final int           unclippedCoor = start ? rec.getUnclippedStart() : rec.getUnclippedEnd();
        final int           alignmentCoor = start ? rec.getAlignmentStart() : rec.getAlignmentEnd();

        if ( !flowOrder.isValid() ) {
            return unclippedCoor;
        } else if ( certain && FLOW_SKIP_FIRST_N_FLOWS != 0 ) {
            final byte[] bases = rec.getReadBases();
            byte        hmerBase = start ? bases[0] : bases[bases.length - 1];
            int         hmersLeft = FLOW_SKIP_FIRST_N_FLOWS;      // number of hmer left to trim

            // advance flow order to base
            while ( flowOrder.current() != hmerBase ) {
                flowOrder.advance();
                hmersLeft--;
            }

            int         hmerSize;
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
            final int     coor = unclippedCoor + (start ? hmerSize : -hmerSize);
            return USE_UNPAIRED_CLIPPED_END
                    ? (start ? Math.max(coor, alignmentCoor) : Math.min(coor, alignmentCoor))
                    : coor;
        } else if ( FLOW_Q_IS_KNOWN_END ? isAdapterClipped(rec) : isAdapterClippedWithQ(rec) ) {
            return unclippedCoor;
        } else if ( !certain && isQualityClipped(rec) ) {
            return END_INSIGNIFICANT;
        } else if (USE_UNPAIRED_CLIPPED_END) {
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
