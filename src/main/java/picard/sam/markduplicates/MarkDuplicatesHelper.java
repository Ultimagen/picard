package picard.sam.markduplicates;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import picard.sam.markduplicates.util.ReadEndsForMarkDuplicates;

public interface MarkDuplicatesHelper {

    void generateDuplicateIndexes(boolean useBarcodes, boolean indexOpticalDuplicates);
    ReadEndsForMarkDuplicates buildReadEnds(SAMFileHeader header,long index, final SAMRecord rec, boolean useBarcodes);
    short getReadDuplicateScore(SAMRecord rec, ReadEndsForMarkDuplicates pairedEnds);
}
