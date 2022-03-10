package picard.sam;

import htsjdk.samtools.SAMReadGroupRecord;
import picard.sam.markduplicates.util.LibraryIdGenerator;

import java.util.List;

public class DuplicationMetricsFactory {

    public static DuplicationMetrics createForLibrary(String library, List<SAMReadGroupRecord> readGroups) {

        // scan read groups, find read group associated with library
        for ( final SAMReadGroupRecord readGroup : readGroups ) {
            if ( library.equals(LibraryIdGenerator.getReadGroupLibraryName(readGroup)) ) {
                return createForLibrary(readGroup);
            }
        }

        // if here library not found, create default
        return new DuplicationMetrics();
    }

    // create a DuplicationMetrics for a specific read group
    public static DuplicationMetrics createForLibrary(SAMReadGroupRecord readGroup) {

        // create based on the presence of flow order
        if ( readGroup.getFlowOrder() == null ) {
            return new DuplicationMetrics();
        } else {
            return new FlowBasedDuplicationMetrics();
        }
    }

    public static DuplicationMetrics createDefault() {
        return new DuplicationMetrics();
    }
}
