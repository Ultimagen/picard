/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package picard.sam.markduplicates;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IterableAdapter;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

/**
 * This class defines the individual test cases to run. The actual running of the test is done
 * by MarkDuplicatesTester (see getTester).
 */
public class MarkDuplicatesForFlowTest  {
    protected static String TEST_BASE_NAME = null;
    protected static File TEST_DATA_DIR = null;
    final static String   FLOW_ORDER = "TGCA";


    @BeforeClass
    public void setUp() {
        TEST_BASE_NAME = "MarkDuplicatesForFlow";
        TEST_DATA_DIR = new File("testdata/picard/sam/MarkDuplicates");
    }

    protected AbstractMarkDuplicatesCommandLineProgramTester getTester() {
        return new MarkDuplicatesForFlowTester();
    }

    @Test
    public void testUSE_END_IN_UNPAIRED_READS() {
        AbstractMarkDuplicatesCommandLineProgramTester tester;

        // End location is not significant
        tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(76);
        tester.addMappedFragment(0, 12, false, 50);
        tester.getSamRecordSetBuilder().setReadLength(74);
        tester.addMappedFragment(0, 12, true, 50);
        tester.addArg("USE_END_IN_UNPAIRED_READS=false");
        tester.getSamRecordSetBuilder().getHeader().getReadGroups().get(0).setFlowOrder(FLOW_ORDER);
        tester.runTest();

        // End location is significant
        tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(76);
        tester.addMappedFragment(0, 12, false, 50);
        tester.getSamRecordSetBuilder().setReadLength(74);
        tester.addMappedFragment(0, 12, false, 50);
        tester.addArg("USE_END_IN_UNPAIRED_READS=true");
        tester.getSamRecordSetBuilder().getHeader().getReadGroups().get(0).setFlowOrder(FLOW_ORDER);
        tester.runTest();
    }

    @Test
    public void testUSE_UNPAIRED_CLIPPED_END() {
        AbstractMarkDuplicatesCommandLineProgramTester tester;

        // Do not use clipped locations (meaning, use unclipped)
        tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(76);
        tester.addMappedFragment(0, 12, false, "1S76M", 50);
        tester.getSamRecordSetBuilder().setReadLength(74);
        tester.addMappedFragment(0, 12, false, "74M", 50);
        tester.addArg("USE_UNPAIRED_CLIPPED_END=false");
        tester.getSamRecordSetBuilder().getHeader().getReadGroups().get(0).setFlowOrder(FLOW_ORDER);
        tester.runTest();

        // Use clipped locations (meaning, use clipped)
        tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(76);
        tester.addMappedFragment(0, 12, false, "1S76M", 50);
        tester.getSamRecordSetBuilder().setReadLength(74);
        tester.addMappedFragment(0, 12, true, "74M", 50);
        tester.addArg("USE_UNPAIRED_CLIPPED_END=true");
        tester.getSamRecordSetBuilder().getHeader().getReadGroups().get(0).setFlowOrder(FLOW_ORDER);
        tester.runTest();

        // Use clipped locations (meaning, use clipped)
        tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(76);
        tester.addMappedFragment(0, 12, true, "1S76M1S", 50);
        tester.getSamRecordSetBuilder().setReadLength(78);
        tester.addMappedFragment(0, 11, false, "78M", 50);
        tester.addArg("USE_UNPAIRED_CLIPPED_END=false");
        tester.addArg("USE_END_IN_UNPAIRED_READS=true");
        tester.getSamRecordSetBuilder().getHeader().getReadGroups().get(0).setFlowOrder(FLOW_ORDER);
        tester.runTest();

        // Use clipped locations (meaning, use clipped)
        tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(76);
        tester.addMappedFragment(0, 12, false, "1S76M1S", 50);
        tester.getSamRecordSetBuilder().setReadLength(78);
        tester.addMappedFragment(0, 11, false, "78M", 50);
        tester.addArg("USE_UNPAIRED_CLIPPED_END=true");
        tester.addArg("USE_END_IN_UNPAIRED_READS=true");
        tester.getSamRecordSetBuilder().getHeader().getReadGroups().get(0).setFlowOrder(FLOW_ORDER);
        tester.runTest();
    }

    @Test
    public void testFLOW_SKIP_FIRST_N_FLOWS() {
        AbstractMarkDuplicatesCommandLineProgramTester tester;
        SAMRecord[]     records;

        // Do not use clipped locations (meaning, use unclipped)
        tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(76);
        tester.addMappedFragment(0, 12, false, "76M", 50);
        tester.addMappedFragment(0, 12, true, "76M", 50);
        records = tester.getSamRecordSetBuilder().getRecords().toArray(new SAMRecord[0]);
        System.arraycopy("ACGTT".getBytes(), 0, records[0].getReadBases(), 0, 5);
        System.arraycopy("TTGCA".getBytes(), 0, records[0].getReadBases(), records[0].getReadBases().length - 5, 4);
        System.arraycopy("ACGGT".getBytes(), 0, records[1].getReadBases(), 0, 5);
        System.arraycopy("TGGCA".getBytes(), 0, records[1].getReadBases(), records[1].getReadBases().length - 5, 4);
        tester.addArg("USE_END_IN_UNPAIRED_READS=true");
        tester.addArg("FLOW_SKIP_FIRST_N_FLOWS=0");
        tester.getSamRecordSetBuilder().getHeader().getReadGroups().get(0).setFlowOrder(FLOW_ORDER);
        tester.runTest();

        // Do not use clipped locations (meaning, use unclipped)
        tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(76);
        tester.addMappedFragment(0, 12, false, "76M", 50);
        tester.addMappedFragment(0, 12, false, "76M", 50);
        records = tester.getSamRecordSetBuilder().getRecords().toArray(new SAMRecord[0]);
        System.arraycopy("ACGTT".getBytes(), 0, records[0].getReadBases(), 0, 5);
        System.arraycopy("TTGCA".getBytes(), 0, records[0].getReadBases(), records[0].getReadBases().length - 5, 4);
        System.arraycopy("CCGGT".getBytes(), 0, records[1].getReadBases(), 0, 5);
        System.arraycopy("TGGCA".getBytes(), 0, records[1].getReadBases(), records[1].getReadBases().length - 5, 4);
        tester.addArg("USE_END_IN_UNPAIRED_READS=true");
        tester.addArg("FLOW_SKIP_FIRST_N_FLOWS=3");
        tester.getSamRecordSetBuilder().getHeader().getReadGroups().get(0).setFlowOrder(FLOW_ORDER);
        tester.runTest();
    }

    @Test
    public void testFLOW_QUALITY_SUM_STRATEGY() {
        AbstractMarkDuplicatesCommandLineProgramTester tester;
        SAMRecord[]     records;

        //  normal sum
        tester = new MarkDuplicatesForFlowTester(DuplicateScoringStrategy.ScoringStrategy.SUM_OF_BASE_QUALITIES);
        tester.getSamRecordSetBuilder().setReadLength(76);
        tester.addMappedFragment(0, 12, true, "76M", 50);
        tester.addMappedFragment(0, 12, false, "76M", 50);
        records = tester.getSamRecordSetBuilder().getRecords().toArray(new SAMRecord[0]);
        records[0].setAttribute("tp", new int[76]);
        records[1].setAttribute("tp", new int[76]);
        records[0].getBaseQualities()[1] = 25; // dip inside AAA
        System.arraycopy("AAAC".getBytes(), 0, records[0].getReadBases(), 0, 4);
        System.arraycopy("AACC".getBytes(), 0, records[1].getReadBases(), 0, 4);
        tester.addArg("FLOW_QUALITY_SUM_STRATEGY=false");
        tester.getSamRecordSetBuilder().getHeader().getReadGroups().get(0).setFlowOrder(FLOW_ORDER);
        tester.runTest();

        // flow (homopolymer based) sum
        tester = new MarkDuplicatesForFlowTester(DuplicateScoringStrategy.ScoringStrategy.SUM_OF_BASE_QUALITIES);
        tester.getSamRecordSetBuilder().setReadLength(76);
        tester.addMappedFragment(0, 12, false, "76M", 50);
        tester.addMappedFragment(0, 12, true, "76M", 50);
        records = tester.getSamRecordSetBuilder().getRecords().toArray(new SAMRecord[0]);
        records[0].setAttribute("tp", new int[76]);
        records[1].setAttribute("tp", new int[76]);
        records[0].getBaseQualities()[1] = 25; // dip inside AAA
        System.arraycopy("AAAC".getBytes(), 0, records[0].getReadBases(), 0, 4);
        System.arraycopy("AACC".getBytes(), 0, records[1].getReadBases(), 0, 4);
        tester.addArg("FLOW_QUALITY_SUM_STRATEGY=true");
        tester.getSamRecordSetBuilder().getHeader().getReadGroups().get(0).setFlowOrder(FLOW_ORDER);
        tester.runTest();
    }
}
