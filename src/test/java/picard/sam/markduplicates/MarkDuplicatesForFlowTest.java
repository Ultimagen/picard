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

    private interface TesterModifier {
        void modify(AbstractMarkDuplicatesCommandLineProgramTester tester);
    }

    private static class TestRecordInfo {
        int             length;
        int             alignmentStart;
        String          cigar;
        boolean         isDuplicate;
        String          startMod;
        String          endMod;

        TestRecordInfo(final int length, final int alignmentStart, final String cigar, final boolean isDuplicate,
                       final String startMod, final String endMod) {
            this.length = length;
            this.alignmentStart = alignmentStart;
            this.cigar = cigar;
            this.isDuplicate = isDuplicate;
            this.startMod = startMod;
            this.endMod = endMod;
        }
    }

    @DataProvider(name ="forFlowDataProvider")
    public Object[][] dataProvider() {
        return new Object[][] {
                // testUSE_END_IN_UNPAIRED_READS: End location is not significant
                {
                    null,
                    new TestRecordInfo[] {
                            new TestRecordInfo(76, 12, null, false, null, null),
                            new TestRecordInfo(74, 12, null, true, null, null)
                    },
                    new String[] { "USE_END_IN_UNPAIRED_READS=false" }, null
                },

                // testUSE_END_IN_UNPAIRED_READS: End location is significant
                {
                        null,
                        new TestRecordInfo[] {
                                new TestRecordInfo(76, 12, null, false, null, null),
                                new TestRecordInfo(74, 12, null, false, null, null)
                        },
                        new String[] { "USE_END_IN_UNPAIRED_READS=true" }, null
                },

                // testUSE_UNPAIRED_CLIPPED_END: Do not use clipped locations (meaning, use unclipped)
                {
                        null,
                        new TestRecordInfo[] {
                                new TestRecordInfo(76, 12, "1S76M", false, null, null),
                                new TestRecordInfo(74, 12, "74M", false, null, null)
                        },
                        new String[] { "USE_UNPAIRED_CLIPPED_END=false" }, null
                },

                // testUSE_UNPAIRED_CLIPPED_END: Use clipped locations (meaning, use clipped)
                {
                        null,
                        new TestRecordInfo[] {
                                new TestRecordInfo(76, 12, "1S76M", false, null, null),
                                new TestRecordInfo(74, 12, "74M", true, null, null)
                        },
                        new String[] { "USE_UNPAIRED_CLIPPED_END=true" }, null
                },

                // testUSE_UNPAIRED_CLIPPED_END: Use clipped locations (meaning, use clipped)
                {
                        null,
                        new TestRecordInfo[] {
                                new TestRecordInfo(76, 12, "1S76M1S", true, null, null),
                                new TestRecordInfo(78, 11, "78M", false, null, null)
                        },
                        new String[] { "USE_UNPAIRED_CLIPPED_END=false", "USE_END_IN_UNPAIRED_READS=true" }, null
                },

                // testUSE_UNPAIRED_CLIPPED_END: Use clipped locations (meaning, use clipped)
                {
                        null,
                        new TestRecordInfo[] {
                                new TestRecordInfo(76, 12, "1S76M1S", false, null, null),
                                new TestRecordInfo(78, 11, "78M", false, null, null)
                        },
                        new String[] { "USE_UNPAIRED_CLIPPED_END=true", "USE_END_IN_UNPAIRED_READS=true" }, null
                },

                // testFLOW_SKIP_FIRST_N_FLOWS: Do not use clipped locations (meaning, use unclipped)
                {
                        null,
                        new TestRecordInfo[] {
                                new TestRecordInfo(76, 12,"76M", false, "ACGTT", "TTGCA"),
                                new TestRecordInfo(76, 12, "76M", true, "ACGGT", "TGGCA")
                        },
                        new String[] { "USE_END_IN_UNPAIRED_READS=true", "FLOW_SKIP_FIRST_N_FLOWS=0" }, null
                },

                // testFLOW_SKIP_FIRST_N_FLOWS: Do not use clipped locations (meaning, use unclipped)
                {
                        null,
                        new TestRecordInfo[] {
                                new TestRecordInfo(76, 12,"76M", false, "ACGTT", "TTGCA"),
                                new TestRecordInfo(76, 12, "76M", false, "CCGGT", "TGGCA")
                        },
                        new String[] { "USE_END_IN_UNPAIRED_READS=true", "FLOW_SKIP_FIRST_N_FLOWS=3" }, null
                },

                // testFLOW_QUALITY_SUM_STRATEGY: normal sum
                {
                        DuplicateScoringStrategy.ScoringStrategy.SUM_OF_BASE_QUALITIES,
                        new TestRecordInfo[] {
                                new TestRecordInfo(76, 12,"76M", true, "AAAC", null),
                                new TestRecordInfo(76, 12, "76M", false, "AACC", null)
                        },
                        new String[] { "FLOW_QUALITY_SUM_STRATEGY=false" },
                        new TesterModifier() {
                            @Override
                            public void modify(final AbstractMarkDuplicatesCommandLineProgramTester tester) {
                                final SAMRecord[] records = tester.getSamRecordSetBuilder().getRecords().toArray(new SAMRecord[0]);
                                records[0].setAttribute("tp", new int[76]);
                                records[1].setAttribute("tp", new int[76]);
                                records[0].getBaseQualities()[1] = 25; // dip inside AAA
                            }
                        }
                },

                // testFLOW_QUALITY_SUM_STRATEGY: flow (homopolymer based) sum
                {
                        DuplicateScoringStrategy.ScoringStrategy.SUM_OF_BASE_QUALITIES,
                        new TestRecordInfo[] {
                                new TestRecordInfo(76, 12,"76M", false, "AAAC", null),
                                new TestRecordInfo(76, 12, "76M", true, "AACC", null)
                        },
                        new String[] { "FLOW_QUALITY_SUM_STRATEGY=true" },
                        new TesterModifier() {
                            @Override
                            public void modify(final AbstractMarkDuplicatesCommandLineProgramTester tester) {
                                final SAMRecord[] records = tester.getSamRecordSetBuilder().getRecords().toArray(new SAMRecord[0]);
                                records[0].setAttribute("tp", new int[76]);
                                records[1].setAttribute("tp", new int[76]);
                                records[0].getBaseQualities()[1] = 25; // dip inside AAA
                            }
                        }
                },

        };
    }

    @Test(dataProvider = "forFlowDataProvider")
    public void testForFlow(final DuplicateScoringStrategy.ScoringStrategy scoringStrategy, final TestRecordInfo[] recInfos, final String[] params, TesterModifier modifier) {

        // get tester, build records
        final AbstractMarkDuplicatesCommandLineProgramTester tester =
                scoringStrategy == null ? getTester() : new MarkDuplicatesForFlowTester(scoringStrategy);
        for ( final TestRecordInfo info : recInfos ) {
            tester.getSamRecordSetBuilder().setReadLength(info.length);
            if ( info.cigar != null ) {
                tester.addMappedFragment(0, info.alignmentStart, info.isDuplicate, info.cigar, 50);
            } else {
                tester.addMappedFragment(0, info.alignmentStart, info.isDuplicate, 50);
            }
        }

        // modify records
        final SAMRecord[] records = tester.getSamRecordSetBuilder().getRecords().toArray(new SAMRecord[0]);
        for ( int i = 0 ; i < records.length ; i++ ) {
            final SAMRecord       rec = records[i];
            final TestRecordInfo  info = recInfos[i];

            if ( info.startMod != null ) {
                System.arraycopy(info.startMod.getBytes(), 0, rec.getReadBases(), 0, info.startMod.length());
            }
            if ( info.endMod != null ) {
                System.arraycopy(info.endMod.getBytes(), 0, rec.getReadBases(), rec.getReadBases().length - info.endMod.length(), info.endMod.length());
            }
        }

        // add parames, set flow order
        for ( final String param : params ) {
            tester.addArg(param);
        }
        tester.getSamRecordSetBuilder().getHeader().getReadGroups().get(0).setFlowOrder(FLOW_ORDER);

        // further modify tester
        if ( modifier != null )
            modifier.modify(tester);

        // run test
        tester.runTest();
    }
}
