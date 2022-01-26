/*
 * The MIT License
 *
 * Copyright (c) 2015 The Broad Institute
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

package picard.analysis;

import htsjdk.samtools.metrics.MetricsFile;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.function.Consumer;

/**
 * Created by kbergin on 11/23/15.
 */
public class CollectQualityYieldMetricsTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/sam/");

    public String getCommandLineProgramName() {
        return CollectQualityYieldMetrics.class.getSimpleName();
    }

    @DataProvider(name = "testSimpleExamplesDataProvider")
    public Object[][] testSimpleExamplesDataProvider() {
        return new Object[][] {
                {
                        "collect_quality_yield_test1.sam",
                        new Consumer<CollectQualityYieldMetrics.QualityYieldMetrics>() {
                            @Override
                            public void accept(final CollectQualityYieldMetrics.QualityYieldMetrics metrics) {
                                Assert.assertEquals(metrics.TOTAL_READS, 3);
                                Assert.assertEquals(metrics.PF_READS, 2);
                                Assert.assertEquals(metrics.READ_LENGTH, 10);
                                Assert.assertEquals(metrics.TOTAL_BASES, 30);
                                Assert.assertEquals(metrics.PF_BASES, 20);
                                Assert.assertEquals(metrics.Q20_BASES, 30);
                                Assert.assertEquals(metrics.PF_Q20_BASES, 20);
                                Assert.assertEquals(metrics.Q30_BASES, 5);
                                Assert.assertEquals(metrics.PF_Q30_BASES, 5);
                                Assert.assertEquals(metrics.Q20_EQUIVALENT_YIELD, 40);
                                Assert.assertEquals(metrics.PF_Q20_EQUIVALENT_YIELD, 27);
                                Assert.assertEquals(metrics.READ_LENGTH_AVG_Q_ABOVE_30, 0);
                                Assert.assertEquals(metrics.READ_LENGTH_AVG_Q_ABOVE_25, 0);
                            }
                        }
                },
        };
    }

    @Test(dataProvider = "testSimpleExamplesDataProvider")
    public void testSimpleExamples(final String inputFilename, Consumer<CollectQualityYieldMetrics.QualityYieldMetrics> metricsConsumer) throws IOException {
        final File input = new File(TEST_DATA_DIR, "CollectQualityYieldMetrics/" + inputFilename);
        final File outfile   = File.createTempFile("test", ".quality_yield_metrics");
        outfile.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
        };

        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<CollectQualityYieldMetrics.QualityYieldMetrics, ?> output = new MetricsFile<>();
        output.read(new FileReader(outfile));

        Assert.assertEquals(output.getMetrics().size(),1);

        final CollectQualityYieldMetrics.QualityYieldMetrics metrics = output.getMetrics().get(0);
        metricsConsumer.accept(metrics);
    }

    @Test
    public void testIntegration() throws IOException {
        final File input = new File(TEST_DATA_DIR, "insert_size_metrics_test.sam");
        final File outfile   = File.createTempFile("test", ".quality_yield_metrics");
        outfile.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
        };

        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<CollectQualityYieldMetrics.QualityYieldMetrics, ?> output = new MetricsFile<>();
        output.read(new FileReader(outfile));

        Assert.assertEquals(output.getMetrics().size(),1);

        final CollectQualityYieldMetrics.QualityYieldMetrics metrics = output.getMetrics().get(0);
        Assert.assertEquals(metrics.TOTAL_READS, 52);
        Assert.assertEquals(metrics.PF_READS, 52);
        Assert.assertEquals(metrics.READ_LENGTH, 101);
        Assert.assertEquals(metrics.TOTAL_BASES, 5252);
        Assert.assertEquals(metrics.PF_BASES, 5252);
        Assert.assertEquals(metrics.Q20_BASES, 3532);
        Assert.assertEquals(metrics.PF_Q20_BASES, 3532);
        Assert.assertEquals(metrics.Q30_BASES, 3145);
        Assert.assertEquals(metrics.PF_Q30_BASES, 3145);
        Assert.assertEquals(metrics.Q20_EQUIVALENT_YIELD, 6497);
        Assert.assertEquals(metrics.PF_Q20_EQUIVALENT_YIELD, 6497);
        Assert.assertEquals(metrics.READ_LENGTH_AVG_Q_ABOVE_30, 0);
        Assert.assertEquals(metrics.READ_LENGTH_AVG_Q_ABOVE_25, 0);
    }

    @Test
    public void testSubsample() throws IOException {
        final File input = new File(TEST_DATA_DIR, "subsample.bam");
        final File outfile   = File.createTempFile("test", ".quality_yield_metrics");
        outfile.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "USE_ORIGINAL_QUALITIES=false"
        };

        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<CollectQualityYieldMetrics.QualityYieldMetrics, ?> output = new MetricsFile<>();
        output.read(new FileReader(outfile));

        Assert.assertEquals(output.getMetrics().size(),1);

        final CollectQualityYieldMetrics.QualityYieldMetrics metrics = output.getMetrics().get(0);
        Assert.assertEquals(metrics.TOTAL_READS, 56);
        Assert.assertEquals(metrics.PF_READS, 56);
        Assert.assertEquals(metrics.READ_LENGTH, 285);
        Assert.assertEquals(metrics.TOTAL_BASES, 15983);
        Assert.assertEquals(metrics.PF_BASES, 15983);
        Assert.assertEquals(metrics.Q20_BASES, 15494);
        Assert.assertEquals(metrics.PF_Q20_BASES, 15494);
        Assert.assertEquals(metrics.Q30_BASES, 14786);
        Assert.assertEquals(metrics.PF_Q30_BASES, 14786);
        Assert.assertEquals(metrics.Q20_EQUIVALENT_YIELD, 30589);
        Assert.assertEquals(metrics.PF_Q20_EQUIVALENT_YIELD, 30589);
        Assert.assertEquals(metrics.READ_LENGTH_AVG_Q_ABOVE_30, 101);
        Assert.assertEquals(metrics.READ_LENGTH_AVG_Q_ABOVE_25, 195);
    }

    @Test
    public void testMulti() throws IOException {
        final File input = new File(TEST_DATA_DIR, "insert_size_metrics_test.sam");
        final File reference = new File(CollectAlignmentSummaryMetricsTest.TEST_DATA_DIR, "summary_alignment_stats_test.fasta");
        final File outfile   = File.createTempFile("test", ".quality_yield_metrics");
        outfile.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + reference.getAbsolutePath(),
                "METRIC_ACCUMULATION_LEVEL=" + MetricAccumulationLevel.ALL_READS.name(),
                "PROGRAM=null",
                "PROGRAM=" + CollectMultipleMetrics.Program.CollectQualityYieldMetrics.name()
        };

        Assert.assertEquals(runPicardCommandLine(CollectMultipleMetrics.class.getSimpleName(), args), 0);

        final MetricsFile<CollectQualityYieldMetrics.QualityYieldMetrics, ?> output = new MetricsFile<>();
        output.read(new FileReader(outfile + ".quality_yield_metrics"));

        Assert.assertEquals(output.getMetrics().size(),1);

        final CollectQualityYieldMetrics.QualityYieldMetrics metrics = output.getMetrics().get(0);
        Assert.assertEquals(metrics.TOTAL_READS, 52);
        Assert.assertEquals(metrics.PF_READS, 52);
        Assert.assertEquals(metrics.READ_LENGTH, 101);
        Assert.assertEquals(metrics.TOTAL_BASES, 5252);
        Assert.assertEquals(metrics.PF_BASES, 5252);
        Assert.assertEquals(metrics.Q20_BASES, 3532);
        Assert.assertEquals(metrics.PF_Q20_BASES, 3532);
        Assert.assertEquals(metrics.Q30_BASES, 3145);
        Assert.assertEquals(metrics.PF_Q30_BASES, 3145);
        Assert.assertEquals(metrics.Q20_EQUIVALENT_YIELD, 6497);
        Assert.assertEquals(metrics.PF_Q20_EQUIVALENT_YIELD, 6497);
        Assert.assertEquals(metrics.READ_LENGTH_AVG_Q_ABOVE_30, 0);
        Assert.assertEquals(metrics.READ_LENGTH_AVG_Q_ABOVE_25, 0);
    }

    @Test
    public void testMerge() {
        CollectQualityYieldMetrics.QualityYieldMetrics      m1 = createTestQualityYieldMetrics();
        CollectQualityYieldMetrics.QualityYieldMetrics      m2 = createTestQualityYieldMetrics();

        // create test data for averaging fields
        m2.READ_LENGTH_AVG_Q_ABOVE_30 *= 3;
        m2.READ_LENGTH_AVG_Q_ABOVE_25 *= 3;

        m1.merge(m2);

        Assert.assertEquals(m1.TOTAL_READS, m2.TOTAL_READS * 2);
        Assert.assertEquals(m1.PF_READS, m2.PF_READS * 2);
        Assert.assertEquals(m1.READ_LENGTH, m2.READ_LENGTH);
        Assert.assertEquals(m1.PF_BASES, m2.PF_BASES * 2);
        Assert.assertEquals(m1.Q20_BASES, m2.Q20_BASES * 2);
        Assert.assertEquals(m1.PF_Q20_BASES, m2.PF_Q20_BASES * 2);
        Assert.assertEquals(m1.Q30_BASES, m2.Q30_BASES * 2);
        Assert.assertEquals(m1.PF_Q30_BASES, m2.PF_Q30_BASES * 2);
        Assert.assertEquals(m1.Q20_EQUIVALENT_YIELD, m2.Q20_EQUIVALENT_YIELD * 2);
        Assert.assertEquals(m1.READ_LENGTH_AVG_Q_ABOVE_30, m2.READ_LENGTH_AVG_Q_ABOVE_30 * 2 / 3);
        Assert.assertEquals(m1.READ_LENGTH_AVG_Q_ABOVE_25, m2.READ_LENGTH_AVG_Q_ABOVE_25 * 2 / 3);
    }

    private CollectQualityYieldMetrics.QualityYieldMetrics createTestQualityYieldMetrics() {

        CollectQualityYieldMetrics.QualityYieldMetrics      m = new CollectQualityYieldMetrics.QualityYieldMetrics();

        m.TOTAL_READS = 52;
        m.PF_READS = 52;
        m.READ_LENGTH = 101;
        m.TOTAL_BASES = 5252;
        m.PF_BASES = 5252;
        m.Q20_BASES = 3532;
        m.PF_Q20_BASES = 3532;
        m.Q30_BASES = 3145;
        m.PF_Q30_BASES = 3145;
        m.Q20_EQUIVALENT_YIELD = 6497;
        m.PF_Q20_EQUIVALENT_YIELD = 6497;
        m.READ_LENGTH_AVG_Q_ABOVE_30 = 100;
        m.READ_LENGTH_AVG_Q_ABOVE_25 = 200;

        return m;
    }


}
