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
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;

/**
 * Created by dror27 on 21/12/23.
 */
public class CollectQualityYieldMetricsFlowTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/sam/");

    public String getCommandLineProgramName() {
        return CollectQualityYieldMetricsFlow.class.getSimpleName();
    }

    @Test
    public void test() throws IOException {
        final File input = new File(TEST_DATA_DIR, "insert_size_metrics_test_flow.sam");
        final File outfile   = File.createTempFile("test", ".quality_yield_metrics_flow_space");
        outfile.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "INCLUDE_BQ_HISTOGRAM=TRUE",
                "--flow-fill-empty-bins-value", "0.000005",
        };

        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<CollectQualityYieldMetricsFlow.QualityYieldMetricsFlow, ?> output = new MetricsFile<>();
        output.read(new FileReader(outfile));

        Assert.assertEquals(output.getMetrics().size(),1);

        final CollectQualityYieldMetricsFlow.QualityYieldMetricsFlow metrics = output.getMetrics().get(0);
        Assert.assertEquals(metrics.TOTAL_READS, 56);
        Assert.assertEquals(metrics.PF_READS, 56);
        Assert.assertEquals(metrics.READ_LENGTH, 375);
        Assert.assertEquals(metrics.TOTAL_FLOWS, 21053);
        Assert.assertEquals(metrics.PF_FLOWS, 21053);
        Assert.assertEquals(metrics.Q20_FLOWS, 20671);
        Assert.assertEquals(metrics.PF_Q20_FLOWS, 20671);
        Assert.assertEquals(metrics.Q30_FLOWS, 20256);
        Assert.assertEquals(metrics.PF_Q30_FLOWS, 20256);
        Assert.assertEquals(metrics.Q20_EQUIVALENT_YIELD, 55247);
        Assert.assertEquals(metrics.PF_Q20_EQUIVALENT_YIELD, 55247);
    }
}
