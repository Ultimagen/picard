/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;
import picard.flow.FlowBasedArgumentCollection;
import picard.flow.FlowBasedRead;
import picard.flow.FlowBasedReadUtils;
import picard.util.help.HelpConstants;

import java.io.File;
import java.util.Arrays;

/**
 * Command line program to calculate quality yield metrics for flow based read files
 *
 * @author Dror Kessler
 */


@CommandLineProgramProperties(
        summary = CollectQualityYieldMetricsFlow.USAGE_SUMMARY + CollectQualityYieldMetricsFlow.USAGE_DETAILS,
        oneLineSummary = CollectQualityYieldMetricsFlow.USAGE_SUMMARY,
        programGroup = DiagnosticsAndQCProgramGroup.class
)
@ExperimentalFeature
public class CollectQualityYieldMetricsFlow extends SinglePassSamProgram {
    private QualityYieldMetricsCollectorFlow collector = null;
    public Histogram<Integer> qualityHistogram = new Histogram<>("KEY", "QUAL_COUNT");
    public Histogram<Integer> flowAvgQualityHistogram = new Histogram<>("KEY", "AVG_FLOW_QUAL");
    public int[] flowAvgQualityCount = null;
    public long[] flowAvgQualitySum = null;

    static final String USAGE_SUMMARY = "Collect metrics about reads that pass quality thresholds from flow based read files.  ";
    static final String USAGE_DETAILS = "This tool evaluates the overall quality of reads within a bam file containing one read group. " +
            "The output indicates the total numbers of flows within a read group that pass a minimum base quality score threshold " +
            "<h4>Usage Example:</h4>" +
            "<pre>" +
            "java -jar picard.jar CollectQualityYieldMetricsFlow \\<br /> " +
            "      I=input.bam \\<br /> " +
            "      O=quality_yield_metrics.txt \\<br />" +
            "</pre>" +
            "<hr />";

    @Argument(doc = "If true, include bases from secondary alignments in metrics. Setting to true may cause double-counting " +
            "of bases if there are secondary alignments in the input file.")
    public boolean INCLUDE_SECONDARY_ALIGNMENTS = false;

    @Argument(doc = "If true, include bases from supplemental alignments in metrics. Setting to true may cause double-counting " +
            "of bases if there are supplemental alignments in the input file.")
    public boolean INCLUDE_SUPPLEMENTAL_ALIGNMENTS = false;

    @ArgumentCollection(doc = "flow based args")
    public FlowBasedArgumentCollection fbargs = new FlowBasedArgumentCollection();

    /**
     * Ensure that we get all reads regardless of alignment status.
     */
    @Override
    protected boolean usesNoRefReads() {
        return true;
    }

    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {
        IOUtil.assertFileIsWritable(OUTPUT);
        this.collector = new QualityYieldMetricsCollectorFlow(INCLUDE_SECONDARY_ALIGNMENTS, INCLUDE_SUPPLEMENTAL_ALIGNMENTS, fbargs);
    }

    @Override
    protected void acceptRead(final SAMRecord rec, final ReferenceSequence ref) {
        this.collector.acceptRecord(rec, ref);
    }

    @Override
    protected void finish() {
        final MetricsFile<QualityYieldMetricsFlow, Integer> metricsFile = getMetricsFile();
        this.collector.finish();
        this.collector.addMetricsToFile(metricsFile);
        metricsFile.addHistogram(qualityHistogram);
        metricsFile.addHistogram(flowAvgQualityHistogram);
        metricsFile.write(OUTPUT);
    }

    public class QualityYieldMetricsCollectorFlow {
        // If true, include bases from secondary alignments in metrics. Setting to true may cause double-counting
        // of bases if there are secondary alignments in the input file.
        private final boolean includeSecondaryAlignments;

        // If true, include bases from supplemental alignments in metrics. Setting to true may cause double-counting
        // of bases if there are supplemental alignments in the input file.
        public final boolean includeSupplementalAlignments;

        // The metrics to be accumulated
        private final QualityYieldMetricsFlow metrics;

        // flow based args
        private final FlowBasedArgumentCollection fbargs;

        public QualityYieldMetricsCollectorFlow(final boolean includeSecondaryAlignments,
                                                final boolean includeSupplementalAlignments,
                                                final FlowBasedArgumentCollection fbargs) {
            this.includeSecondaryAlignments = includeSecondaryAlignments;
            this.includeSupplementalAlignments = includeSupplementalAlignments;
            this.fbargs = fbargs;
            this.metrics = new QualityYieldMetricsFlow();
        }

        public void acceptRecord(final SAMRecord rec, final ReferenceSequence ref) {
            if (!this.includeSecondaryAlignments && rec.isSecondaryAlignment()) return;
            if (!this.includeSupplementalAlignments && rec.getSupplementaryAlignmentFlag()) return;

            // convert to a flow based read
            FlowBasedReadUtils.ReadGroupInfo info = FlowBasedReadUtils.getReadGroupInfo(rec.getHeader(), rec);
            FlowBasedRead fread = new FlowBasedRead(rec, info.flowOrder, info.maxClass, fbargs);

            final int length = rec.getReadLength();
            metrics.TOTAL_READS++;
            metrics.TOTAL_FLOWS += fread.getKey().length;

            final boolean isPfRead = !rec.getReadFailsVendorQualityCheckFlag();
            if (isPfRead) {
                metrics.PF_READS++;
                metrics.PF_FLOWS += fread.getKey().length;
            }

            // get flow quals
            final byte[] quals = getFlowQualities(fread);

            // make sure flowQualityCount/Sum are large enough
            if ( flowAvgQualityCount == null ) {
                flowAvgQualityCount = new int[quals.length];
                flowAvgQualitySum = new long[quals.length];
            } else if ( flowAvgQualityCount.length < quals.length ) {
                flowAvgQualityCount = Arrays.copyOf(flowAvgQualityCount, quals.length);
                flowAvgQualitySum = Arrays.copyOf(flowAvgQualitySum, quals.length);
            }

            // add up quals, and quals >= 20
            int flow = 0;
            for (final int qual : quals) {
                metrics.Q20_EQUIVALENT_YIELD += qual;

                if (qual >= 30) {
                    metrics.Q20_FLOWS++;
                    metrics.Q30_FLOWS++;
                } else if (qual >= 20) {
                    metrics.Q20_FLOWS++;
                }

                if (isPfRead) {
                    metrics.PF_Q20_EQUIVALENT_YIELD += qual;
                    if (qual >= 30) {
                        metrics.PF_Q20_FLOWS++;
                        metrics.PF_Q30_FLOWS++;
                    } else if (qual >= 20) {
                        metrics.PF_Q20_FLOWS++;
                    }
                }

                // enter quality into histograms
                qualityHistogram.increment(qual);
                flowAvgQualityCount[flow]++;
                flowAvgQualitySum[flow] += qual;

                // advance
                flow++;
            }

        }

        private byte[] getFlowQualities(FlowBasedRead fread) {

            double[] errorProbs = computeErrorProb(fread);
            byte[] quals = new byte[errorProbs.length];
            for ( int i = 0 ; i < errorProbs.length ; i++ ) {
                quals[i] = (byte)Math.ceil(-10 * Math.log10(errorProbs[i]));
            }
            return quals;
        }

        /*
         * compute error probability vector for a read
         *
         * The vector has one element for each flow key, representing the probability complementing the call-probability to 1
         * This is further complicated by the optional presence of a genome-prior database, which provides factoring for
         * each hmer length (on a base basis)
         */
        private double[] computeErrorProb(final FlowBasedRead flowRead) {

            final int[] key = flowRead.getKey();
            final byte[] flowOrder = flowRead.getFlowOrderArray();

            final double[] probCol = new double[flowRead.getMaxHmer() + 1];
            double[] result = new double[key.length];

            for ( int i = 0 ; i < key.length ; i++ ) {

                // step 1 - extract column & sum
                double  sum = 0;
                for ( int j = 0 ; j < probCol.length ; j++ ) {
                    sum += (probCol[j] = flowRead.getProb(i, j));
                }

                // step 2 - normalize column
                for ( int j = 0 ; j < probCol.length ; j++ ) {
                    probCol[j] /= sum;
                }

                // assign normalized result
                result[i] = 1 - probCol[Math.min(key[i], flowRead.getMaxHmer())];
            }

            return result;
        }

        public void finish() {
            metrics.Q20_EQUIVALENT_YIELD = metrics.Q20_EQUIVALENT_YIELD / 20;
            metrics.PF_Q20_EQUIVALENT_YIELD = metrics.PF_Q20_EQUIVALENT_YIELD / 20;
            metrics.calculateDerivedFields();

            for (int i = 0; i < flowAvgQualityCount.length ; i++ ) {
                flowAvgQualityHistogram.increment(i, flowAvgQualityCount[i] != 0 ? ((double) flowAvgQualitySum[i] / flowAvgQualityCount[i]) : 0.0);
            }
        }

        public void addMetricsToFile(final MetricsFile<QualityYieldMetricsFlow, Integer> metricsFile) {
            metricsFile.addMetric(metrics);
        }
    }

    /**
     * A set of metrics used to describe the general quality of a BAM file
     */
    @DocumentedFeature(groupName = HelpConstants.DOC_CAT_METRICS, summary = HelpConstants.DOC_CAT_METRICS_SUMMARY)
    public static class QualityYieldMetricsFlow extends MergeableMetricBase {

        public QualityYieldMetricsFlow() {
        }

        /**
         * The total number of reads in the input file
         */
        @MergeByAdding
        public long TOTAL_READS = 0;

        /**
         * The number of reads that are PF - pass filter
         */
        @MergeByAdding
        public long PF_READS = 0;

        /**
         * The average read length of all the reads
         */
        @NoMergingIsDerived
        public int READ_LENGTH = 0;

        /**
         * The total number of flows in all reads
         */
        @MergeByAdding
        public long TOTAL_FLOWS;

        /**
         * The total number of flows in all PF reads
         */
        @MergeByAdding
        public long PF_FLOWS = 0;

        /**
         * The number of flows in all reads that achieve quality score 20 or higher
         */
        @MergeByAdding
        public long Q20_FLOWS = 0;

        /**
         * The number of flows in PF reads that achieve quality score 20 or higher
         */
        @MergeByAdding
        public long PF_Q20_FLOWS = 0;

        /**
         * The number of flows in all reads that achieve quality score 30 or higher
         */
        @MergeByAdding
        public long Q30_FLOWS = 0;

        /**
         * The number of flows in PF reads that achieve quality score 30 or higher
         */
        @MergeByAdding
        public long PF_Q30_FLOWS = 0;

        /**
         * The sum of quality scores of all flows divided by 20
         */
        @MergeByAdding
        public long Q20_EQUIVALENT_YIELD = 0;

        /**
         * The sum of quality scores of all flows in PF reads divided by 20
         */
        @MergeByAdding
        public long PF_Q20_EQUIVALENT_YIELD = 0;

        @Override
        public void calculateDerivedFields() {
            super.calculateDerivedFields();
            this.READ_LENGTH = this.TOTAL_READS == 0 ? 0 : (int) (this.TOTAL_FLOWS / this.TOTAL_READS);
        }

        @Override
        public MergeableMetricBase merge(final MergeableMetricBase other) {
            if (!(other instanceof QualityYieldMetricsFlow)){
                throw new PicardException("Only objects of the same type can be merged");
            }

            final QualityYieldMetricsFlow otherMetric = (QualityYieldMetricsFlow) other;

            super.merge(otherMetric);
            calculateDerivedFields();
            return this;
        }
    }
}
