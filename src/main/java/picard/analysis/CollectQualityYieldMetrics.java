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
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Command line program to calculate quality yield metrics
 *
 * @author Martha Borkan
 */


@CommandLineProgramProperties(
        summary = CollectQualityYieldMetrics.USAGE_SUMMARY + CollectQualityYieldMetrics.USAGE_DETAILS,
        oneLineSummary = CollectQualityYieldMetrics.USAGE_SUMMARY,
        programGroup = DiagnosticsAndQCProgramGroup.class
)
@DocumentedFeature
public class CollectQualityYieldMetrics extends SinglePassSamProgram {
    private QualityYieldMetricsCollector collector = null;

    static final String USAGE_SUMMARY = "Collect metrics about reads that pass quality thresholds and Illumina-specific filters.  ";
    static final String USAGE_DETAILS = "This tool evaluates the overall quality of reads within a bam file containing one read group. " +
            "The output indicates the total numbers of bases within a read group that pass a minimum base quality score threshold and " +
            "(in the case of Illumina data) pass Illumina quality filters as described in the " +
            "<a href='https://www.broadinstitute.org/gatk/guide/article?id=6329'>GATK Dictionary entry</a>. " +
            "<br />" +
            "<h4>Note on base quality score options</h4>" +
            "If the quality score of read bases has been modified in a previous data processing step such as " +
            "<a href='https://www.broadinstitute.org/gatk/guide/article?id=44'>GATK Base Recalibration</a> " +
            "and an OQ tag is available, this tool can be set to use the OQ value instead of the primary quality value for the evaluation. " +
            "<br /><br />" +
            "Note that the default behaviour of this program changed as of November 6th 2015 to no longer include secondary and " +
            "supplemental alignments in the computation. <br />" +
            "<h4>Usage Example:</h4>" +
            "<pre>" +
            "java -jar picard.jar CollectQualityYieldMetrics \\<br /> " +
            "      I=input.bam \\<br /> " +
            "      O=quality_yield_metrics.txt \\<br />" +
            "</pre>" +
            "Please see " +
            "<a href='https://broadinstitute.github.io/picard/picard-metric-definitions.html#CollectQualityYieldMetrics.QualityYieldMetrics'>" +
            "the QualityYieldMetrics documentation</a> for details and explanations of the output metrics." +
            "<hr />";

    @Argument(shortName = StandardOptionDefinitions.USE_ORIGINAL_QUALITIES_SHORT_NAME,
            doc = "If available in the OQ tag, use the original quality scores " +
                    "as inputs instead of the quality scores in the QUAL field.")
    public boolean USE_ORIGINAL_QUALITIES = true;

    @Argument(doc = "If true, include bases from secondary alignments in metrics. Setting to true may cause double-counting " +
            "of bases if there are secondary alignments in the input file.")
    public boolean INCLUDE_SECONDARY_ALIGNMENTS = false;

    @Argument(doc = "If true, include bases from supplemental alignments in metrics. Setting to true may cause double-counting " +
            "of bases if there are supplemental alignments in the input file.")
    public boolean INCLUDE_SUPPLEMENTAL_ALIGNMENTS = false;

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
        this.collector = new QualityYieldMetricsCollector(USE_ORIGINAL_QUALITIES, INCLUDE_SECONDARY_ALIGNMENTS, INCLUDE_SUPPLEMENTAL_ALIGNMENTS);
    }

    @Override
    protected void acceptRead(final SAMRecord rec, final ReferenceSequence ref) {
        this.collector.acceptRecord(rec, ref);
    }

    @Override
    protected void finish() {
        final MetricsFile<QualityYieldMetrics, Integer> metricsFile = getMetricsFile();
        this.collector.finish();
        this.collector.addMetricsToFile(metricsFile);
        metricsFile.write(OUTPUT);
    }

    public static class QualityYieldMetricsCollector {
        // If true, include bases from secondary alignments in metrics. Setting to true may cause double-counting
        // of bases if there are secondary alignments in the input file.
        private final boolean useOriginalQualities;

        // If true, include bases from secondary alignments in metrics. Setting to true may cause double-counting
        // of bases if there are secondary alignments in the input file.
        private final boolean includeSecondaryAlignments;

        // If true, include bases from supplemental alignments in metrics. Setting to true may cause double-counting
        // of bases if there are supplemental alignments in the input file.
        public final boolean includeSupplementalAlignments;

        private final List<Double> qualityAccumulator = new ArrayList<>();
        private final List<Integer> qualityCount = new ArrayList<>();

        private double[] qual2prob = {1.00000000e+00, 7.94328235e-01, 6.30957344e-01, 5.01187234e-01,
                3.98107171e-01, 3.16227766e-01, 2.51188643e-01, 1.99526231e-01,
                1.58489319e-01, 1.25892541e-01, 1.00000000e-01, 7.94328235e-02,
                6.30957344e-02, 5.01187234e-02, 3.98107171e-02, 3.16227766e-02,
                2.51188643e-02, 1.99526231e-02, 1.58489319e-02, 1.25892541e-02,
                1.00000000e-02, 7.94328235e-03, 6.30957344e-03, 5.01187234e-03,
                3.98107171e-03, 3.16227766e-03, 2.51188643e-03, 1.99526231e-03,
                1.58489319e-03, 1.25892541e-03, 1.00000000e-03, 7.94328235e-04,
                6.30957344e-04, 5.01187234e-04, 3.98107171e-04, 3.16227766e-04,
                2.51188643e-04, 1.99526231e-04, 1.58489319e-04, 1.25892541e-04,
                1.00000000e-04, 7.94328235e-05, 6.30957344e-05, 5.01187234e-05,
                3.98107171e-05, 3.16227766e-05, 2.51188643e-05, 1.99526231e-05,
                1.58489319e-05, 1.25892541e-05, 1.00000000e-05, 7.94328235e-06,
                6.30957344e-06, 5.01187234e-06, 3.98107171e-06, 3.16227766e-06,
                2.51188643e-06, 1.99526231e-06, 1.58489319e-06, 1.25892541e-06};
        // The metrics to be accumulated
        private final QualityYieldMetrics metrics = new QualityYieldMetrics();

        public QualityYieldMetricsCollector(final boolean useOriginalQualities,
                                            final boolean includeSecondaryAlignments,
                                            final boolean includeSupplementalAlignments) {
            this.useOriginalQualities = useOriginalQualities;
            this.includeSecondaryAlignments = includeSecondaryAlignments;
            this.includeSupplementalAlignments = includeSupplementalAlignments;
        }

        public void acceptRecord(final SAMRecord rec, final ReferenceSequence ref) {
            if (!this.includeSecondaryAlignments && rec.isSecondaryAlignment()) return;
            if (!this.includeSupplementalAlignments && rec.getSupplementaryAlignmentFlag()) return;

            final int length = rec.getReadLength();
            metrics.TOTAL_READS++;
            metrics.TOTAL_BASES += length;

            final boolean isPfRead = !rec.getReadFailsVendorQualityCheckFlag();
            if (isPfRead) {
                metrics.PF_READS++;
                metrics.PF_BASES += length;
            }

            final byte[] quals;
            if (this.useOriginalQualities) {
                byte[] tmp = rec.getOriginalBaseQualities();
                if (tmp == null) tmp = rec.getBaseQualities();
                quals = tmp;
            } else {
                quals = rec.getBaseQualities();
            }

            // add up quals, and quals >= 20
            for (final int qual : quals) {
                metrics.Q20_EQUIVALENT_YIELD += qual;

                if (qual >= 30) {
                    metrics.Q20_BASES++;
                    metrics.Q30_BASES++;
                } else if (qual >= 20) {
                    metrics.Q20_BASES++;
                }

                if (isPfRead) {
                    metrics.PF_Q20_EQUIVALENT_YIELD += qual;
                    if (qual >= 30) {
                        metrics.PF_Q20_BASES++;
                        metrics.PF_Q30_BASES++;
                    } else if (qual >= 20) {
                        metrics.PF_Q20_BASES++;
                    }
                }
            }
            if (rec.getReadNegativeStrandFlag()) {
                reverseArray(quals);
            }
            int count = 0 ;
            for(final int qual: quals ) {
                Double prob = qual2prob[qual];
                if (count < qualityAccumulator.size()){
                    qualityAccumulator.set(count, qualityAccumulator.get(count) + prob);
                    qualityCount.set(count, qualityCount.get(count)+1);
                } else {
                    qualityAccumulator.add(count, prob);
                    qualityCount.add(count,1);
                }
                count ++;
            }
        }

        public void finish() {
            metrics.Q20_EQUIVALENT_YIELD = metrics.Q20_EQUIVALENT_YIELD / 20;
            metrics.PF_Q20_EQUIVALENT_YIELD = metrics.PF_Q20_EQUIVALENT_YIELD / 20;

            metrics.calculateDerivedFields();
            metrics.RLQ30 = calculateLQ(30);
            metrics.RLQ25 = calculateLQ(25);

        }

        public void addMetricsToFile(final MetricsFile<QualityYieldMetrics, Integer> metricsFile) {
            metricsFile.addMetric(metrics);
        }

        private int calculateLQ(final double threshold){
            double error_prob_threshold = Math.pow(10, -(double)threshold/10);
            final int OFFSET = 10;
            int cur_result = 0 ;
            List<Double> result = new ArrayList<>();
            for (int i = OFFSET; i < qualityAccumulator.size(); i++ ){
                if (qualityCount.get(i) < 50){
                    break;
                }
                result.add(qualityAccumulator.get(i)/qualityCount.get(i));
            }
            double cumsum = 0;
            for (int i = 0; i < result.size(); i++){
                cumsum += result.get(i);
                if (cumsum/(i+1) < error_prob_threshold) {
                    cur_result = i+1+OFFSET;
                }
            }
            return cur_result;
        }

        private void reverseArray(final byte[] array) {
            for (int i=0, j=array.length-1; i<j; ++i, --j) {
                final byte tmp = array[i];
                array[i] = array[j];
                array[j] = tmp;
            }
        }

    }

    /**
     * A set of metrics used to describe the general quality of a BAM file
     */
    public static class QualityYieldMetrics extends MergeableMetricBase {

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
         * The average read length of all the reads (will be fixed for a lane)
         */
        @NoMergingIsDerived
        public int READ_LENGTH = 0;

        /**
         * The total number of bases in all reads
         */
        @MergeByAdding
        public long TOTAL_BASES;

        /**
         * The total number of bases in all PF reads
         */
        @MergeByAdding
        public long PF_BASES = 0;

        /**
         * The number of bases in all reads that achieve quality score 20 or higher
         */
        @MergeByAdding
        public long Q20_BASES = 0;

        /**
         * The number of bases in PF reads that achieve quality score 20 or higher
         */
        @MergeByAdding
        public long PF_Q20_BASES = 0;

        /**
         * The number of bases in all reads that achieve quality score 30 or higher
         */
        @MergeByAdding
        public long Q30_BASES = 0;

        /**
         * The number of bases in PF reads that achieve quality score 30 or higher
         */
        @MergeByAdding
        public long PF_Q30_BASES = 0;

        /**
         * The sum of quality scores of all bases divided by 20
         */
        @MergeByAdding
        public long Q20_EQUIVALENT_YIELD = 0;

        /**
         * The sum of quality scores of all bases in PF reads divided by 20
         */
        @MergeByAdding
        public long PF_Q20_EQUIVALENT_YIELD = 0;

        @Override
        public void calculateDerivedFields() {
            super.calculateDerivedFields();

            this.READ_LENGTH = this.TOTAL_READS == 0 ? 0 : (int) (this.TOTAL_BASES / this.TOTAL_READS);
        }
        /** The average read length until the average base quality is above 30 */
        public long RLQ30 = 0;

        /** The average read length until the average base quality is above 25 */
        public long RLQ25 = 0;

    }

}
