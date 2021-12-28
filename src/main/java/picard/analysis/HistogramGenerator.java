package picard.analysis;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.QualityUtil;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public final class HistogramGenerator {
    final boolean useOriginalQualities;
    int maxLengthSoFar = 0;
    double[] firstReadTotalsByCycle  = new double[maxLengthSoFar];
    double[] firstReadTotalProbsByCycle  = new double[maxLengthSoFar];
    long[]   firstReadCountsByCycle  = new long[maxLengthSoFar];
    double[] secondReadTotalsByCycle = new double[maxLengthSoFar];
    double[] secondReadTotalProbsByCycle = new double[maxLengthSoFar];
    long[]   secondReadCountsByCycle = new long[maxLengthSoFar];

    public HistogramGenerator(final boolean useOriginalQualities) {
        this.useOriginalQualities = useOriginalQualities;
    }

    public void addRecord(final SAMRecord rec) {
        final byte[] quals = (useOriginalQualities ? rec.getOriginalBaseQualities() : rec.getBaseQualities());
        if (quals == null) return;

        final int length = quals.length;
        final boolean rc = rec.getReadNegativeStrandFlag();
        ensureArraysBigEnough(length+1);

        for (int i=0; i<length; ++i) {
            final int cycle = rc ? length-i : i+1;

            if (rec.getReadPairedFlag() && rec.getSecondOfPairFlag()) {
                secondReadTotalsByCycle[cycle] += quals[i];
                secondReadTotalProbsByCycle[cycle] += QualityUtil.getErrorProbabilityFromPhredScore(quals[i]);
                secondReadCountsByCycle[cycle] += 1;
            }
            else {
                firstReadTotalsByCycle[cycle] += quals[i];
                firstReadTotalProbsByCycle[cycle] += QualityUtil.getErrorProbabilityFromPhredScore(quals[i]);
                firstReadCountsByCycle[cycle] += 1;
            }
        }
    }

    private void ensureArraysBigEnough(final int length) {
        if (length > maxLengthSoFar) {
            firstReadTotalsByCycle  = Arrays.copyOf(firstReadTotalsByCycle, length);
            firstReadTotalProbsByCycle  = Arrays.copyOf(firstReadTotalProbsByCycle, length);
            firstReadCountsByCycle  = Arrays.copyOf(firstReadCountsByCycle, length);
            secondReadTotalsByCycle = Arrays.copyOf(secondReadTotalsByCycle , length);
            secondReadTotalProbsByCycle = Arrays.copyOf(secondReadTotalProbsByCycle , length);
            secondReadCountsByCycle = Arrays.copyOf(secondReadCountsByCycle, length);
            maxLengthSoFar = length;
        }
    }

    Histogram<Integer> getMeanQualityHistogram() {
        final String label = useOriginalQualities ? "MEAN_ORIGINAL_QUALITY" : "MEAN_QUALITY";
        final Histogram<Integer> meanQualities = new Histogram<Integer>("CYCLE", label);

        int firstReadLength = 0;

        for (int cycle=0; cycle < firstReadTotalsByCycle.length; ++cycle) {
            if (firstReadTotalsByCycle[cycle] > 0) {
                meanQualities.increment(cycle, firstReadTotalsByCycle[cycle] / firstReadCountsByCycle[cycle]);
                firstReadLength = cycle;
            }
        }

        for (int i=0; i< secondReadTotalsByCycle.length; ++i) {
            if (secondReadCountsByCycle[i] > 0) {
                final int cycle = firstReadLength + i;
                meanQualities.increment(cycle, secondReadTotalsByCycle[i] / secondReadCountsByCycle[i]);
            }
        }

        return meanQualities;
    }

    boolean isEmpty() {
        return maxLengthSoFar == 0;
    }

    public int calculateLQ(final int threshold, int read_in_pair){
        double error_prob_threshold = QualityUtil.getErrorProbabilityFromPhredScore(threshold);
        final int OFFSET = 10;
        int cur_result = 1 ;
        List<Double> result = new ArrayList<>();
        double[] accumulator;
        long[] counts;
        if (read_in_pair == 1){
            accumulator = firstReadTotalProbsByCycle;
            counts = firstReadCountsByCycle;
        } else {
            accumulator = secondReadTotalProbsByCycle;
            counts = secondReadCountsByCycle;
        }
        for (int i = OFFSET+1; i < accumulator.length; i++ ){
            if (counts[i] < 25){
                break;
            }
            result.add(accumulator[i]/counts[i]);
        }

        double cumsum = 0;
        for (int i = 0; i < result.size(); i++){
            cumsum += result.get(i);
            if (cumsum/(i+1) < error_prob_threshold) {
                cur_result = i+1+OFFSET;
            }
        }
        return cur_result-1;
    }

}
