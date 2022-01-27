/*
 * The MIT License
 *
 * Copyright (c) 2022 The Broad Institute
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

    public int calculateLQ(final int threshold, int readInPair, int spanningWindowLength){
        double errorProbThreshold = QualityUtil.getErrorProbabilityFromPhredScore(threshold);
        final int OFFSET = 10;
        int cur_result = 1 ;
        List<Double> result = new ArrayList<>();
        List<Long> weights = new ArrayList<>();
        double[] accumulator;
        long[] counts;
        if (readInPair == 1){
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
            weights.add(counts[i]);
        }
        applySpanningWindowMean(result, weights, spanningWindowLength);
        return longestHighQuality(result,errorProbThreshold);
    }

    private void applySpanningWindowMean(List<Double> vector, List<Long> weights, final int spanLength){
        List<Double> tmp = new ArrayList<>(vector);
        for (int i = 0; i < vector.size(); i++){
            double tmpEr = 0;
            long tmpWeight = 0;
            for (int j = Math.max(i-spanLength,0); j < Math.min(i+spanLength+1, vector.size()); j++){
                tmpEr += tmp.get(j)*weights.get(j);
                tmpWeight += weights.get(j);
            }
            vector.set(i, tmpEr/tmpWeight);
        }
    }

    private int longestHighQuality(List<Double> averageErrorProbabilities, double errorProbThreshold){
        int curStart = 0;
        int curEnd = 0;
        int curBestIntervalLength = 0;

        while ( curEnd < averageErrorProbabilities.size() ) {
            if (averageErrorProbabilities.get(curEnd) <= errorProbThreshold) {
                curEnd++;
            } else {
                if ((curEnd - curStart) > curBestIntervalLength) {
                    curBestIntervalLength = curEnd - curStart;
                }
                curStart = curEnd + 1;
                curEnd = curStart;
            }
        }
        if ((curEnd - curStart) > curBestIntervalLength) {
            curBestIntervalLength = curEnd - curStart;
        }
        return curBestIntervalLength;
    }
}
