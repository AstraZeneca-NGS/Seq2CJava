package com.astrazeneca.seq2c;

import org.apache.commons.math3.analysis.function.Log;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.commons.math3.util.Precision;

import java.io.*;
import java.util.*;

/**
 * Normalize the coverage from targeted sequencing to CNV log2 ratio. The algorithm assumes the medium is diploid, thus not suitable
 * for homogeneous samples (e.g. parent-child).
 */

public class Cov2lr {
    /**
     * Help utils for statistics
     */
    private Log log = new Log();
    private Median median = new Median();
    private Mean mean = new Mean();
    private final double LOG_OF_TWO = log.value(2);
    /**
     * Indicate this is amplicon or exon based calling. By default, it'll aggregate at gene level.
     */
    private boolean amplicon;

    /**
     * Array of control samples
     */
    private String[] controlSamples;

    /**
     * Statistics from first input file Key = name of sample Value = number of aligned reads
     */
    private Map<String, Long> mappingReads;

    /**
     * Normalization factor for each sample calculated as mean number of reads / number of reads per sample Key = name of sample
     * Value = factor
     */
    private Map<String, Double> factor;
    /**
     * Map of samples Key = key for string, consists of sample name, gene name, start, end and length Value = map of sample names
     * and sample objects
     */
    private Map<String, Map<String, Sample>> samples;

    /**
     * The failed factor for individual amplicons. If (the 80th percentile of an amplicon depth)/(the global median depth) is less
     * than the argument, the amplicon is considered failed and won't be used in calculation. Default: 0.2.
     */
    private final double FAILEDFACTOR = 0.2;

    private final String covFile;

    /**
     * Constructor reads the input files and constructs genes, samples, factor maps
     *
     * @param amplicon       = determines the aggregation level (gene or record)
     * @param controlSamples = multiple controls are allowed, which are separated by ":"
     */

//    public Cov2lr(boolean amplicon, Map<String, Long> stat, Collection<Gene> genesArray, String controlSamples) {
//        init(amplicon, stat);
////        readCoverage(genesArray, this.samples);
//        initFactor(controlSamples);
//    }
    public Cov2lr(boolean amplicon, Map<String, Long> stat, String covFile, String controlSamples) {
        init(amplicon, stat);
//        readCoverageFile(covFile, this.samples);
        initFactor(controlSamples);
        this.covFile = covFile;

    }

    private void init(boolean amplicon, Map<String, Long> stat) {
        this.amplicon = amplicon;
        this.mappingReads = stat;
        this.samples = new HashMap<>();
    }


    private void initFactor(String controlSamples) {
        this.factor = new LinkedHashMap<>();
        setFactor();
        if (controlSamples != null && !controlSamples.trim().isEmpty()) {
            this.controlSamples = controlSamples.split(":");
        }
    }

    private void readCoverageFile(String covFile, Map<String, Map<String, Sample>> samples) {
        try (BufferedReader reader = new BufferedReader(new FileReader(covFile))) {
            String line;
            while ((line = reader.readLine()) != null) {
                //skip lines with words "Sample, Whole, Control, Undertermined
                if (line.contains("Sample") || line.contains("Whole") || line.contains("Control") ||
                        line.contains("Undertermined")) continue;

                String[] sampleLines = line.split("\\s+");
                if (sampleLines.length == 8) {

                    String sample = sampleLines[0];
                    String gene = sampleLines[1];
                    String chr = sampleLines[2];
                    long start = Long.parseLong(sampleLines[3]);
                    long end = Long.parseLong(sampleLines[4]);
                    long len = Long.parseLong(sampleLines[6]);

                    double depth = Double.parseDouble(sampleLines[7]);

                    String key = amplicon ? join(" ", gene, chr, Long.toString(start), Long.toString(end), Long.toString(len)) : gene;
                    addSample(samples, key, sample, chr, start, end, len, Precision.round(depth, 2), gene);
                }
            }
        } catch (IOException e) {
            System.err.println("Cannot open file " + covFile);
        } catch (NumberFormatException e) {
            System.err.println("Cannot parse long number");
        }
    }

    /**
     * Main method, makes the statistics calculation according to the algorithm print result to the Standart output
     */

    public List<Sample> doWork() {

        Set<String> samp = new HashSet<>();
//        readCoverageFile(covFile, this.samples);

        long start, end;
        //      start = System.nanoTime();
        //       double medDepth = getMedDepth(samp);
        //       end = System.nanoTime();

//        System.out.println("expected " + medDepth);
//        System.out.println("elapsed " + (end - start) / 1_000_000);

//        start = System.nanoTime();
        Map<String, Map<String, Double>> norm1 = new LinkedHashMap<>();
        double medDepth = getMedDepthStream(samp, norm1);
        //       end = System.nanoTime();
//
//        System.out.println(samples.size());
//        System.out.println(norm1.size());
//
//        System.out.println("actual " + medDepth);
//        System.out.println("elapsed " + (end - start) / 1_000_000);

//        start = System.nanoTime();
//        List<Double> goodDepth = splitQualitySamples(samp, medDepth);
//        double medDepthGood = median.evaluate(toDoubleArray(goodDepth));
//        end = System.nanoTime();

//        System.out.println("expected " + medDepthGood);
//        System.out.println("elapsed " + (end - start) / 1_000_000);

//        start = System.nanoTime();
        List<Double> gooddepthStream = splitQualitySamplesStream(samp, medDepth, norm1);
        double medDepthGood = median.evaluate(toDoubleArray(gooddepthStream));
//        System.out.println("actual " + medDepthGood);
////        System.out.println("elapsed " + (end - start) / 1_000_000);
////
//        System.out.println(samples.size());
//        System.out.println(norm1.size());
//
//        Map<String, Double> factor2 = getFactor2(medDepthGood);
        Map<String, Double> factorStream = getFactor2Stream(medDepthGood, norm1);
//
//        Map<String, Double> sampleMedian = getSampleMedian(samp);
        Map<String, Double> sampleMedianStream = getSampleMedianStream(samp, norm1);

//        System.out.println("expected " + sampleMedian);
//        System.out.println("actual " + sampleMedianStream);
//
//        setNormStream(medDepth, factor2, sampleMedian, norm1);
        return setNormStream(medDepth, factorStream, sampleMedianStream, norm1);
    }

    private Map<String, Map<String, Double>> readCoverage() {
        SampleIterator iterator = new SampleIterator(covFile, amplicon);
        Map<String, Map<String, Double>> coverage = new LinkedHashMap<>();
        Sample sample = iterator.nextSample();
        while (sample != null) {
            String gene = sample.getName();
            Map<String, Double> genesToNorm = coverage.get(gene);
            if (genesToNorm == null) {
                genesToNorm = new HashMap<>();
            }
            genesToNorm.put(sample.getSample(), sample.getCov());
            coverage.put(gene, genesToNorm);
            sample = iterator.nextSample();
        }
        return coverage;
    }

    private Map<String, Long> readStat(String fileName) {
        Map<String, Long> map = new LinkedHashMap<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(fileName))) {
            String line;
            while ((line = reader.readLine()) != null) {
                String[] samples = line.split("\\s+");
                if (samples.length == 2) {
                    map.put(samples[0], Long.parseLong(samples[1]));
                }
            }
        } catch (IOException e) {
            System.err.println("Cannot open file " + fileName);
        } catch (NumberFormatException e) {
            System.err.println("Cannot parse long number");
        }
        return map;
    }

    private void readCoverage(Collection<Gene> genesArray, Map<String, Map<String, Sample>> samples) {
        for (Gene gn : genesArray) {
            if (gn.getTag().contains("Whole")) {
                continue;
            }
            // skip lines with words "Sample, Whole, Control, Undertermined
            String sample = gn.getSample();
            String gene = gn.getName();
            String chr = gn.getChr();
            String tag = gn.getTag();
            long start = gn.getStart();
            long end = gn.getEnd();
            long len = gn.getLen();
            double depth = gn.getDepth();
//            addGene(genes, sample, gene, chr, start, end, tag, len, Precision.round(depth, 2));
            String key = amplicon ? join(" ", gene, chr, Long.toString(start), Long.toString(end), Long.toString(len)) : gene;
            addSample(samples, key, sample, chr, start, end, len, Precision.round(depth, 2), gene);
        }
    }

    private void addGene(Map<String, Gene> map, String sample, String gene, String chr, long start, long end, String tag, long len, double depth) {

        if (map.containsKey(gene)) {
            Gene geneObj = map.get(gene);
            // move start
            if (geneObj.getStart() > start) {
                geneObj.setStart(start);
            }
            // move end
            if (geneObj.getEnd() < end) {
                geneObj.setEnd(end);
            }
            // update length
            geneObj.setLen(geneObj.getLen() + len);
        } else {
            Gene geneObj = new Gene(sample, gene, chr, start, end, tag, len, depth);
            map.put(gene, geneObj);
        }
    }

    private void addSample(Map<String, Map<String, Sample>> map, String key, String sample, String chr, Long start, long end, long len, double depth, String gene) {
        if (map.containsKey(key)) {
            Map<String, Sample> sampleMap = map.get(key);
            if (sampleMap.containsKey(sample)) {
                Sample sampleObj = sampleMap.get(sample);
                sampleObj.addLen(len);
                sampleObj.addCov(depth);
            } else {
                Sample sampleObj = new Sample(key, sample, chr, start, end, gene, len, depth);
                sampleMap.put(sample, sampleObj);
            }
        } else {
            Sample sampleObj = new Sample(key, sample, chr, start, end, gene, len, depth);
            Map<String, Sample> sampleMap = new LinkedHashMap<>();
            sampleMap.put(sample, sampleObj);
            map.put(key, sampleMap);
        }
    }

    private void setFactor() {
        double[] array = new double[mappingReads.size()];
        int i = 0;
        for (Long l : mappingReads.values()) {
            array[i++] = l;
        }
        double meanReads = mean.evaluate(array);

        for (Map.Entry<String, Long> entry : mappingReads.entrySet()) {
            factor.put(entry.getKey(), meanReads / entry.getValue());
        }

    }


    private ArrayList<Sample> printResult() {
        boolean useControlSamples = isUseControlSamples();
        ArrayList<Sample> sampleArr = new ArrayList<>();
        for (Map.Entry<String, Map<String, Sample>> entry : samples.entrySet()) {
            for (Map.Entry<String, Sample> entrySample : entry.getValue().entrySet()) {
                Sample sample = entrySample.getValue();
                if (useControlSamples) {
                    sample = addControlSamples(sample);
                }
                sampleArr.add(sample);
            }
        }
        return sampleArr;
    }

    public boolean isUseControlSamples() {
        return controlSamples != null && controlSamples.length > 0;
    }

    private Sample addControlSamples(Sample sample) {
        List<Double> list = new LinkedList<>();
        for (String s : controlSamples) {
            for (Map<String, Sample> entry : samples.values()) {
                Sample cs = entry.get(s);
                if (cs != null) {
                    list.add(cs.getNorm1b());
                }
            }
        }
        if (!list.isEmpty()) {
            double meanVal = mean.evaluate(toDoubleArray(list));

            sample.setNorm3s(meanVal == 0 ? sample.getNorm1b() / meanVal / log.value(2) : 0);
        }
        return sample;
    }

    private List<Sample> setNorm(double medDepth, Map<String, Double> factor2, Map<String, Double> sampleMedian) {
        List<Sample> sampleResult = new ArrayList<>();
        boolean useControlSamples = isUseControlSamples();

        for (Map.Entry<String, Map<String, Sample>> entry : samples.entrySet()) {
            String key = entry.getKey();
            double fact2 = factor2.get(key);
            for (Sample sample : entry.getValue().values()) {
                double norm1 = sample.getNorm1();
                double smplMed = sampleMedian.get(sample.getSample());
                sample.setNorm1b(Precision.round(norm1 * fact2 + 0.1, 2));
                sample.setNorm2(Precision.round(medDepth != 0 ? log.value((norm1 * fact2 + 0.1) / medDepth) / LOG_OF_TWO : 0, 2));
                sample.setNorm3(Precision.round(smplMed != 0 ? log.value((norm1 * fact2 + 0.1) / smplMed) / LOG_OF_TWO : 0, 2));
                if (useControlSamples) {
                    sample = addControlSamples(sample);
                }
                sampleResult.add(sample);
            }
        }
        return sampleResult;
    }

    private List<Sample> setNormStream(double medDepth, Map<String, Double> factor2, Map<String, Double> sampleMedian, Map<String, Map<String, Double>> norm) {
        List<Sample> sampleResult = new ArrayList<>();
        try {
            PrintWriter writer = new PrintWriter("temp.txt");
            boolean useControlSamples = isUseControlSamples();
            SampleIterator iterator = new SampleIterator(covFile, amplicon);
            Sample sample = iterator.nextSample();
            while (sample != null) {
                String key = sample.getName();
                if (!factor2.containsKey(key)) {
                    sample = iterator.nextSample();
                    continue;
                }
                double fact2 = factor2.get(key);
                double norm1 = norm.get(key).get(sample.getSample());
                double smplMed = sampleMedian.get(sample.getSample());
                sample.setNorm1b(Precision.round(norm1 * fact2 + 0.1, 2));
                sample.setNorm2(Precision.round(medDepth != 0 ? log.value((norm1 * fact2 + 0.1) / medDepth) / LOG_OF_TWO : 0, 2));
                sample.setNorm3(Precision.round(smplMed != 0 ? log.value((norm1 * fact2 + 0.1) / smplMed) / LOG_OF_TWO : 0, 2));
                if (useControlSamples) {
                    sample = addControlSamples(sample);
                }
//                sampleResult.add(sample);
                writer.write(sample.getResultString() + "\n");
                writer.flush();
                sample = iterator.nextSample();
            }


        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        return sampleResult;
    }

    private Map<String, Double> getSampleMedian(Set<String> samp) {
        Map<String, Double> sampleMedian = new LinkedHashMap<>();

        for (String s : samp) {
            List<Double> list = new LinkedList<>();
            for (Map.Entry<String, Map<String, Sample>> entry : samples.entrySet()) {
                list.add(entry.getValue().get(s).getNorm1());
            }
            sampleMedian.put(s, median.evaluate(toDoubleArray(list)));
        }
        return sampleMedian;
    }

    private Map<String, Double> getSampleMedianStream(Set<String> samp, Map<String, Map<String, Double>> norm) {
        Map<String, Double> sampleMedian = new LinkedHashMap<>();

        for (String s : samp) {
            List<Double> list = new LinkedList<>();
            for (Map.Entry<String, Map<String, Double>> entry : norm.entrySet()) {
                list.add(entry.getValue().get(s));
            }
            sampleMedian.put(s, median.evaluate(toDoubleArray(list)));
        }
        return sampleMedian;
    }

    private Map<String, Double> getFactor2(final double medDepth) {
        Map<String, Double> factor2 = new LinkedHashMap<>();
        double result = 0.0;
        for (Map.Entry<String, Map<String, Sample>> entry : samples.entrySet()) {
            String key = entry.getKey();
            double[] norms1 = new double[entry.getValue().size()];
            int idx = 0;
            for (Sample sample : entry.getValue().values()) {
                norms1[idx++] = sample.getNorm1();
            }
            double median = new Median().evaluate(norms1);
            if (median != 0) {
                factor2.put(key, medDepth / median);
                result += medDepth / median;
            } else {
                factor2.put(key, 0.0);
            }
        }
        System.out.println("expected " + result);
        System.out.println("expected size " + factor2.size());
        return factor2;
    }

    private Map<String, Double> getFactor2Stream(final double medDepth, Map<String, Map<String, Double>> norm) {
        Map<String, Double> factor2 = new HashMap<>();
        double result = 0.0;
        for (Map.Entry<String, Map<String, Double>> entry : norm.entrySet()) {
            String key = entry.getKey();
            double[] norms1 = new double[entry.getValue().size()];
            int idx = 0;
            for (Double normValue : entry.getValue().values()) {
                norms1[idx++] = normValue;
            }
            double median = new Median().evaluate(norms1);
            if (median != 0) {
                factor2.put(key, medDepth / median);
                result += medDepth / median;
            } else {
                factor2.put(key, 0.0);
            }
        }
        System.out.println("actual " + result);
        System.out.println("actual size " + factor2.size());
        return factor2;
    }

    private List<Double> splitQualitySamples(Set<String> samp, double medDepth) {
        List<Double> gooddepth = new LinkedList<>();
        Set<String> bad = new HashSet<>();
        for (Map.Entry<String, Map<String, Sample>> entry : samples.entrySet()) {
            List<Double> temp = new LinkedList<>();
            double kp80 = filterData(samp, entry.getValue(), temp);
            if (kp80 < medDepth * FAILEDFACTOR) {
                bad.add(entry.getKey());
            } else {
                gooddepth.addAll(temp);
            }
        }
        for (String s : bad) {
            samples.remove(s);
        }
        return gooddepth;
    }

    private List<Double> splitQualitySamplesStream(Set<String> samp, double medDepth, Map<String, Map<String, Double>> norm) {
        List<Double> gooddepth = new LinkedList<>();
        Set<String> bad = new HashSet<>();
        for (Map.Entry<String, Map<String, Double>> entry : norm.entrySet()) {
            List<Double> temp = new LinkedList<>();
            double kp80 = filterDataStream(samp, entry.getValue(), temp);
            if (kp80 < medDepth * FAILEDFACTOR) {
                bad.add(entry.getKey());
            } else {
                gooddepth.addAll(temp);
            }
        }
        for (String s : bad) {
            norm.remove(s);
        }
        return gooddepth;
    }

    private double getMedDepth(Set<String> samp) {
        List<Double> depth = new ArrayList<>();
        double result = 0.0;
        for (Map.Entry<String, Map<String, Sample>> entry : samples.entrySet()) {
            Map<String, Sample> sampleMap = entry.getValue();

            for (Map.Entry<String, Sample> sampleEntry : sampleMap.entrySet()) {
                Sample sample = sampleEntry.getValue();
                double norm1 = Precision.round((sample.getCov() * factor.get(sample.getSample())), 2);
                sample.setNorm1(norm1);
                depth.add(norm1);
                samp.add(sample.getSample());
                result += norm1;
            }
        }
//        System.out.println("expected sum " + result);
//        System.out.println("expected size " + depth.size());
        return median.evaluate(toDoubleArray(depth));
    }

    private double getMedDepthStream(Set<String> samp, Map<String, Map<String, Double>> norm) {
        List<Double> depth = new ArrayList<>();
        double result = 0.0;
        SampleIterator iterator = new SampleIterator(covFile, amplicon);
        Sample sample = iterator.nextSample();
        while (sample != null) {
            double norm1 = Precision.round((sample.getCov() * factor.get(sample.getSample())), 2);
            depth.add(norm1);
            putInMap(norm, sample, norm1);
            result += norm1;
            samp.add(sample.getSample());
            sample = iterator.nextSample();
        }
//        System.out.println("actual sum " + result);
//        System.out.println("actual size " + depth.size());
        return median.evaluate(toDoubleArray(depth));
    }


    private void putInMap(Map<String, Map<String, Double>> norm, Sample sample, double norm1) {
        String gene = sample.getName();
        Map<String, Double> genesToNorm = norm.get(gene);
        if (genesToNorm == null) {
            genesToNorm = new HashMap<>();
        }
        genesToNorm.put(sample.getSample(), norm1);
        norm.put(gene, genesToNorm);

    }

    private double filterData(Set<String> sampleSet, Map<String, Sample> map, List<Double> list) {
        for (String sample : sampleSet) {
            if (map.containsKey(sample)) {
                list.add(map.get(sample).getNorm1());
            }
        }
        double[] result = toDoubleArray(list);
        Percentile perc = new Percentile(80).withEstimationType(Percentile.EstimationType.R_5);
        return perc.evaluate(result);
    }

    private double filterDataStream(Set<String> sampleSet, Map<String, Double> map, List<Double> list) {
        for (String sample : sampleSet) {
            if (map.containsKey(sample)) {
                list.add(map.get(sample));
            }
        }
        double[] result = toDoubleArray(list);
        Percentile perc = new Percentile(80).withEstimationType(Percentile.EstimationType.R_5);
        return perc.evaluate(result);
    }

    private double[] toDoubleArray(List<Double> list) {
        double[] array = new double[list.size()];
        int i = 0;
        for (Double d : list) {
            array[i++] = d;
        }
        return array;
    }

    private static String join(String delim, Object... array) {
        StringBuilder builder = new StringBuilder();
        for (Object s : array) {
            if (builder.length() > 0) {
                builder.append(delim);
            }
            builder.append(s);
        }
        return builder.toString();
    }

}