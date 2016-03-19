package com.astrazeneca.seq2c;

import org.apache.commons.math3.analysis.function.Log;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.commons.math3.util.Precision;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

/**
 * Normalize the coverage from targeted sequencing to CNV log2 ratio. The algorithm assumes the medium is diploid, thus not suitable
 * for homogeneous samples (e.g. parent-child).
 */

public class Cov2lr {
    final static int BUFFER_SIZE = 1000;
    private final String coverageFile;

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
    private BlockingQueue samples;

    /**
     * The failed factor for individual amplicons. If (the 80th percentile of an amplicon depth)/(the global median depth) is less
     * than the argument, the amplicon is considered failed and won't be used in calculation. Default: 0.2.
     */
    private final double FAILEDFACTOR = 0.2;

    /**
     * Constructor reads the input files and constructs genes, samples, factor maps
     *
     * @param amplicon       = determines the aggregation level (gene or record)
     * @param controlSamples = multiple controls are allowed, which are separated by ":"
     */

    public Cov2lr(boolean amplicon, Map<String, Long> stat, String covFile, String controlSamples) {
        this.coverageFile = covFile;
        init(amplicon, stat);
        initFactor(controlSamples);
    }

    private void init(boolean amplicon, Map<String, Long> stat) {
        this.amplicon = amplicon;
        this.mappingReads = stat;
        BlockingQueue<Sample> queue = new LinkedBlockingQueue<>(BUFFER_SIZE);
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
                    String tag = "";
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

    public ArrayList<Sample> doWork() {

        Set<String> samples = new HashSet<>();

        Map<String, Double> norm1Map = new HashMap<>();
        double medDepth = getMedDepth(samples, norm1Map);

        List<Double> gooddepth = splitQualitySamples(samples, medDepth);

        medDepth = median.evaluate(toDoubleArray(gooddepth));

        Map<String, Double> factor2 = getFactor2(medDepth);

        Map<String, Double> sampleMedian = getSampleMedian(samples);

        setNorm(medDepth, factor2, sampleMedian);

        return printResult();
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

    public void doTest() {

        try (FileWriter writer = new FileWriter("java_debug.txt")) {
            writer.write("Factor:" + "\n");
            for (Map.Entry<String, Double> entry : factor.entrySet()) {
                writer.write(entry.getKey() + "\t" + entry.getValue() + "\n");
                writer.flush();
            }
            writer.write("*********************" + "\n");
            writer.write("Loc hash:" + "\n");
            writer.flush();
            long start = 0;
            long end = 0;
            long len = 0;


            writer.write("start = " + start + " ");
            writer.write("end = " + end + " ");
            writer.write("len = " + len + "\n");
            writer.write("*********************" + "\n");
            writer.write("Data hash:" + "\n");
            writer.flush();
            long n = samples.size();
            long lenS = 0;
            long cov = 0;
            for (Map.Entry<String, Map<String, Sample>> entry : samples.entrySet()) {
                for (Map.Entry<String, Sample> sampleEntry : entry.getValue().entrySet()) {
                    lenS += sampleEntry.getValue().getLen();
                    cov += sampleEntry.getValue().getCov();
                }
            }
            writer.write("data len = " + lenS + " ");
            writer.write("cov = " + cov + " ");
            writer.write("numK = " + n + "\n");
            writer.write("*********************" + "\n");

        } catch (IOException e) {
            System.err.println(e.getLocalizedMessage());
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

    private void setNorm(double medDepth, Map<String, Double> factor2, Map<String, Double> sampleMedian) {
        for (Map.Entry<String, Map<String, Sample>> entry : samples.entrySet()) {
            String key = entry.getKey();
            double fact2 = factor2.get(key);
            for (Sample sample : entry.getValue().values()) {
                double norm1 = sample.getNorm1();
                double smplMed = sampleMedian.get(sample.getSample());
                sample.setNorm1b(Precision.round(norm1 * fact2 + 0.1, 2));
                sample.setNorm2(Precision.round(medDepth != 0 ? log.value((norm1 * fact2 + 0.1) / medDepth) / LOG_OF_TWO : 0, 2));
                sample.setNorm3(Precision.round(smplMed != 0 ? log.value((norm1 * fact2 + 0.1) / smplMed) / LOG_OF_TWO : 0, 2));
            }
        }
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

    private Map<String, Double> getFactor2(final double medDepth) {
        Map<String, Double> factor2 = new LinkedHashMap<>();
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
            } else {
                factor2.put(key, 0.0);
            }

        }
        return factor2;
    }

    private List<Double> splitQualitySamples(Set<String> samp, double medDepth) {
        List<Double> gooddepth = new LinkedList<>();
        Set<String> bad = new HashSet<>();
        SampleIterator iterator = new SampleIterator(coverageFile, amplicon);
        Sample sample = iterator.nextSample();

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

    private double getMedDepth(Set<String> samp, Map<String, Double> norm1Map) {
        List<Double> depth = new ArrayList<>();

        SampleIterator iterator = new SampleIterator(coverageFile, amplicon);
        Sample sample = iterator.nextSample();
        while (sample != null) {
            double norm1 = Precision.round((sample.getCov() * factor.get(sample.getSample())), 2);
            norm1Map.put(sample.getKey(), norm1);
            depth.add(norm1);
            samp.add(sample.getSample());
            sample = iterator.nextSample();
        }
        return median.evaluate(toDoubleArray(depth));
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