package com.astrazeneca.seq2c;

import com.astrazeneca.seq2c.input.FileDataIterator;
import com.astrazeneca.seq2c.input.FileStoredDataFactory;
import com.astrazeneca.seq2c.input.Sample;
import com.astrazeneca.seq2c.input.SampleFactory;
import htsjdk.samtools.util.CloseableIterator;
import org.apache.commons.math3.analysis.function.Log;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.commons.math3.util.Precision;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.*;

/**
 * Normalize the coverage from targeted sequencing to CNV log2 ratio. The algorithm assumes the medium is diploid, thus not suitable
 * for homogeneous genes (e.g. parent-child).
 */
//perl version: cov2lr.pl
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
     * Array of control genes
     */
    private String[] controlSamples;

    /**
     * Statistics from first input file Key = name of sample Value = number of aligned reads
     */
    private Map<String, Long> mappingReads;

    /**
     * SampleStatistics factor for each sample calculated as mean number of reads / number of reads per sample Key = name of sample
     * Value = factor
     */
    private Map<String, Double> factor;

    /**
     * key - {@link Sample#name}
     * value - map: key - sample, value - norm1
     */
    private Map<String, Map<String, Double>> geneStatistics;
    /**
     * The failed factor for individual amplicons. If (the 80th percentile of an amplicon depth)/(the global median depth) is less
     * than the argument, the amplicon is considered failed and won't be used in calculation. Default: 0.2.
     */
    private final double FAILEDFACTOR = 0.2;

    private final String covFile;
    private final String tempFile;

    private final Map<String, Locus> locusMap;
    private FileStoredDataFactory<Sample> factory;

    /**
     * Constructor reads the input files and constructs genes, genes, factor maps
     *
     * @param amplicon       = determines the aggregation level (gene or record)
     * @param controlSamples = multiple controls are allowed, which are separated by ":"
     */

    public Cov2lr(boolean amplicon, Map<String, Long> stat, String covFile, String controlSamples, String tmpFile) {
        init(amplicon, stat);
        initFactor(controlSamples);
        this.covFile = covFile;
        this.tempFile = tmpFile;
        this.locusMap = new HashMap<>();
        this.factory = new SampleFactory(amplicon);
    }

    private void init(boolean amplicon, Map<String, Long> stat) {
        this.amplicon = amplicon;
        this.mappingReads = stat;
        this.geneStatistics = new HashMap<>();
    }

    private void initFactor(String controlSamples) {
        this.factor = new LinkedHashMap<>();
        setFactor();
        if (controlSamples != null && !controlSamples.trim().isEmpty()) {
            this.controlSamples = controlSamples.split(":");
        }
    }

    /**
     * Main method, makes the statistics calculation according to the algorithm print result to the Standart output
     */
    public void doWork() {

        Set<String> sampleNames = new HashSet<>();
        //perl version: 63 str
        double medDepth = readCoverageAndGetMedDepth(sampleNames);

        Set<String> badGenes = new HashSet<>();
        //perl version: 76 str
        List<Double> gooddepth = splitQualitySamples(sampleNames, medDepth, badGenes);

        //Re-adjust median depth using only those from good amplicons/genes
        //perl version: 87 str
        double medDepthGood = median.evaluate(toDoubleArray(gooddepth));
        //Gene/amplicon factor
        //perl version: 89 str
        Map<String, Double> factor = getFactor2(medDepthGood, badGenes);
        //perl version: 96 str
        Map<String, Double> sampleMedian = getSampleMedian(sampleNames, badGenes);
        //perl version: 107 str
        double controlSamplesMean = prepareControlSamples(factor);
        //perl version: 116 str
        setNormalization(medDepth, factor, sampleMedian, controlSamplesMean, badGenes);

        this.geneStatistics.clear();
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

    public boolean useControlSamples() {
        return controlSamples != null && controlSamples.length > 0;
    }

    /**
     * Evaluate mean value from control samples
     * @param factor2 factor from {@link Cov2lr#getFactor2(double, Set)}
     * @return mean quantity of (norm1 * fact2 + 0.1)
     */
    //perl version: 107 str
    double prepareControlSamples(Map<String, Double> factor2) {
        if (!useControlSamples()) return 0.0;
        List<Double> list = new LinkedList<>();
        for (String s : controlSamples) {
            for (Map.Entry<String, Map<String, Double>> entry : geneStatistics.entrySet()) {

                Double norm1 = entry.getValue().get(s);
                if (norm1 != null) {
                    if (factor2.containsKey(entry.getKey())) {
                        double fact2 = factor2.get(entry.getKey());
                        list.add(Precision.round(norm1 * fact2 + 0.1, 2));
                    }
                }
            }
        }
        double meanVal = 0.0;
        if (!list.isEmpty()) {
            meanVal = mean.evaluate(toDoubleArray(list));
        }
        return meanVal;
    }

    /**
     * Print results in statistics file.
     * @param medDepth median depth from {@link Cov2lr#readCoverageAndGetMedDepth(Set)}
     * @param factor2 factor from {@link Cov2lr#getFactor2(double, Set)}
     * @param sampleMedian sample median from {@link Cov2lr#getSampleMedian(Set, Set)}
     * @param controlSamplesMean mean value from control samples
     * @param badGenes set of {@link Sample#name} that failed
     * @return set of samples that should be printed in statistics file and print them
     */
    //perl version: 116 str
    List<Sample> setNormalization(double medDepth,
                                  Map<String, Double> factor2,
                                  Map<String, Double> sampleMedian,
                                  double controlSamplesMean,
                                  Set<String> badGenes) {
        List<Sample> sampleResult = new ArrayList<>();
        try {
            PrintWriter writer = getPrintWriter();
            boolean useControlSamples = useControlSamples();
            CloseableIterator<Sample> iterator = getFileDataIterator();
            int count = 0;
            while (iterator.hasNext()) {
                Sample sample = iterator.next();
                if (badGenes.contains(sample.getName())) { //perl version: 108 str
                    continue;
                }
                String key = sample.getName();
                if (!factor2.containsKey(key)) {
                    continue;
                }
                count++;
                double fact2 = factor2.get(key);
                double norm1 = geneStatistics.get(key).get(sample.getSample());
                double smplMed = sampleMedian.get(sample.getSample());
                sample.setNorm1b(Precision.round(norm1 * fact2 + 0.1, 2)); //perl version: 110 str
                sample.setNorm2(Precision.round(medDepth != 0 ? log.value((norm1 * fact2 + 0.1) / medDepth) / LOG_OF_TWO : 0, 2)); //perl version: 111 str
                sample.setNorm3(Precision.round(smplMed != 0 ? log.value((norm1 * fact2 + 0.1) / smplMed) / LOG_OF_TWO : 0, 2)); //perl version: 112 str
                //perl version: 122 str
                if (useControlSamples) {
                    sample.setNorm3s(Precision.round(controlSamplesMean != 0 ? sample.getNorm1b() / controlSamplesMean / log.value(2) : 0, 2));
                }
                sampleResult.add(sample);
                writer.write(sample.getResultString() + "\n");
                writer.flush();
                addLocus(sample); //perl version: 120 str
            }
            iterator.close();
            writer.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        return sampleResult;
    }

    PrintWriter getPrintWriter() throws FileNotFoundException {
        return new PrintWriter(tempFile);
    }

    //perl version: 120 str
    void addLocus(Sample sample) {
        Locus locus = new Locus(sample.getGene(), sample.getStart(), sample.getEnd(), sample.getChr());
        long len = sample.getEnd() - sample.getStart() + 1;
        locus.setLength(len);
        if (locusMap.containsKey(sample.getGene())) {
            Locus previous = locusMap.get(sample.getGene());
            locus.addLength(previous.length);
            locus.shiftStart(previous.start);
            locus.shiftEnd(previous.end);
        }
        locusMap.put(sample.getGene(), locus);
    }

    /**
     * Calculate median for each sample from coverage file
     * @param sampleNames set of sample names from coverage file
     * @param badGenes set of {@link Sample#name} that failed
     * @return median by sample
     */
    //perl version: 96 str
    Map<String, Double> getSampleMedian(Set<String> sampleNames, Set<String> badGenes) {
        Map<String, Double> sampleMedian = new LinkedHashMap<>();
        for (String s : sampleNames) {
            List<Double> list = new LinkedList<>();
            for (Map.Entry<String, Map<String, Double>> gene : geneStatistics.entrySet()) {
                if (badGenes.contains(gene.getKey())) {
                    continue;
                }
                list.add(gene.getValue().get(s));
            }
            sampleMedian.put(s, median.evaluate(toDoubleArray(list)));
        }
        return sampleMedian;
    }

    /**
     * Gene/amplicon factor
     * @param medDepth median of result {@link Cov2lr#splitQualitySamples(Set, double, Set)}
     * @param badGenes set of {@link Sample#name} filled in {@link Cov2lr#splitQualitySamples(Set, double, Set)}
     * @return map: key - {@link Sample#name}, value - medDepth / median
     */
    //perl version: 89 str
    Map<String, Double> getFactor2(final double medDepth, Set<String> badGenes) {
        Map<String, Double> factor2 = new HashMap<>();
        for (Map.Entry<String, Map<String, Double>> gene : geneStatistics.entrySet()) {
            String key = gene.getKey();
            if (badGenes.contains(key)) {
                continue;
            }
            double[] norms1 = new double[gene.getValue().size()];
            int idx = 0;
            for (Double normValue : gene.getValue().values()) {
                norms1[idx++] = normValue;
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

    /**
     * Split up genes/amplicons from coverage file
     * @param samp set of sample names from coverage file
     * @param medDepth median of result {@link Cov2lr#splitQualitySamples(Set, double, Set)}
     * @param badGenes set of {@link Sample#name} filled in {@link Cov2lr#splitQualitySamples(Set, double, Set)} that failed
     * @return List of good depth - calculated norm1 in according to 80% percentile
     */
    //perl version: 76 str
    List<Double> splitQualitySamples(Set<String> samp, double medDepth, Set<String> badGenes) {
        List<Double> gooddepth = new LinkedList<>();
        for (Map.Entry<String, Map<String, Double>> gene : geneStatistics.entrySet()) {
            List<Double> temp = new LinkedList<>();
            double kp80 = filterDataStream(samp, gene.getValue(), temp); //perl version: 80 str
            if (kp80 < medDepth * FAILEDFACTOR) {
                badGenes.add(gene.getKey());
            } else {
                gooddepth.addAll(temp);
            }
        }

        return gooddepth;
    }

    /**
     * Read samples from coverage file, put them in map, calculates median depth and fill in map {@link Cov2lr#geneStatistics}
     * @param samp set of sample names from coverage file; empty
     * @return median depth
     */
    //perl version: 63 str
    double readCoverageAndGetMedDepth(Set<String> samp) {
        List<Double> depth = new ArrayList<>();
        CloseableIterator<Sample> iterator = getFileDataIterator();
        while (iterator.hasNext()) {
            Sample sample = iterator.next();
            double norm1 = Precision.round((sample.getCov() * factor.get(sample.getSample())), 2);
            depth.add(norm1);
            samp.add(sample.getSample());
            putInMap(sample, norm1); //perl version: 68 str
        }
        iterator.close();
        return median.evaluate(toDoubleArray(depth));
    }

    CloseableIterator<Sample> getFileDataIterator() {
        return new FileDataIterator(covFile, factory);
    }

    //perl version: 68 str
    // fill in geneStatistics
    void putInMap(Sample sample, double norm1) {
        String gene = sample.getName();
        Map<String, Double> genesToNorm = geneStatistics.get(gene);
        if (genesToNorm == null) {
            genesToNorm = new HashMap<>();
        }
        genesToNorm.put(sample.getSample(), norm1);
        geneStatistics.put(gene, genesToNorm);
    }

    /**
     * Evaluates 80% percentile
     * @param sampleSet set of sample names from coverage file
     * @param map value from geneStatistics map
     * @param list list of numbers from geneStatistics by sample
     * @return 80% percentile
     */
    //perl version: 80 str
    double filterDataStream(Set<String> sampleSet, Map<String, Double> map, List<Double> list) {
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

    public Map<String, Locus> getLocusMap() {
        return locusMap;
    }

    Map<String, Map<String, Double>> getGeneStatistics() {
        return geneStatistics;
    }
}