package seq2c;

import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 *
 * @author Petr_Rastegaev
 */
import org.apache.commons.math3.analysis.function.Log;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.commons.math3.util.Precision;

/**
 *  Normalize the coverage from targeted sequencing to CNV log2 ratio.  The algorithm assumes the medium
 *     is diploid, thus not suitable for homogeneous samples (e.g. parent-child).
 */


public class Cov2lr {
    /**
     * Help utils for statistics
    */
    private Log log = new Log();
    private Median median = new Median();
    private Mean mean = new Mean();
    /**
      *  Indicate this is amplicon or exon based calling.  By default, it'll aggregate at gene level.
     */
    private boolean amplicon;

    /**
     * Use the control sample(s)
     */
    private boolean useControlSamples;

    /**
      *  Array of control samples
     */
    private String[] controlSamples;


    /**
     *   Statistics from first input file
     *   Key = name of sample
     *   Value = number of aligned reads
     */
    private Map<String, Long> mappingReads;

    /**
      *  Normalization factor for each sample
      *  calculated as mean number of reads / number of reads per sample
      *  Key = name of sample
      *  Value = factor
     */
    private Map<String, Double> factor;
    /**
     *  Map of genes from second input file
     *  Key = name of gene
     *  Value = gene object
     */
    private Map<String, Gene> genes;
    /**
     *  Map of samples
     *  Key = key for string, consists of sample name, gene name, start, end and length
     *  Value = map of sample names and sample objects
     */
    private Map<String, Map<String, Sample>> samples;

    /**
     *  The failed factor for individual amplicons.
     *  If (the 80th percentile of an amplicon depth)/(the global median depth)
     *  is less than the argument, the amplicon is considered failed and won't be used in calculation.  Default: 0.2.
     */
    private final double FAILEDFACTOR = 0.2;

    /**
     * Constructor reads the input files and constructs genes, samples, factor maps
     * @param amplicon  =   determines the aggregation level (gene or record)
     * @param fileStat  =   A file containing # of mapped or sequenced reads for samples.  At least two columns.
     *                      First is the sample name, 2nd is the number of mapped or sequenced reads.
     * @param fileCov   =   The coverage output file from checkCov.pl script.  Can also take from standard in or more than
     *                      one file.
     * @param useControlSamples = use the control sample(s)
     * @param controlSamples     = multiple controls are allowed, which are separated by ":"
     *
     */

    public Cov2lr(boolean amplicon, Map<String, Long> stat, Collection<Gene> genesArray, boolean useControlSamples, String controlSamples) {
        this.amplicon = amplicon;
        this.useControlSamples = false;
        this.mappingReads = stat;
        this.genes = new LinkedHashMap<>();
        this.samples = new LinkedHashMap<>();
        readCoverage(genesArray, this.genes, this.samples);
        this.factor = new LinkedHashMap<>();
        setFactor();
        this.useControlSamples = useControlSamples;
        if (useControlSamples) {
            this.controlSamples = controlSamples.split(":");
        }
    }


    /**
     * Main method, makes the statistics calculation according to the algorithm
     * print result to the Standart output
     */

    public ArrayList<Sample> doWork() {
        double[] depth = new double[samples.size()];

        Set<String> samp = new HashSet<>();

        double medDepth = getMedDepth(depth, samp);

        List<String> bad = new LinkedList<>();
        List<Double> gooddepth = new LinkedList<>();

        splitQualitySamples(samp, medDepth, bad, gooddepth);

        medDepth = median.evaluate(toDoubleArray(gooddepth));

        Map<String, Double> factor2 = new LinkedHashMap<>();

        setFactor2(medDepth, bad, factor2);

        Map<String, Double> sampleMedian = getSampleMedian(samp, bad);

        setNorm(medDepth, bad, factor2, sampleMedian);

        return printResult(bad);
    }


//    private Map<String, Long> readStat(String fileName) {
//        Map<String, Long> map = new LinkedHashMap<>();
//        try (BufferedReader reader= new BufferedReader(new FileReader(fileName))) {
//            String line;
//            while ((line = reader.readLine()) != null) {
//                String[] samples = line.split("\\s+");
//                if (samples.length == 2) {
//                    map.put(samples[0], Long.parseLong(samples[1]));
//                }
//            }
//        } catch (IOException e) {
//            System.err.println("Cannot open file " + fileName);
//        } catch (NumberFormatException e) {
//            System.err.println("Cannot parse long number");
//        }
//        return map;
//    }

    private void readCoverage(Collection<Gene> genesArray, Map<String, Gene> genes, Map<String, Map<String, Sample>> samples) {

        for (Gene gn : genesArray) {
            if (gn.getTag().contains("Whole")) {
                continue;
            }

            String sample = gn.getSample();
            String gene = gn.getName();
            String chr = gn.getChr();
            String tag = gn.getTag();
            long start = gn.getStart();
            long end = gn.getEnd();
            long len = gn.getLen();
            double depth = gn.getDepth();
            addGene(genes, sample, gene, chr, start, end, tag, len, depth);
            String key = amplicon ? join(gene, chr, start, end, len) : gene;
            addSample(samples, key, sample, chr, start, end, len, depth, gene);

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
                //System.out.println("addcov: "+ gene +" "+ depth);
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
            for (Map.Entry<String, Double> entry: factor.entrySet()) {
                writer.write(entry.getKey() + "\t" + entry.getValue() + "\n");
                writer.flush();
            }
            writer.write("*********************" + "\n");
            writer.write("Loc hash:" + "\n");
            writer.flush();
            long start = 0;
            long end = 0;
            long len = 0;
            for (Map.Entry<String, Gene> entry : genes.entrySet()) {
                start += entry.getValue().getStart();
                end += entry.getValue().getEnd();
                len += entry.getValue().getLen();

                //writer.write(entry.getKey() + "\t" + entry.getValue() + "\n");
                //writer.flush();
            }

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

    private ArrayList<Sample> printResult(List<String> bad) {
        ArrayList<Sample> sampleArr = new ArrayList<>();
        for (Map.Entry<String, Map<String, Sample>> entry : samples.entrySet()) {
            if (!bad.contains(entry.getKey())) {
                for (Map.Entry<String, Sample> entrySample : entry.getValue().entrySet()) {
                    Sample sample = entrySample.getValue();
//                    String title = amplicon ? entry.getKey() : sample.getTitle(genes.get(sample.getGene()));
//                    String result = sample.getResultString(title);

                    if (useControlSamples) {
                        sample = addControlSamples(sample);
                    }
                    sampleArr.add(sample);
                    //System.out.println(result);
                }
            }
        }
    return sampleArr;
    }

    private Sample addControlSamples(Sample sample) {
        List<Double> list = new LinkedList<>();
        for (String s: controlSamples) {
            for (Map.Entry<String, Map<String, Sample>> entry : samples.entrySet()) {
                list.add(entry.getValue().get(s).getNorm1b());
            }
        }
        if (!list.isEmpty()) {
            double meanVal = mean.evaluate(toDoubleArray(list));

            sample.setNorm3s( meanVal == 0 ? sample.getNorm1b() / meanVal / log.value(2) : 0);
        }
        return sample;
    }

    private void setNorm(double medDepth, List<String> bad, Map<String, Double> factor2, Map<String, Double> sampleMedian) {
        for (Map.Entry<String, Map<String, Sample>> entry : samples.entrySet()) {
            if (!bad.contains(entry.getKey())) {
                for (Map.Entry<String, Sample> entrySample : entry.getValue().entrySet()) {
                    Sample sample = entrySample.getValue();
                    double norm1 = sample.getNorm1();
                    double fact2 = factor2.get(entry.getKey());
                    double smplMed = sampleMedian.get(sample.getSample());
                    sample.setNorm1b(Precision.round(norm1 * fact2 + 0.1, 2));
                    sample.setNorm2(Precision.round(medDepth != 0 ? log.value((norm1 * fact2 + 0.1) / medDepth) / log.value(2) : 0, 2));
                    sample.setNorm3(Precision.round(smplMed != 0 ? log.value((norm1 * fact2 + 0.1) / smplMed) / log.value(2) : 0, 2));

                }
            }
        }
    }

    private Map<String, Double> getSampleMedian(Set<String> samp, List<String> bad) {
        Map<String, Double> sampleMedian = new LinkedHashMap<>();

        for (String s : samp) {
            List<Double> list = new LinkedList<>();
            for (Map.Entry<String, Map<String, Sample>> entry : samples.entrySet()) {
                if (!bad.contains(entry.getKey())) {
                    list.add(entry.getValue().get(s).getNorm1());
                }
            }
            sampleMedian.put(s, median.evaluate(toDoubleArray(list)));
        }
        return sampleMedian;
    }

    private void setFactor2(double medDepth, List<String> bad, Map<String, Double> factor2) {
        for (Map.Entry<String, Map<String, Sample>> entry : samples.entrySet()) {
            if (!bad.contains(entry.getKey())) {
                List<Double> list = new LinkedList<>();
                for (Map.Entry<String, Sample> entrySample : entry.getValue().entrySet()) {
                    Sample sample = entrySample.getValue();
                    list.add(sample.getNorm1());
                }
                double median = new Median().evaluate(toDoubleArray(list));
                if (median != 0) {
                    factor2.put(entry.getKey(), medDepth/median);
                } else {
                    factor2.put(entry.getKey(), 0.0);
                }
            }
        }
    }

    private void splitQualitySamples(Set<String> samp, double medDepth, List<String> bad, List<Double> gooddepth) {
        for (Map.Entry<String, Map<String, Sample>> entry : samples.entrySet()) {
            List<Double> temp = new LinkedList<>();
            double kp80 = filterData(samp, entry.getValue(), temp);
            if (kp80 < medDepth * FAILEDFACTOR) {
                bad.add(entry.getKey());
            } else {
                gooddepth.addAll(temp);
            }
        }
    }

    private double getMedDepth(double[] depth, Set<String> samp) {
        int i = 0;
        for (Map.Entry<String, Map<String, Sample>> entry: samples.entrySet()) {
            Map<String, Sample> sampleMap =  entry.getValue();

            for(Map.Entry<String, Sample> sampleEntry : sampleMap.entrySet()) {
                Sample sample = sampleEntry.getValue();
                double norm1 = Precision.round((sample.getCov() * factor.get(sample.getSample())), 2);
                //System.out.println("norm1" + sample.getName() + " " + sample.getCov() + " " + factor.get(sample.getSample()));
                sample.setNorm1(norm1);
                depth[i++] = norm1;
                samp.add(sample.getSample());
            }
        }
        return median.evaluate(depth);
    }


    private double filterData(Set<String> sampleSet, Map<String, Sample> map, List<Double> list) {
        for (String sample : sampleSet) {
            if (map.containsKey(sample)) {
                list.add(map.get(sample).getNorm1());
            }
        }
        return new Percentile().evaluate(toDoubleArray(list), 80);
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

