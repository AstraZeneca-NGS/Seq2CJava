package com.astrazeneca.seq2c;

import com.astrazeneca.seq2c.input.FileDataIterator;
import com.astrazeneca.seq2c.input.Sample;
import com.astrazeneca.seq2c.input.SampleStatistics;
import com.astrazeneca.seq2c.input.StatisticsFactory;
import org.apache.commons.cli.*;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.util.Precision;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Lr2gene {

    private Map<String, Locus> genes;
    private double MINMAD = 10;
    private double MINDIFF = 0.7;
    private double PVALUE = 0.00001;
    private double AMP = 1.5;
    private double DEL = -2.00;
    private double EXONDIFF = 1.25;
    private double MAXRATE = 0.1; // maximum
    private double MAXCNT = 5;
    private double MINTTDIFF = 0.7;
    private double TTPVALUE = 0.000001;
    private double MINBPEXONS = 8;
    private double MINSEGS = 1;
    private boolean useControl;
    private static final Lr2gene DEFAULTS = new Lr2gene();
    private static final double EPSILON = 1E-8;

    private static final Pattern SLASH_D = Pattern.compile("\\d");

    private String tempFile;

    private static final java.util.Comparator<double[]> DIS_COMPARATOR = new java.util.Comparator<double[]>() {
        @Override
        public int compare(double[] a, double[] b) {
            return Double.compare(b[0], a[0]);
        }
    };

    Lr2gene() {
    }

    public Lr2gene(Map<String, Locus> locusMap, String file) {
        this.genes = locusMap;
        this.tempFile = file;
    }

    public void process() {
        int j = 0;
        FileDataIterator<SampleStatistics> iterator = new FileDataIterator(tempFile, new StatisticsFactory());
        Map<String, Map<String, Sig>> sampleSigs = new HashMap<>();
        Map<String, Integer> geneSigs = new HashMap<>();

        while (iterator.hasNext()) {
            SampleStatistics statistics = iterator.next();
            String sample = statistics.getSample();
            Map<String, Sig> results = new HashMap<>();
            sampleSigs.put(sample, results);
            Map<String, List<Sample>> gq = statistics.getGenes();
            for (Map.Entry<String, List<Sample>> entry2 : gq.entrySet()) {
                String gene = entry2.getKey();
                List<Sample> sq2amparr = entry2.getValue();

                Collections.sort(sq2amparr, new Comparator<Sample>() {
                    @Override
                    public int compare(Sample ints, Sample otherInts) {
                        return Double.compare(ints.getStart(), otherInts.getStart());
                    }
                });

                ArrayList<Double> lr = new ArrayList<>(sq2amparr.size());
                for (Sample sqarr : sq2amparr) {
                    if (useControl) {
                        lr.add(sqarr.getNorm3s());
                    } else {
                        lr.add(sqarr.getNorm3());
                    }
                }
                double lr_med;
                if (lr.size() > 1) {
                    lr_med = StatUtils.percentile(convertDoubles(lr), 50);// median
                } else {
                    lr_med = lr.get(0);
                }
                Sig sig = checkBP(sq2amparr);
                if (sig == null || Double.compare(sig.getSig(), -1.0) == 0) {
                    if (lr_med >= AMP) {
                        sig = new Sig(0.0, 0.0, lr.size(), "Whole", lr_med, "Amp", lr.size(), 0.0, "ALL", lr_med);
                    } else if (lr_med <= DEL) {
                        sig = new Sig(0.0, 0.0, lr.size(), "Whole", lr_med, "Del", lr.size(), 0.0, "ALL", lr_med);
                    } else {
                        sig = new Sig(-1, 0, 0, "", 0, "", lr.size(), 0, "", 0);
                    }
                }
                Matcher matcher = SLASH_D.matcher(sig.getSigseg());
                if (!sig.getSigseg().isEmpty() && matcher.find()) {
                    String[] exons = sig.getSigseg().split(",");
                    long estart = sq2amparr.get(Integer.valueOf(exons[0]) - 1).getStart();
                    long eend = sq2amparr.get(Integer.valueOf(exons[exons.length - 1]) - 1).getEnd();
                    sig.setSigseg(sig.getSigseg() + "(" + estart + "-" + eend + ")");
                }

                sig.setLrMed(lr_med);

                if (sig.getBp().equals("BP")) {
                    String key = (gene + " " + sig.getSigseg()).intern();
                    if (geneSigs.containsKey(key)) {
                        geneSigs.put(key, geneSigs.get(key) + 1);
                    } else {
                        geneSigs.put(key, 1);
                    }
                }
                results.put(gene, sig);
            }
        }
        iterator.close();
        printResults(sampleSigs, geneSigs);
    }

    private void printResults(Map<String, Map<String, Sig>> sampleSigs, Map<String, Integer> geneSigs) {
        int j = 0;
        for (Map.Entry<String, Map<String, Sig>> entry : sampleSigs.entrySet()) {
            String sample = entry.getKey();
            int samplesNumber = sampleSigs.size();
            for (Map.Entry<String, Sig> sigs : entry.getValue().entrySet()) {
                String gene = sigs.getKey();
                Sig sig = sigs.getValue();
                Locus locus = genes.get(gene);
                StringBuilder locStr = new StringBuilder();

                locStr.append(sample).append("\t").append(gene).append("\t").append(locus.getName()).append("\t").append(Precision.round(sig.getLrMed(), 3)).append("\t");
                if (sig.getSig() != -1) {
                    double bpPercent = 0.0;
                    int bpNumber = 0;
                    if (sig.getBp().equals("BP")) {
                        String key = (gene + " " + sig.getSigseg()).intern();
                        bpNumber = geneSigs.get(key);
                        bpPercent = Precision.round(bpNumber/(double)samplesNumber, 3);
                    }
                    if (Double.compare(bpPercent, 0.0) != 0 && Double.compare(bpPercent, MAXRATE) > 0
                            && bpNumber > MAXCNT) {
                        locStr.append("\t\t\t\t\t").append(sig.getTotal());
                    } else {
                        locStr.append(sig.getName().append("\t").append(bpNumber).append("\t").append(bpPercent));
                    }
                } else {
                    locStr.append("\t\t\t\t").append(sig.getTotal());
                }
                if (j == 0)
                    System.out.println("Sample\tGene\tChr\tStart\tEnd\tLength\tLog2ratio\tSig\tBP_Whole\tAmp_Del\tAb_Seg\tTotal_Seg\tAb_log2ratio\tLog2r_Diff\tAb_Seg_Loc\tAb_Samples\tAb_Samples_Pcnt");
                j++;


                System.out.println(locStr);
            }
        }
    }

    Sig checkBP(List<Sample> segs) {
        if (segs.size() < 4) return null;

        double[][] arr = new double[segs.size()][3];
        double[] lr = new double[segs.size()];
        int i = 0;
        for (Sample sqarr : segs) {
            arr[i][0] = sqarr.getStart();
            arr[i][1] = useControl ? sqarr.getNorm3s() : sqarr.getNorm3();
            arr[i][2] = i+1;
            i++;
        }
        i = 0;
        for (double[] l : arr) {
            lr[i] = l[1];
            i++;
        }

        double[] bps = getBPS(lr);
        double minbp = 1;
        double maxmd = 1;
        Sig sig = new Sig();
        Sig sigbp = new Sig();
        Sig sigmd = new Sig();
        for (double bp : bps) {

            List<double[]> bm = new ArrayList<>();
            List<double[]> up = new ArrayList<>();
            List<Double> lrup = new ArrayList<>();
            List<Double> lrbm = new ArrayList<>();
            List<Integer> upsegArr = new ArrayList<>();
            List<Integer> bmsegArr = new ArrayList<>();
            i = 0;
            for (double[] intarr : arr) {
                if (intarr[1] > bp) {
                    up.add(intarr);
                    lrup.add(intarr[1]);
                    upsegArr.add(i + 1);
                } else {
                    bm.add(intarr);
                    lrbm.add(intarr[1]);
                    bmsegArr.add(i + 1);
                }
                i++;
            }
            double lrupm = StatUtils.percentile(convertDoubles(lrup), 50);
            double lrbmm = StatUtils.percentile(convertDoubles(lrbm), 50);
            String cn = Double.compare(lrbmm, -0.35) < 0 ? "Del" : (lrupm > 0.35 && Math.abs(lrbmm) < Math.abs(lrupm) ? "Amp" : "NA");
            if ("NA".equals(cn)) {
                continue;
            }
            double[] bmiscArr = isConsecutive(bm);
            double[] upiscArr = isConsecutive(up);
            if (bmiscArr[0] != 0) {
                if (bmiscArr[1] != -1) {
                    int ti = 0;
                    for (int k = 0; k < up.size(); k++) {
                        if (up.get(k)[2] == bmiscArr[1])
                            ti = k;
                    }
                    int index = (int)bmiscArr[1];
                    index = index < 0 ? bm.size() - 1 + index : index;
                    if (index > bm.size()) {
                        bm.add(up.remove(ti));
                    } else {
                        bm.add(index, up.remove(ti));
                    }
                }
                Sig tmpSig = getCalls(bm, up);
                double[] issig = isSig(bm, up);
                if (issig[0] >= 0 && issig[0] < minbp) {
                    sigbp.addSig(issig[0]);
                    sigbp.addSdiff(issig[1]);
                    sigbp.update(tmpSig);
                    minbp = issig[0];
                } else if (issig[0] > maxmd) {
                    sigmd.addSig(issig[0]);
                    sigmd.addSdiff(issig[1]);
                    sigmd.update(tmpSig);
                    maxmd = issig[0];
                }

            } else if (upiscArr[0] != 0) {
                if (upiscArr[1] != -1) {
                    int ti = 0;
                    for (int k = 0; k < bm.size(); k++) {
                        if (bm.get(k)[2] == upiscArr[1])
                            ti = k;
                    }
                    int index = (int) upiscArr[1];
                    index = index < 0 ? up.size() - 1 + index : index;
                    if (index >= up.size()) {
                        up.add(bm.remove(ti));
                    } else {
                        up.add(index, bm.remove(ti));
                    }
                }
                Sig tmpSig = getCalls(up, bm);
                double[] issig = isSig(up, bm);
                if (issig[0] >= 0 && issig[0] < minbp) {
                    sigbp.addSig(issig[0]);
                    sigbp.addSdiff(issig[1]);
                    sigbp.update(tmpSig);
                    minbp = issig[0];
                } else if (issig[0] > maxmd) {
                    sigmd.addSig(issig[0]);
                    sigmd.addSdiff(issig[1]);
                    sigmd.update(tmpSig);
                    maxmd = issig[0];
                }
            }
        }
        if (sigbp.getSig() != -1) {
            return sigbp;
        }
        if (sigmd.getSig() != -1) {
            return sigmd;
        }
        if (sig.getSig() == -1) {
            sig = findBP(lr);
        }

        if (sig.getSig() != -1) {
            sig.setBp("BP");
        }
        sig.setTotal(arr.length);
        return sig;
    }

    private double[] isConsecutive(List<double[]> ref) {

        double skip = 0;
        double si = -1;
        double sii = -1;
        double[] result = new double[3];
        for (int i = 1; i < ref.size(); i++) {
            skip += ref.get(i)[2] - ref.get(i - 1)[2] - 1;
            if (ref.get(i)[2] - ref.get(i - 1)[2] == 2) {
                si = ref.get(i)[2] - 1;
                sii = i;
            }
        }
        if ((skip == 0) || (skip == 1 && ref.size() >= 10)) {
            result[0] = 1;
        } else {
            result[0] = 0;
        }
        result[1] = si;
        result[2] = sii;
        return result;
    }

    private double[] getBPS(double[] lr) {
        double[] lrSorted = new double[lr.length];
        System.arraycopy(lr, 0, lrSorted, 0, lr.length);
        Arrays.sort(lrSorted);

        double[][] dis = new double[lrSorted.length-1][3];
        for (int i = 1; i < lrSorted.length; i++) {
            dis[i-1][0] = lrSorted[i] - lrSorted[i-1];
            dis[i-1][1] = lrSorted[i];
            dis[i-1][2] = lrSorted[i-1];
        }

        java.util.Arrays.sort(dis, DIS_COMPARATOR);
        ArrayList<Double> bpsArr = new ArrayList<>(dis.length);
        for (double[] bp : dis) {
            if (Precision.compareTo(bp[0], 0.1, EPSILON) < 0)
                break;
            bpsArr.add((bp[1] + bp[2]) / 2);

        }
        return convertDoubles(bpsArr);
    }

    private Sig findBP(double[] lr) {
        if (lr.length < 15) {
            Sig sig = new Sig(-1.0, 0, 0, "", 0, "", 0, 0, "", 0);
            return sig;
        }
        double minp = 1;
        int bpi = 0;
        double siglr = 0;
        String cn = "NA";
        double mindiff = 0;
        String sigseg = "";

        for (int i = (int) MINBPEXONS; i < lr.length - (int) MINBPEXONS; i++) {
            double[] lr_x = new double[i];
            double[] lr_y = new double[lr.length - i];

            for (int k = 0; k < i; k++) {
                lr_x[k] = lr[k];
            }
            Arrays.sort(lr_x);
            for (int k = i; k < lr.length; k++) {
                lr_y[k-i] = lr[k];
            }
            Arrays.sort(lr_y);

            double bpleft = StatUtils.mean(lr_x);
            double bpright = StatUtils.mean(lr_y);
            if ((bpleft > bpright && lr_x[1] < lr_y[lr_y.length - 1]) || (bpleft < bpright && lr_y[1] < lr_x[lr_x.length - 1])) {
                continue;
            }
            double p = TTestStatistics.getT(lr_x, lr_y);
            double[] sigseg1 = new double[i];
            double[] sigseg2 = new double[lr.length - i];
            int j = 0;
            for (int k = 1; k <= i; k++) {
                sigseg1[j] = k;
                j++;
            }
            j = 0;
            for (int k = i + 1; k <= lr.length; k++) {
                sigseg2[j] = k;
                j++;
            }
            double diff = StatUtils.mean(lr_x) - StatUtils.mean(lr_y);
            if ((p < minp || ((p > 0 && minp / p < 10 && Math.abs(diff) > mindiff) || (p == 0 && Math.abs(diff) > mindiff))) && ((p < TTPVALUE && Math.abs(diff) > MINTTDIFF) || (p < 0.001 && Math.abs(diff) >= MINTTDIFF && (Math.abs(bpleft) > 0.80 || Math.abs(bpright) > 0.80)))) {
                minp = p;
                bpi = Math.abs(bpleft) > Math.abs(bpright) ? i : (lr.length - i + 1);
                siglr = Math.abs(bpleft) > Math.abs(bpright) ? bpleft : bpright;
                sigseg = Math.abs(bpleft) > Math.abs(bpright) ? joinDouble(sigseg1, ",") : joinDouble(sigseg2, ",");
                cn = Math.abs(bpleft) > Math.abs(bpright) ? (bpleft < -0.5 ? "Del" : (bpleft > 0.5 ? "Amp" : "NA")) : (bpright < -0.5 ? "Del" : (bpright > 0.5 ? "Amp" : "NA"));
                mindiff = Math.abs(diff);
            }
        }
        Sig sig;
        if (minp < 1) {
            sig = new Sig(minp, minp, bpi, "", siglr, cn, 0, mindiff, sigseg, mindiff);
        } else {
            sig = new Sig(-1.0, 0, 0, "", 0, "", 0, 0, "", 0);
        }
        return sig;

    }

    private Sig getCalls(List<double[]> bm, List<double[]> up) {
        double[] tlr1 = new double[bm.size()];
        double[] ti1 = new double[bm.size()];
        double[] tlr2 = new double[up.size()];
        double[] ti2 = new double[up.size()];
        int i = 0;
        for (double[] b : bm) {
            tlr1[i] = b[1];
            ti1[i] = b[2];
            i++;
        }
        i = 0;
        for (double[] u : up) {
            tlr2[i] = u[1];
            ti2[i] = u[2];
            i++;
        }
        i = 0;
        String cn;
        String ti;
        int segs;
        double mean;
        double mean1 = StatUtils.mean(tlr1);
        double mean2 = StatUtils.mean(tlr2);
        if (Math.abs(mean1) > Math.abs(mean2)) {
            cn = mean1 < -0.35 ? "Del" : (Double.compare(mean1, 0.35) > 0 ? "Amp" : "NA");
            segs = tlr1.length;
            mean = mean1;
            ti = joinDouble(ti1, ",");
        } else {
            cn = mean2 < -0.35 ? "Del" : (Double.compare(mean2, 0.35) > 0 ? "Amp" : "NA");
            segs = tlr2.length;
            mean = mean2;
            ti = joinDouble(ti2, ",");
        }
        Sig sig = new Sig(0.0, 0.0, segs, "BP", mean, cn, bm.size() + up.size(), 0.0, ti, 0.0);
        return sig;

    }

    private double[] isSig(List<double[]> bm, List<double[]> up) {

        double[] bm_x = new double[bm.size()];
        double[] up_y = new double[up.size()];
        double[] result = new double[2];
        int i = 0;
        for (double[] b : bm) {
            bm_x[i] = b[1];
            i++;
        }
        i = 0;
        for (double[] u : up) {
            up_y[i] = u[1];
            i++;
        }

        if (bm_x.length >= 3 && up_y.length >= 3) {
            double p = TTestStatistics.getT(bm_x, up_y);
            double diff = StatUtils.mean(bm_x) - StatUtils.mean(up_y);
            if ((p < PVALUE && Math.abs(diff) >= MINDIFF) || (p < 0.001 && Math.abs(diff) >= MINDIFF && (Math.abs(StatUtils.mean(bm_x)) > 0.8 || Math.abs(StatUtils.mean(up_y)) > 0.8))) {
                result[0] = p;
                result[1] = Precision.round(Math.abs(diff), 3);
                return result;
            }
        } else if (bm_x.length >= MINSEGS && up_y.length >= 3) {
            double med = StatUtils.percentile(up_y, 50);
            double[] d = new double[up_y.length];
            i = 0;
            for (double m : up_y) {
                d[i] = Math.abs(m - med);
                i++;
            }
            double mad = StatUtils.percentile(d, 50);
            if (mad == 0)
                mad += 0.1;
            double[] t = new double[bm_x.length];
            i = 0;
            for (double m : bm_x) {
                t[i] = (m - med) / mad;
                i++;
            }
            double mean = StatUtils.mean(t);
            double sum = StatUtils.sum(t);
            double diff = Math.abs(StatUtils.mean(bm_x) - StatUtils.mean(up_y));
            if (Precision.compareTo(Math.abs(sum), MINMAD, EPSILON) > 0 && Precision.compareTo(diff, EXONDIFF, EPSILON) > 0) {
                result[0] = Math.abs(mean);
                result[1] = diff;
                return result;
            }

        } else if (up_y.length >= MINSEGS && bm_x.length >= 3) {

            double med = StatUtils.percentile(bm_x, 50);

            double[] d = new double[bm_x.length];
            i = 0;
            for (double m : bm_x) {
                d[i] = Math.abs(m - med);
                i++;
            }
            double mad = StatUtils.percentile(d, 50);
            if (mad == 0)
                mad += 0.1;
            double[] t = new double[up_y.length];
            i = 0;
            for (double m : up_y) {
                t[i] = (m - med) / mad;
                i++;
            }
            double mean = StatUtils.mean(t);
            double sum = StatUtils.sum(t);
            double diff = Math.abs(StatUtils.mean(bm_x) - StatUtils.mean(up_y));
            if (Precision.compareTo(Math.abs(sum), MINMAD, EPSILON) > 0 && Precision.compareTo(diff, EXONDIFF, EPSILON) > 0) {
                result[0] = Math.abs(mean);
                result[1] = diff;
                return result;
            }

        } else {
            result[0] = -1;
            result[1] = 0;
        }
        result[0] = -1;
        result[1] = 0;
        return result;
    }

    private String joinDouble(double[] doubles, String delim) {
        int i = 1;
        String result = new String();
        for (double d : doubles) {
            String r = String.valueOf((int)d);
            if (i < doubles.length) {
                r += delim;
            }
            result += r;
            i++;
        }

        return result;
    }

    public boolean isUseControl() {
        return useControl;
    }

    public void setUseControl(boolean useControl) {
        this.useControl = useControl;
    }

    private double[] convertDoubles(List<Double> doubles) {
        double[] ret = new double[doubles.size()];
        Iterator<Double> iterator = doubles.iterator();
        int i = 0;
        while (iterator.hasNext()) {
            ret[i] = iterator.next().doubleValue();
            i++;
        }
        return ret;
    }

    public void init(CommandLine cmd) throws ParseException {
        MINMAD = getDoubleValue(cmd, "M", DEFAULTS.MINMAD);
        MINDIFF = getDoubleValue(cmd, "d", DEFAULTS.MINDIFF);
        PVALUE = getDoubleValue(cmd, "p", DEFAULTS.PVALUE);
        AMP = getDoubleValue(cmd, "A", DEFAULTS.AMP);
        DEL = getDoubleValue(cmd, "D", DEFAULTS.DEL);
        EXONDIFF = getDoubleValue(cmd, "E", DEFAULTS.EXONDIFF);
        MAXRATE = getDoubleValue(cmd, "R", DEFAULTS.MAXRATE);
        MAXCNT = (int) getDoubleValue(cmd, "N", DEFAULTS.MAXCNT);
        MINTTDIFF = getDoubleValue(cmd, "t", DEFAULTS.MINTTDIFF);
        TTPVALUE = getDoubleValue(cmd, "P", DEFAULTS.TTPVALUE);
        MINBPEXONS = getDoubleValue(cmd, "e", DEFAULTS.MINBPEXONS);
        if (cmd.hasOption("c")) {
            useControl = true;
        }
    }

    private static void help(Options options) {
        HelpFormatter formater = new HelpFormatter();
        formater.setOptionComparator(null);
        formater.printHelp(142,
                Lr2gene.class.getSimpleName() + " [-aPH] [-c control] [-F float] [-s min_amplicon_#] [-A float] [-D float] mapping_reads coverage.txt\n\n"
                        + "The " + Lr2gene.class.getSimpleName() + " program will convert a coverage file to copy number profile.",

                "\nArguments are:\n"
                        + "mapping_reads: Required. A file containing # of mapped or sequenced reads for genes.  At least two columns.\n"
                        + "               First is the sample name, 2nd is the number of mapped or sequenced reads.\n"
                        + "coverage.txt:  The coverage output file from checkCov.pl script.  Can also take from standard in or more than\n"
                        + "               one file.\n\n",
                options,
                "\nAUTHOR\n"
                        + "       Written by Zhongwu Lai, AstraZeneca, Boston, USA\n\n"
                        + "REPORTING BUGS\n"
                        + "       Report bugs to zhongwu@yahoo.com\n\n"
                        + "COPYRIGHT\n"
                        + "       This is free software: you are free to change and redistribute it.  There is NO WARRANTY, to the extent permitted by law.");

        System.exit(0);

    }

    public static CommandLine buildCommandLine(String[] args) throws ParseException {
        Options options = getOptions();
        CommandLineParser parser = new BasicParser();
        CommandLine cmd = parser.parse(options, args);
        return cmd;
    }

    public static Options getOptions() {
        Options options = new Options();

        options.addOption(OptionBuilder.withArgName("float")
                .hasArg()
                .withType(Number.class)
                .withDescription("When considering partial deletions less than 3 exons/amplicons, the minimum MAD value, in addition to -d,\n"
                        + "before considering it to be amplified or deleted.  Default: " + DEFAULTS.MINMAD)
                .isRequired(false)
                .create("M"));
        options.addOption(OptionBuilder
                .withDescription("Indidate that control sample is used for normalization")
                .isRequired(false)
                .create("c"));
        options.addOption(OptionBuilder.withArgName("float")
                .hasArg()
                .withType(Number.class)
                .withDescription("When considering >=3 exons deletion/amplification within a gene, the minimum differences between the log2 of two segments.\n"
                        + "Default: " + DEFAULTS.MINDIFF)
                .isRequired(false)
                .create("d"));
        options.addOption(OptionBuilder.withArgName("float (0-1)")
                .hasArg()
                .withType(Number.class)
                .withDescription("The p-value for t-test when consecutive exons/amplicons are >= 3.  Default: " + String.format("%f", DEFAULTS.PVALUE))
                .isRequired(false)
                .create("p"));
        options.addOption(OptionBuilder.withArgName("float")
                .hasArg()
                .withType(Number.class)
                .withDescription("Minimum log2 ratio for a whole gene to be considered amplified.  Default: " + DEFAULTS.AMP)
                .isRequired(false)
                .create("A"));
        options.addOption(OptionBuilder.withArgName("float")
                .hasArg()
                .withType(Number.class)
                .withDescription("Minimum log2 ratio for a whole gene to be considered deleted.  Default: " + DEFAULTS.DEL)
                .isRequired(false)
                .create("D"));
        options.addOption(OptionBuilder.withArgName("float")
                .hasArg()
                .withType(Number.class)
                .withDescription("Minimum mean log2 ratio difference for <3 exon deletion/amplification to be called.  Default: " + DEFAULTS.EXONDIFF)
                .isRequired(false)
                .create("E"));
        options.addOption(OptionBuilder.withArgName("float (0-1)")
                .hasArg()
                .withType(Number.class)
                .withDescription("If a breakpoint has been detected more than \"float\" fraction of genes, it's considered false positive and removed.\n"
                        + "Default: " + DEFAULTS.MAXRATE + ", or 10%.  Use in combination with -N")
                .isRequired(false)
                .create("R"));
        options.addOption(OptionBuilder.withArgName("int")
                .hasArg()
                .withType(Number.class)
                .withDescription("If a breakpoint has been detected more than \"int\" genes, it's considered false positives and removed.\n"
                        + "Default: " + (int) DEFAULTS.MAXCNT + ".  Use in combination with -R.")
                .isRequired(false)
                .create("N"));
        options.addOption(OptionBuilder.withArgName("float")
                .hasArg()
                .withType(Number.class)
                .withDescription("When considering breakpoint in the middle of a gene, the minimum differences between the log2 of two segments."
                        + "Default: " + DEFAULTS.MINTTDIFF)
                .isRequired(false)
                .create("t"));
        options.addOption(OptionBuilder.withArgName("float (0-1)")
                .hasArg()
                .withType(Number.class)
                .withDescription("The p-value for t-test when the breakpoint is in the middle with min exons/amplicons >= [-e].  Default: " + String.format("%f", DEFAULTS.TTPVALUE))
                .isRequired(false)
                .create("P"));
        options.addOption(OptionBuilder.withArgName("float")
                .hasArg()
                .withType(Number.class)
                .withDescription("When considering breakpoint in the middle of a gene, the minimum number of exons.  Default: " + DEFAULTS.MINBPEXONS)
                .isRequired(false)
                .create("e"));

        return options;
    }

    private static double getDoubleValue(CommandLine cmd, String opt, double defaultValue) throws ParseException {
        Object value = cmd.getParsedOptionValue(opt);
        return value == null ? defaultValue : ((Number) value).doubleValue();
    }

}