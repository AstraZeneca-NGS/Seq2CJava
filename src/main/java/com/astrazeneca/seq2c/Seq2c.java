package com.astrazeneca.seq2c;

import org.apache.commons.cli.*;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

/**
 * The Seq2C starting point class: contains {@link Seq2c#main(String[])} method that initializes all modules and runs the
 * program. Also, this class contains methods for printing output information into file.
 */
public class Seq2c {

    private static final String TEMP_FILE = "statistics.txt";

    public static void main(String[] args) throws Exception {

        Options options = buildCmdOptions();
        CommandLine cmd = null;
        try {
            cmd = new BasicParser().parse(options, args);
        } catch (ParseException e) {
            System.out.println(e.getMessage());
            help(options);
        }

        String[] cmdArgs = cmd.getArgs();
        if (cmd.getArgList().isEmpty() || cmd.hasOption('h') || cmdArgs.length < 3) {
            help(options);
        }
        String sam2bamFile = cmdArgs[0];
        final String bedFile = cmdArgs[1];
        final String covFile = cmdArgs[2];
        String control = "";

        if (cmdArgs.length > 3) {
            control = cmdArgs[3];
        }

        Dispatcher.init(getThreadsCount(cmd));
        try {

            int part = getRunPart(cmd);
            // first part is not processed for -r == 2
            if (part != 2) {

                final Map<String, String> sam2bam = Bam2Reads.parseFile(sam2bamFile);


                boolean append = false;
                for (final Map.Entry<String, String> entry : sam2bam.entrySet()) {
                    Seq2cov sec2cov = new Seq2cov(bedFile, entry.getValue(), entry.getKey());
                    Collection<Gene> cov = sec2cov.process();
                    //print coverage to the file
                    printCoverage(cov, covFile, append);
                    append = true;
                }
                //if -r = 1, we finish here
                if (part == 1) {
                    return;
                }
            }

            Map<String, Long> stat = Bam2Reads.printStatsToFile(sam2bamFile);
            printStat(stat, "read_stats.txt");

            Cov2lr cov2lr = new Cov2lr(true, stat, covFile, control, TEMP_FILE);
            cov2lr.doWork();

            Lr2gene lr2gene = new Lr2gene(cov2lr.getLocusMap(), TEMP_FILE);
            lr2gene.init(cmd);
            lr2gene.setUseControl(cov2lr.useControlSamples());
            lr2gene.process();

        } finally {
            Dispatcher.shutdown();
        }
    }

    static void printStat(Map<String, Long> sam2bam, String readStatFile) {
        try {
            FileWriter writer = new FileWriter(readStatFile);
            for (Map.Entry<String, Long> entry : sam2bam.entrySet()) {
                writer.write(entry.getKey() + "\t" + entry.getValue() + "\n");
                writer.flush();
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    //perl version: seq2cov 130 str
    static void printCoverage(Collection<Gene> sqrlist, String covFile, boolean append) {
        try (FileWriter writer = new FileWriter(covFile, append)) {
            writer.write("Sample\tGene\tChr\tStart\tEnd\tTag\tLength\tMeanDepth\n");
            writer.flush();
            for (Gene gene : sqrlist) {
                StringBuilder str = new StringBuilder();
                str.append(gene.getSample()).append("\t");
                str.append(gene.getName()).append("\t");
                str.append(gene.getChr()).append("\t");
                str.append(gene.getStart()).append("\t");
                str.append(gene.getEnd()).append("\t");
                str.append(gene.getTag()).append("\t");
                str.append(gene.getLen()).append("\t");
                str.append(String.format("%.2f", gene.getDepth())).append("\n");
                writer.write(str.toString());
                writer.flush();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void help(Options options) {
        HelpFormatter formater = new HelpFormatter();
        formater.setOptionComparator(null);
        formater.printHelp(150,
                Seq2c.class.getSimpleName() + " sample2bam.txt regions.bed coverage.txt [sample_name1[:sample_name2]] [-options]\n",
                "Arguments are:\n"
                        + "sample2bam.txt:      Required. A file containing list of genes and bam files, each sample at new line.  At least two columns. \nFirst is the sample name, 2nd is the bam file name.\n"
                        + "regions.bed:         Required. A bed file with regions of interest. At least 4 columns.\n"
                        + "coverage.txt:        Required. A file for writing coverage output.\n"
                        + "sample_name:         Control sample names. For multiple controls, separate them using ':'.\n\n",
                options, "");
        System.exit(0);
    }

    static Options buildCmdOptions() {
        Options options = Lr2gene.getOptions();
        options.addOption(OptionBuilder.withArgName("number of threads")
                .hasOptionalArg()
                .withType(Number.class)
                .withDescription("Manages the number of threads that do the work.\n" +
                        "        If this parameter is missing, then the mode is one-thread.\n" +
                        "        If you add the `-i` parameter, the number of threads equals to the number of processor cores.\n" +
                        "        The parameter '-i threads' sets the number of threads explicitly.")
                .isRequired(false)
                .create("i"));
        options.addOption(OptionBuilder.withArgName("run part")
                .hasOptionalArg()
                .withType(Number.class)
                .withDescription("Manages the separate launch of two parts of code.\n" +
                        "        For launching the first  part of code option '–r 1' (analogue to seq2cov.pl script in Perl) is used.\n" +
                        "        For launching the second part of code option '–r 2' (analogue to bam2reads.pl, cov2lr.pl and lr2gene.pl scripts in Perl).\n" +
                        "        Without option '–r' all code will run from start to end.")
                .isRequired(false)
                .create("r"));
        options.addOption(OptionBuilder.withArgName("help")
                .hasOptionalArg()
                .withType(Number.class)
                .withDescription("Prints usage to console.")
                .isRequired(false)
                .create("h"));
        return options;
    }

    static int getThreadsCount(CommandLine cmd) throws ParseException {
        int threads = 0;
        if (cmd.hasOption("i")) {
            Object value = cmd.getParsedOptionValue("i");
            if (value == null) {
                threads = Runtime.getRuntime().availableProcessors();
            } else {
                threads = ((Number)value).intValue();
            }
        }
        return threads;

    }

    static int getRunPart(CommandLine cmd) throws ParseException  {
        int part = 0;
        if (cmd.hasOption("r")) {
            Object value = cmd.getParsedOptionValue("r");
            if (value != null) {
                part = ((Number)value).intValue();
            }
        }
        return part;
    }
}
