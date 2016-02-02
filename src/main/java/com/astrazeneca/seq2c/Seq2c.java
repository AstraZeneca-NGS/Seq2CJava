package com.astrazeneca.seq2c;

import org.apache.commons.cli.*;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;

public class Seq2c {

    public static void main(String[] args) throws Exception {

        CommandLine cmd = new BasicParser().parse(buildCmdOptions(), args);

        String[] cmdArgs = cmd.getArgs();
        String sam2bamFile = cmdArgs[0];
        final String bedFile = cmdArgs[1];
        final String covFile = cmdArgs[2];
        String control = "";

        if (cmdArgs.length > 2) {
            control = cmdArgs[2];
        }

        Dispatcher.init(getThreadsCount(cmd));
        try {

            final Map<String, String> sam2bam = Bam2Reads.parseFile(sam2bamFile);

            final Collection<Gene> sqrlist = new ArrayList<>();
            boolean rewrite = false;
            for (final Map.Entry<String, String> entry : sam2bam.entrySet()) {
                Seq2cov sec2cov = new Seq2cov(bedFile, entry.getValue(), entry.getKey());
                Collection<Gene> cov = sec2cov.process();
                //print coverage to the file
                printCoverage(cov, covFile, rewrite);
                rewrite = true;
                sqrlist.addAll(cov);
            }


            Map<String, Long> stat = Bam2Reads.printStatsToFile(sam2bamFile);

            Cov2lr cov2lr = new Cov2lr(true, stat, sqrlist, control);
            List<Sample> cov = cov2lr.doWork();

            Lr2gene lr2gene = new Lr2gene(cov);
            lr2gene.init(cmd);
            lr2gene.setUseControl(cov2lr.isUseControlSamples());
            lr2gene.process();

        } finally {
            Dispatcher.shutdown();
        }
    }

    private static void printCoverage(Collection<Gene> sqrlist, String covFile, boolean rewrite) {
        try (FileWriter writer = new FileWriter(covFile, rewrite)) {
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

    private static Options buildCmdOptions() {
        Options options = Lr2gene.getOptions();
        options.addOption(OptionBuilder.withArgName("number of threads")
                .hasOptionalArg()
                .withType(Number.class)
                .isRequired(false)
                .create("i"));
        return options;
    }

    private static int getThreadsCount(CommandLine cmd) throws ParseException {
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

}
