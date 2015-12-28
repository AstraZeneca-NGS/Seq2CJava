package com.astrazeneca.seq2c;

import java.util.*;

import org.apache.commons.cli.*;

/**
 *
 * @author Petr_Rastegaev
 */
public class Seq2c {

    public static void main(String[] args) throws Exception {


        CommandLine cmd = new BasicParser().parse(buildCmdOptions(), args);

        String[] cmdArgs = cmd.getArgs();
        String sam2bamFile = cmdArgs[0];
        final String bedFile = cmdArgs[1];
        String control = "";

        if (cmdArgs.length > 2) {
            control = cmdArgs[2];
        }

        Dispatcher.init(getThreadsCount(cmd));
        try {

            final Map<String, String> sam2bam = Bam2Reads.parseFile(sam2bamFile);

            final Collection<Gene> sqrlist = new ArrayList<>();
            for (final Map.Entry<String, String> entry : sam2bam.entrySet()) {
                Seq2cov sec2cov = new Seq2cov(bedFile, entry.getValue(), entry.getKey());
                Collection<Gene> cov = sec2cov.process();
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
