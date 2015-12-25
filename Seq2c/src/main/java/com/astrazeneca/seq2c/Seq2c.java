package com.astrazeneca.seq2c;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author Petr_Rastegaev
 */
public class Seq2c {

    public static void main(String[] args) throws Exception {

        String sam2bamFile = args[0];
        final String bedFile = args[1];
        String control = "";
        String seq2copt = "";

        if (args.length > 2) {
            if(!args[2].startsWith("-")) {
                control = args[2];
            }
            StringBuilder builder = new StringBuilder();
            for (int i = 2; i < args.length; i++) {
                if (builder.length() > 0) {
                    builder.append(' ');
                }
                builder.append(args[i]);
            }
            seq2copt = builder.toString();
        }

        final boolean controlFlag = control.isEmpty() ? false : true;

        final int threadsCnt = getThreadsCount(seq2copt);

        Dispatcher.init(threadsCnt);
        try {

            final Map<String, String> sam2bam = Bam2Reads.parseFile(sam2bamFile);

            final Collection<Gene> sqrlist = new ArrayList<>();
            for (final Map.Entry<String, String> entry : sam2bam.entrySet()) {
                Seq2cov sec2cov = new Seq2cov(bedFile, entry.getValue(), entry.getKey());
                Collection<Gene> cov = sec2cov.process();
                sqrlist.addAll(cov);
            }

            Map<String, Long> stat = Bam2Reads.printStatsToFile(sam2bamFile);

            Cov2lr cov2lr = new Cov2lr(true, stat, sqrlist, controlFlag, control);
            List<Sample> cov = cov2lr.doWork();

            new Lr2gene(cov, seq2copt, controlFlag).run();

        } finally {
            Dispatcher.shutdown();
        }

    }

    private static final Pattern threadsOpt = Pattern.compile("-i\\s*(\\d+)?");
    private static int getThreadsCount(String opts) {
        if(opts.isEmpty()) {
            return 0;
        }

        Matcher matcher = threadsOpt.matcher(opts);
        if (matcher.find()) {
            String num = matcher.group(1);
            if (num == null) {
                return Runtime.getRuntime().availableProcessors();
            }
            return Integer.parseInt(num);
        }
        return 0;

    }

}
