package seq2c;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;

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
            control = args[2];
        }

        if (args.length > 3) {
            seq2copt = args[3];
        }

        boolean controlFlag = control.isEmpty() ? false : true;

        final int threadsCnt = Runtime.getRuntime().availableProcessors();

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

}
