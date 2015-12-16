package seq2c;

import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
/**
 *
 * @author Petr_Rastegaev
 */
public class Seq2c {

    /**
     * @param args the command line arguments
     * @throws IOException
     * @throws InterruptedException
     */
    public static void main(String[] args) throws IOException, InterruptedException {

        String sam2bamFile =  args[0]; // "D:\\Samles\\seq2c\\sample2bam1.txt";
        final String bedFile = args[1];//"D:\\Samles\\seq2c\\Illumina_TruSeq_Exome.bed";
        String control = "";
        String seq2copt = "";

        if(args.length > 2){
            control = args[2];
        }

        if(args.length > 3){
            seq2copt = args[3];
        }

        boolean controlFlag = control.isEmpty() ? false : true;

        //Seq2cov seq2cov = new Seq2cov("D:\\VarDict\\AURA\\GeneRead_Lung.primers.bed","D:\\VarDict\\AURA\\AURA_26-ready.bam","AURA_26");


//        int threads = 1;

//        ExecutorService es = Executors.newFixedThreadPool(threads);

        Map<String, String> sam2bam = Bam2Reads.parseFile(sam2bamFile);

        final Collection<Gene> sqrlist = new ConcurrentLinkedQueue<>();

        long time = System.currentTimeMillis();
        for (final Map.Entry<String, String> entry : sam2bam.entrySet()) {
//            es.submit(new Runnable() {
//                @Override
//                public void run() {
//                    try {
                        sqrlist.addAll(new Seq2cov(bedFile, entry.getValue(), entry.getKey()).process());
//                    } catch (IOException e) {
//                        e.printStackTrace();
//                        throw new RuntimeException(e);
//                    }
//                }
//            });
        }

//        es.shutdown();
//        es.awaitTermination(1, TimeUnit.DAYS);

//        System.err.println("Time: " + (System.currentTimeMillis() - time));


        Map<String, Long> stat = Bam2Reads.printStatsToFile(sam2bamFile);

        Cov2lr cov2lr = new Cov2lr(true, stat, sqrlist, controlFlag, control);

        List<Sample> cov = cov2lr.doWork();

        new Lr2gene(cov, seq2copt, controlFlag).run();

    }

}
