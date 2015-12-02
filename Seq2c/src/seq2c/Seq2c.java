/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package seq2c;
import java.util.*;
/**
 *
 * @author Petr_Rastegaev
 */
public class Seq2c {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        String sam2bamFile =  args[0]; // "D:\\Samles\\seq2c\\sample2bam1.txt";
        String bedFile = args[1];//"D:\\Samles\\seq2c\\Illumina_TruSeq_Exome.bed";
        String control = "";
        String seq2copt = "";
        
        if(args.length > 2){
            control = args[2];
        }
        if(args.length > 3){
            seq2copt = args[3];
        }
        
        boolean controlFlag = false;
        if(!control.isEmpty()) controlFlag = true;
        //Seq2cov seq2cov = new Seq2cov("D:\\VarDict\\AURA\\GeneRead_Lung.primers.bed","D:\\VarDict\\AURA\\AURA_26-ready.bam","AURA_26");
         Map<String, String> sam2bam = Bam2Reads.parseFile(sam2bamFile);
        ArrayList<Gene> sqrlist = new ArrayList();
        for(Map.Entry<String,String> entry: sam2bam.entrySet()){
            Seq2cov seq2cov = new Seq2cov(bedFile,entry.getValue(),entry.getKey());
            ArrayList<Gene> sq = seq2cov.process();
            sqrlist.addAll(sq);
        } 
         Map<String, Long> stat = Bam2Reads.printStatsToFile(sam2bamFile);
         Cov2lr cov2lr = new Cov2lr(true,stat,sqrlist,controlFlag,control);
         ArrayList<Sample> cov = cov2lr.doWork();
         Lr2gene lr2 = new Lr2gene(cov,seq2copt);
         lr2.run();
        //System.out.println(Arrays.asList(result));
    }
    
}
