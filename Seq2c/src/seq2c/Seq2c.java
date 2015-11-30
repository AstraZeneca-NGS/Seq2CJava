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
        Seq2cov seq2cov = new Seq2cov("D:\\Samles\\seq2c\\Illumina_TruSeq_Exome.bed","D:\\Samles\\seq2c\\Patient-E7801004-ff-ready.bam","Patient-E7801004-ff");
        //Seq2cov seq2cov = new Seq2cov("D:\\VarDict\\AURA\\GeneRead_Lung.primers.bed","D:\\VarDict\\AURA\\AURA_26-ready.bam","AURA_26");
         ArrayList<Gene> sqrlist = seq2cov.process();
         for(Gene sqr : sqrlist){
             String out = sqr.toString();
             //System.out.println(out);
         }
         Map<String, Long> stat = Bam2Reads.printStatsToFile("D:\\Samles\\seq2c\\sample2bam1.txt");
         Cov2lr cov2lr = new Cov2lr(true,stat,sqrlist,false,"");
         cov2lr.doWork();
        //System.out.println(Arrays.asList(result));
    }
    
}
