/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package seq2c;
import htsjdk.samtools.AbstractBAMFileIndex;
import htsjdk.samtools.BAMIndexMetaData;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReaderFactory;
import java.io.*;
import java.util.*;
import java.util.regex.*;
//import java.math.*;
/**
 *
 * @author Petr_Rastegaev
 */
public class Seq2cov {
    //public static void main(String[] args) 
    //{
       //processBed(fileName);
     //}
    boolean PCRamplbc = false;
    private String bedfileName;
    private String bamfileName;
    private String sampleName;
    private LinkedHashMap<String,Seq2cov_gene> result_bed = new LinkedHashMap();
    private ArrayList<Gene> result_bam = new ArrayList();
    public Seq2cov(String bename, String bamname, String spname){
         bedfileName = bename;
         bamfileName = bamname;
         sampleName = spname;
    }
    public ArrayList<Gene> process(){
        processBed();
        ProcessBam();
        return result_bam;
    }
    
    private LinkedHashMap<String,Seq2cov_gene> processBed(){
    //System.out.println("Reading File from Java code");
       //Name of the file
       try{
          //Create object of FileReader
          FileReader inputFile = new FileReader(bedfileName);
           //Variable to hold the one line data
          try ( //Instantiate the BufferedReader Class
                   BufferedReader bufferReader = new BufferedReader(inputFile)) {
               //Variable to hold the one line data
               String line;
               // Read file line by line and print on the console
               while ((line = bufferReader.readLine()) != null) {
                String[] parts = line.split("\t");
                String chr = parts[0];
                String gene = parts[3];
                int start = Integer.parseInt(parts[1]);
                int end = Integer.parseInt(parts[2]);
                if(parts.length == 8){ 
                    PCRamplbc = true;
                    //System.out.println("Its opt_a!");
                }
                start += 1;
                Integer[]CDS = {start,end};
                Seq2cov_gene gn;
                if(result_bed.containsKey(gene)){
                    gn =  result_bed.get(gene);
                    gn.addCDS(CDS);
                }else{
                    gn = new Seq2cov_gene(chr,gene,CDS);    
                }
                result_bed.put(gene,gn);
                
               }
               //Close the buffer reader
           }
       }catch(Exception e){
          System.out.println("Error :" + e.getMessage() +e.getStackTrace()[0]);
          e.printStackTrace();
       }
       return result_bed;
    }
    
    private ArrayList ProcessBam(){
    try (SamReader sam = SamReaderFactory.makeDefault().open(new File(bamfileName))) {
                    SamReader.Indexing ind = sam.indexing();
                    AbstractBAMFileIndex index = (AbstractBAMFileIndex) ind.getIndex();
                    String bamhdr = sam.getFileHeader().getTextHeader();
                    boolean genome;
                    String regexp = "SN:chr"; //ignore supplementary alignments
                    Pattern pattern = Pattern.compile(regexp);
                    Matcher matcher = pattern.matcher(bamhdr);
                    if(matcher.find()){
                        genome = true; //"hg";
                    } else {
                        genome = false; // "grch";
                    }
     
    for (Seq2cov_gene gn : result_bed.values()) {
        //System.out.println( tchr + " " + gn.CDS[0] + " " + gn.CDS[1]);
            int exoncov;
            int total = 0;
            int geneLength = 0;
            int geneStart = 500000000;
            int geneEnd = 0;
            String tchr = gn.chr;
            regexp = "^chr"; 
            pattern = Pattern.compile(regexp);
            matcher = pattern.matcher(tchr);
            if(genome && !matcher.find()){
                tchr += "chr"; 
            }else if(!genome && matcher.find()){
                tchr = tchr.replace("chr", "");
            }
            for(Integer[] CDS_arr : gn.CDS){
                exoncov = 0;
                int cdsStart = CDS_arr[0];
                int cdsEnd = CDS_arr[1];
                geneLength += cdsEnd - cdsStart + 1;
                if(cdsStart < geneStart){
                    geneStart = cdsStart;
                }
                if(cdsEnd > geneEnd){
                    geneEnd = cdsEnd;
                }
                SAMRecordIterator iterator = sam.query(tchr,CDS_arr[0],CDS_arr[1],false);
                while (iterator.hasNext()) {
                    SAMRecord rec = iterator.next();
                    String dir;
                    String samstring = rec.getSAMString();
                    regexp = "\\tSA:Z:"; //ignore supplementary alignments
                    pattern = Pattern.compile(regexp);
                    matcher = pattern.matcher(samstring);
                    String [] samstringFields = samstring.split("\t");
                    int bt = Integer.parseInt(samstringFields[1]);
                    int btw1 = bt & 0x800;
                    if(matcher.find() && btw1 != 0){
                        continue;
                    }
                    int btw2 = bt & 0x10;
                    if(btw2 != 0){
                        dir = "-";
                    }else{
                        dir = "+";
                    }
                    int start = Integer.parseInt(samstringFields[3]);
                    regexp = "(\\d+)[MD]";
                    pattern = Pattern.compile(regexp);
                    matcher = pattern.matcher(samstringFields[5]);
                    List<Integer> seqs = new ArrayList<Integer>();
                    while (matcher.find()) {
                        //System.out.println(matcher.group(1));
                        seqs.add(Integer.parseInt(matcher.group(1)));
                    }
                    int end = start - 1;
                    for (Integer seq: seqs){
                        end += seq;
                    }
                    //System.out.println(gn.gene + "     " + start + "     " + end );
                    if(PCRamplbc){ //$opt_a
                        double dis = 10;
                        double ovlp = 0.95;
                        regexp = "^(\\d+)S";
                        pattern = Pattern.compile(regexp);
                        matcher = pattern.matcher(samstringFields[5]);
                        regexp = "(\\d+)S$";
                        Pattern pattern2 = Pattern.compile(regexp);
                        Matcher matcher2 = pattern2.matcher(samstringFields[5]);
                        int ts1;
                        int te1;
                        int seqStart = start;
                        int seqEnd = end;
                        if(matcher.find() && "+".equals(dir) ){
                            int add = Integer.parseInt(matcher.group(1));
                            if(seqStart - add > cdsStart){
                                ts1 = seqStart - add;
                            } else {
                                ts1 = cdsStart;
                            }
                            if(seqEnd < cdsEnd){
                                te1 = seqEnd;
                            } else {
                                te1 = cdsEnd;
                            }
                            //System.out.println("entered if " +cdsStart+" "+cdsEnd + tchr);
                            if(!(Math.abs((double)ts1-(double)te1)/((double)seqEnd - (double)seqStart + (double)add) > ovlp)){
                                continue;
                            }
                                    
                        }else if(matcher2.find() && "-".equals(dir)){
                            int add = Integer.parseInt(matcher2.group(1));
                            if(seqStart > cdsStart){
                                ts1 = seqStart;
                            } else {
                                ts1 = cdsStart;
                            }
                            if(seqEnd + add < cdsEnd){
                                te1 = seqEnd + add;
                            } else {
                                te1 = cdsEnd;
                            }
                            //System.out.println("entered elseif " +cdsStart+" "+cdsEnd + tchr);
                            if(!(Math.abs((double)te1-(double)ts1)/((double)seqEnd  + (double)add - (double)seqStart) > ovlp)){
                                continue;
                            }
                        }else{
                            if("=".equals(samstringFields[6]) && !samstringFields[8].isEmpty() && !"0".equals(samstringFields[8])){
                                int pNext = Integer.parseInt(samstringFields[8]);
                                int rNext = Integer.parseInt(samstringFields[7]);
                                if(pNext > 0){
                                    seqStart = seqStart;
                                    seqEnd = seqStart + pNext - 1;
                                } else {
                                   seqStart = rNext;
                                   seqEnd = rNext - pNext - 1; 
                                }
                            }
                                if(seqStart > cdsStart){
                                    ts1 = seqStart;
                                } else {
                                    ts1 = cdsStart;
                                }
                                if(seqEnd < cdsEnd){
                                    te1 = seqEnd;
                                } else {
                                    te1 = cdsEnd;
                                }
                                //System.out.println("entered else " +cdsStart+" "+cdsEnd + tchr);
                                if(!((Math.abs((double)seqStart-(double)cdsStart) <= dis && Math.abs((double)seqEnd-(double)cdsEnd) <= dis) && Math.abs(((double)ts1-(double)te1)/((double)seqEnd-(double)seqStart)) > ovlp)){
                                    continue;
                                }
                            
                        }
                        exoncov++;
                        total++;
                    }else{
                        double alEnd;
                        double alStart;
                        if(cdsEnd > end){
                            alEnd = end;
                        }else{
                            alEnd = cdsEnd;
                        }
                        
                        if(cdsStart > start){
                            alStart = cdsStart;
                        }else{
                            alStart = start;
                        }
                        
                        double alignen = (alEnd - alStart) +1;
                        exoncov += alignen;
                        total += alignen;
                    }

                } //iterator
                iterator.close();
                //System.out.println(exoncov + " " + total);
                Gene sqr;
                if(PCRamplbc){
                    sqr = new Gene(sampleName,gn.gene,tchr,cdsStart,cdsEnd,"Amplicon",cdsEnd - cdsStart + 1,exoncov);
                }else{
                    sqr = new Gene(sampleName,gn.gene,tchr,cdsStart,cdsEnd,"Amplicon",cdsEnd - cdsStart + 1,(double)exoncov/(double)(cdsEnd - cdsStart + 1));
                }
                result_bam.add(sqr);
            
            }//for loop2
            Gene sqr;
            if(PCRamplbc){
                sqr = new Gene(sampleName,gn.gene,tchr,geneStart,geneEnd,"Whole-Gene",geneLength,total);
            }else{
                sqr = new Gene(sampleName,gn.gene,tchr,geneStart,geneEnd,"Whole-Gene",geneLength,(double)total/geneLength);
            }
            result_bam.add(sqr);
        } //for loop1
    } catch (IOException e) {
                    System.err.println("Cannot read file " + bamfileName);
                }               
        //
    return result_bam;
    }
}

