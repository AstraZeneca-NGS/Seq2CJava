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
public class Seq2cov_gene {
    public ArrayList<Integer[]> CDS = new ArrayList();
    public String chr;
    public String gene;
    public Seq2cov_gene(String chr_, String gene_, Integer[] CDS_) {
        chr = chr_;
        gene = gene_;
        CDS.add(CDS_);
    }
    public void addCDS(Integer[] CDS_){
      CDS.add(CDS_);  
    }
}
