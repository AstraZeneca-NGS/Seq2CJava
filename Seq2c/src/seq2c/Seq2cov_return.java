/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package seq2c;

import java.util.ArrayList;

/**
 *
 * @author Petr_Rastegaev
 */
public class Seq2cov_return {
    public ArrayList<Integer[]> CDS = new ArrayList();
    public String Sample;
    public String Gene;
    public String Chr;
    public int Start;
    public int End;
    public String Tag;
    public int Length;
    public double MeanDepth;        
    public Seq2cov_return(String sample_, String gene_, String chr_, int start_, int end_, String tag_, int length_, double mdepth_) {
        Sample = sample_;
        Gene = gene_;
        Chr = chr_;
        Start = start_;
        End = end_;
        Tag = tag_;
        Length = length_;
        MeanDepth = mdepth_;
    }
    public String toCovString(){
        String str = Sample + "\t" + Gene + "\t" + Chr + "\t" + Start + "\t" + End + "\t" + Tag + "\t" + Length + "\t" + String.format("%.2f",MeanDepth);
        return str;
    }
}
