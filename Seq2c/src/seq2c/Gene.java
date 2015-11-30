/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package seq2c;

/**
 * Created by Mariia_Zueva on 11/25/2015.
 */
public class Gene {
    private String sample;
    private String name;
    private String chr;
    private long start;
    private long end;
    private long len;
    private double mdepth;
    public String tag;

    public Gene(String sample, String name, String chr, long start, long end, String tag, long len, double mdepth) {
        this.sample = sample;
        this.name = name;
        this.chr = chr;
        this.start = start;
        this.tag = tag;
        this.end = end;
        this.len = len;
        this.mdepth = mdepth;
    }
    public double getDepth(){
        return this.mdepth;
    }
    public String getSample(){
        return this.sample;
    }
    
    public String getTag(){
        return this.tag;
    }
    
     
    public long getStart() {
        return this.start;
    }

    public void setStart(long start) {
        this.start = start;
    }

    public long getEnd() {
        return end;
    }

    public void setEnd(long end) {
        this.end = end;
    }

    public long getLen() {
        return len;
    }

    public String getName() {
        return this.name;
    }

    public String getChr() {
        return this.chr;
    }

    public void setLen(long len) {
        this.len = len;
    }
    
    public void setMdepth(long mdepth) {
        this.mdepth = mdepth;
    }

    @Override
    public String toString() {
        return sample + "\t" + name + "\t" + chr + "\t" + start + "\t" + end + "\t" + tag + "\t" + len + "\t" + String.format("%.2f",mdepth);
    }
}

