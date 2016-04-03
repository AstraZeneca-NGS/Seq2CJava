package com.astrazeneca.seq2c;

import java.io.Serializable;

public class Sample implements Serializable{

    private String name;
    private String sample;
    private String gene;
    private String chr;
    private int len;
    private int start;
    private int end;
    private double cov;
    private double norm1;
    private double factor2;
    private double norm1b;
    private double norm2;

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Sample sample1 = (Sample) o;

        if (name != null ? !name.equals(sample1.name) : sample1.name != null) return false;
        return sample != null ? sample.equals(sample1.sample) : sample1.sample == null;

    }

    @Override
    public int hashCode() {
        int result = name != null ? name.hashCode() : 0;
        result = 31 * result + (sample != null ? sample.hashCode() : 0);
        return result;
    }

    private double norm3;
    private double norm3s;

    public Sample(String name, String sample, String chr, int start, int end, String gene, int len, double cov) {
        this.name = name;
        this.sample = sample;
        this.start = start;
        this.end = end;
        this.len = len;
        this.cov = cov;
        this.gene = gene;
        this.chr = chr;
    }


    public long getStart() {
        return this.start;
    }

    public long getEnd() {
        return end;
    }

    public String getSample() {
        return sample;
    }

    public String getName() {
        return name;
    }

    public String getGene() {
        return gene;
    }

    public int getLen() {
        return len;
    }

    public void addLen(int len) {
        this.len = this.len + len;
    }

    public double getCov() {
        return cov;
    }

    public void addCov(double cov) {
        this.cov = this.cov + cov;
    }

    @Override
    public String toString() {
        return sample + " " + len + " " + cov;
    }

    public double getNorm1b() {
        return norm1b;
    }

    public void setNorm1b(double norm1b) {
        this.norm1b = norm1b;
    }

    public void setNorm2(double norm2) {
        this.norm2 = norm2;
    }

    public double getNorm3() {
        return norm3;
    }

    public void setNorm3(double norm3) {
        this.norm3 = norm3;
    }

    public String getResultString() {
        String format = "%.2f%n";
        StringBuilder builder = new StringBuilder();
        builder.append(sample).append("\t");
        builder.append(getTitle()).append("\t");
        builder.append(String.format(format, norm3).trim()).append("\t");
        builder.append(String.format(format, norm3s).trim()).append("\t");

        return builder.toString();
    }

    public String getTitle() {
        StringBuilder builder = new StringBuilder();
        builder.append(gene).append("\t");
        builder.append(chr).append("\t");
        builder.append(start).append("\t");
        builder.append(end).append("\t");
        builder.append(len).append("\t");
        return builder.toString();
    }

    /**
     * @return the chr
     */
    public String getChr() {
        return chr;
    }

    /**
     * @return the norm3s
     */
    public double getNorm3s() {
        return norm3s;
    }

    /**
     * @param norm3s
     *            the norm3s to set
     */
    public void setNorm3s(double norm3s) {
        this.norm3s = norm3s;
    }

}
