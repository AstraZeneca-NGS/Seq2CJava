package com.astrazeneca.seq2c;

public class Sample {

    private String name;
    private String sample;
    private String gene;
    private String chr;
    private long len;
    private long start;
    private long end;
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

    public Sample(String name, String sample, String chr, long start, long end, String gene, long len, double cov) {
        this.name = name;
        this.sample = sample;
        this.start = start;
        this.end = end;
        this.len = len;
        this.cov = cov;
        this.gene = gene;
        this.chr = chr;
    }

    public String getKey() {
        return name + ":" + sample;
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

    public String getSample() {
        return sample;
    }

    public String getName() {
        return name;
    }

    public String getGene() {
        return gene;
    }

    public double getFactor2() {
        return factor2;
    }

    public void setFactor2(double factor2) {
        this.factor2 = factor2;
    }

    public long getLen() {
        return len;
    }

    public void addLen(long len) {
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

    public double getNorm1() {
        return norm1;
    }

    public double getNorm1b() {
        return norm1b;
    }

    public void setNorm1b(double norm1b) {
        this.norm1b = norm1b;
    }

    public double getNorm2() {
        return norm2;
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

    public void setNorm1(double norm1) {
        this.norm1 = norm1;
    }

    public String getResultString(String title) {
        String format = "%.2f%n";
        StringBuilder builder = new StringBuilder();
        builder.append(sample).append("\t");
        builder.append(title).append("\t");
        builder.append(cov).append("\t");
        builder.append(String.format(format, norm1).trim()).append("\t");
        builder.append(String.format(format, norm1b).trim()).append("\t");
        builder.append(String.format(format, norm2).trim()).append("\t");
        builder.append(String.format(format, norm3).trim()).append("\t");

        return builder.toString();
    }

    public String getTitle(Gene gene) {
        StringBuilder builder = new StringBuilder();
        builder.append(gene.getName()).append("\t");
        builder.append(gene.getChr()).append("\t");
        builder.append(gene.getStart()).append("\t");
        builder.append(gene.getEnd()).append("\t");
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
     * @param chr
     *            the chr to set
     */
    public void setChr(String chr) {
        this.chr = chr;
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
