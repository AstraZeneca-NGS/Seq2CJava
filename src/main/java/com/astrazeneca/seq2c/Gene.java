package com.astrazeneca.seq2c;

/**
 * Data structure for gene.
 * */
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

    public double getDepth() {
        return this.mdepth;
    }

    public String getSample() {
        return this.sample;
    }

    public String getTag() {
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
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Gene gene = (Gene) o;

        if (getStart() != gene.getStart()) return false;
        if (getEnd() != gene.getEnd()) return false;
        if (getLen() != gene.getLen()) return false;
        if (Double.compare(gene.mdepth, mdepth) != 0) return false;
        if (getSample() != null ? !getSample().equals(gene.getSample()) : gene.getSample() != null)
            return false;
        if (getName() != null ? !getName().equals(gene.getName()) : gene.getName() != null)
            return false;
        if (getChr() != null ? !getChr().equals(gene.getChr()) : gene.getChr() != null)
            return false;
        return getTag() != null ? getTag().equals(gene.getTag()) : gene.getTag() == null;
    }

    @Override
    public int hashCode() {
        int result;
        long temp;
        result = getSample() != null ? getSample().hashCode() : 0;
        result = 31 * result + (getName() != null ? getName().hashCode() : 0);
        result = 31 * result + (getChr() != null ? getChr().hashCode() : 0);
        result = 31 * result + (int) (getStart() ^ (getStart() >>> 32));
        result = 31 * result + (int) (getEnd() ^ (getEnd() >>> 32));
        result = 31 * result + (int) (getLen() ^ (getLen() >>> 32));
        temp = Double.doubleToLongBits(mdepth);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        result = 31 * result + (getTag() != null ? getTag().hashCode() : 0);
        return result;
    }

    @Override
    public String toString() {
        return sample + "\t" + name + "\t" + chr + "\t" + start + "\t" + end + "\t" + tag + "\t" + len + "\t" + String.format("%.2f", mdepth);
    }
}
