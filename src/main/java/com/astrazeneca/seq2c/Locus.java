package com.astrazeneca.seq2c;


class Locus {
    String geneName;
    String chr;
    long start;
    long end;
    long length;

    Locus(String geneName, long start, long end, String chr) {
        this.geneName = geneName;
        this.chr = chr;
        this.start = start;
        this.end = end;
    }

    void setLength(long length) {
        this.length = length;
    }

    void addLength(long length) {
        this.length += length;
    }

    public void shiftStart(long start) {
        if (this.start > start) {
            this.start = start;
        }
    }

    public StringBuilder getName() {
        StringBuilder result = new StringBuilder();
        result.append(chr).append("\t");
        result.append(start).append("\t");
        result.append(end).append("\t");
        result.append(length);
        return result;
    }

    public void shiftEnd(long end) {
        if (this.end < end) {
            this.end = end;
        }
    }
}
