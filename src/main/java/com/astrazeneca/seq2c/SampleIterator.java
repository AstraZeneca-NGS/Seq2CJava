package com.astrazeneca.seq2c;

import org.apache.commons.math3.util.Precision;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Iterator;

public class SampleIterator implements Iterator<Sample> {

    BufferedReader reader;
    String currentLine;
    boolean amplicon;

    public SampleIterator(String fileName, boolean amplicon) {
        this.amplicon = amplicon;
        try {
            this.reader = new BufferedReader(new FileReader(fileName));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }


    @Override
    public boolean hasNext() {
        try {
            if (nextLine() != null) {
                return true;
            } else {
                reader.close();
                return false;
            }
        } catch (IOException e) {
            e.printStackTrace();
            return false;
        }
    }

    private String nextLine() {
        try {
            currentLine = reader.readLine();
        } catch (IOException e) {
            currentLine = null;
            e.printStackTrace();
        }
        return currentLine;
    }

    @Override
    public Sample next() {
        if (currentLine == null) {
            throw new IllegalStateException("Reader has no more lines");
        }

        //skip lines with words "Sample, Whole, Control, Undertermined
        while (currentLine != null && currentLine.contains("Sample") || currentLine.contains("Whole") || currentLine.contains("Control") ||
                currentLine.contains("Undertermined")) {
            currentLine = nextLine();
        }

        String[] sampleLines = currentLine.split("\\s+");
        if (sampleLines.length == 8) {

            String sample = sampleLines[0];
            String gene = sampleLines[1];
            String chr = sampleLines[2];
            long start = Long.parseLong(sampleLines[3]);
            long end = Long.parseLong(sampleLines[4]);
            long len = Long.parseLong(sampleLines[6]);
            double depth = Double.parseDouble(sampleLines[7]);
            String key = amplicon ? join(" ", gene, chr, Long.toString(start), Long.toString(end), Long.toString(len)) : gene;
            Sample sampleObj = getSample(key, sample, chr, start, end, len, Precision.round(depth, 2), gene);
        }
        return null;
    }

    private Sample getSample(String key, String sample, String chr, Long start, long end, long len, double depth, String gene) {

        Sample sampleObj = new Sample(key, sample, chr, start, end, gene, len, depth);
        return sampleObj;
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("Cannot renove from file");
    }

    private static String join(String delim, Object... array) {
        StringBuilder builder = new StringBuilder();
        for (Object s : array) {
            if (builder.length() > 0) {
                builder.append(delim);
            }
            builder.append(s);
        }
        return builder.toString();
    }
}
