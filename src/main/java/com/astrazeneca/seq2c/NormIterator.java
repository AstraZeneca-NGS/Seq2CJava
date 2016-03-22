package com.astrazeneca.seq2c;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

/**
 * Created by Mariia_Zueva on 3/8/2016.
 */
public class NormIterator {

    BufferedReader bufferReader;
    Normalization next;

    public NormIterator(String fileName) {
        try {
            bufferReader = new BufferedReader(new FileReader(fileName));
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void close() {
        try {
            bufferReader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public Normalization nextNorm() {

        while (true) {
            Normalization norm;
            if (next == null) {
                norm = readNextNorm();
            } else {
                norm = next;
            }
            next = readNextNorm();
            while (next != null && norm.equals(next)) {
                addSamples(norm, next);
                next = readNextNorm();
            }
            return norm;
        }
    }

    private void addSamples(Normalization norm, Normalization next) {
        norm.samples.putAll(next.samples);
    }

    private Normalization readNextNorm() {
        String line = getNextLine();
        if (line == null) return null;
        next = readNorm(line);
        return next;
    }

    private Normalization readNorm(String line) {
        if (line == null) return null;

        String[] sampleLines = line.split("\\t");
        if (sampleLines.length != 3) return null;

        String gene = sampleLines[0];
        String sample = sampleLines[1];
        double norm1 = Double.parseDouble(sampleLines[2]);

        Normalization norm = new Normalization();
        norm.gene = gene;
        HashMap<String, Double> samples = new HashMap<>();
        samples.put(sample, norm1);
        norm.samples = samples;
        return norm;
    }


    private String getNextLine() {

        String line = null;
        try {
            while (true) {
                line = bufferReader.readLine();
                if (line == null) {
                    return null;
                }
                if (line.isEmpty()) {
                    return null;
                }
                return line;
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
        return line;
    }


}
