package com.astrazeneca.seq2c;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
        for(Map.Entry<String, List<Sample>> entry : next.genes.entrySet()) {
            String gene = entry.getKey();
            List<Sample> list;
            if (norm.genes.containsKey(gene)) {
                list = norm.genes.get(gene);
            } else {
                list = new ArrayList<>();
                norm.genes.put(gene, list);
            }
            list.addAll(entry.getValue());
        }
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
        if (sampleLines.length < 8) return null;

        String sample = sampleLines[0];
        String gene = sampleLines[1];
        Sample sampleObj = readSample(line);

        Normalization norm = new Normalization();
        norm.sample = sample;
        HashMap<String, List<Sample>> genes = new HashMap<>();
        List<Sample> list = new ArrayList<>();
        list.add(sampleObj);
        genes.put(gene, list);
        norm.genes = genes;
        return norm;
    }

    private Sample readSample(String currentLine) {

        if (currentLine == null) return null;

        String[] sampleLines = currentLine.split("\\s+");
        if (sampleLines.length < 8) return null;

        String sampleName = sampleLines[0];
        String gene = sampleLines[1];
        String chr = sampleLines[2];
        int start = Integer.parseInt(sampleLines[3]);
        int end = Integer.parseInt(sampleLines[4]);
        int len = Integer.parseInt(sampleLines[5]);

        String key = gene;
        Sample sample = new Sample(key, sampleName, chr, start, end, gene, len, 0.0);
        sample.setNorm3(Double.parseDouble(sampleLines[6]));
        sample.setNorm3s(Double.parseDouble(sampleLines[7]));
        return sample;
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
