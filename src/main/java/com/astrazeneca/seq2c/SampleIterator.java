package com.astrazeneca.seq2c;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

/**
 * Created by Mariia_Zueva on 3/8/2016.
 */
public class SampleIterator {

    boolean amplicon;
    BufferedReader bufferReader;
    Sample next;

    public SampleIterator(String fileName, boolean amplicon) {
        try {
            bufferReader = new BufferedReader(new FileReader(fileName));
            this.amplicon = amplicon;
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

    public Sample nextSample() {

        try {

            String line = bufferReader.readLine();
            while (true) {
                if (line == null) {
                    return null;
                }
                if (line.contains("Sample") || line.contains("Whole") || line.contains("Control") ||
                        line.contains("Undertermined")) {
                    line = bufferReader.readLine();
                    continue;
                }
                Sample sample;
                if (next == null) {
                    sample = readSample(line);
                } else {
                    sample = next;
                }
                next = readSample(bufferReader.readLine());
                while (next != null && sample.equals(next)) {
                    sample.addLen(next.getLen());
                    sample.addCov(next.getCov());
                    next = readSample(bufferReader.readLine());
                }
                return sample;
            }
        }catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }

    private Sample readSample(String currentLine) throws IOException {

        if (currentLine == null) return null;

        String[] sampleLines = currentLine.split("\\s+");
        if (sampleLines.length != 8) return null;

        String sampleName = sampleLines[0];
        String gene = sampleLines[1];
        String chr = sampleLines[2];
        long start = Long.parseLong(sampleLines[3]);
        long end = Long.parseLong(sampleLines[4]);
        long len = Long.parseLong(sampleLines[6]);
        double depth = Double.parseDouble(sampleLines[7]);

        String key = amplicon ? join(" ", gene, chr, Long.toString(start), Long.toString(end), Long.toString(len)) : gene;
        Sample sample = new Sample(key, sampleName, chr, start, end, gene, len, depth);
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
}
