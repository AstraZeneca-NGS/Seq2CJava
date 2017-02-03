package com.astrazeneca.seq2c.input;

/**
 * Construct Samples from String.
 * */
public class SampleFactory extends FileStoredDataFactory<Sample> {
    private boolean amplicon;

    public SampleFactory(boolean amplicon) {
        super();
        this.amplicon = amplicon;
    }

    @Override
    public Sample createObjectFromLine(String line) {
        if (line == null) return null;

        String[] sampleLines = line.split("\\s+");
        if (sampleLines.length != 8) return null;

        String sampleName = sampleLines[0].intern();
        String gene = sampleLines[1].intern();
        String chr = sampleLines[2].intern();
        int start = Integer.parseInt(sampleLines[3]);
        int end = Integer.parseInt(sampleLines[4]);
        int len = Integer.parseInt(sampleLines[6]);
        double depth = Double.parseDouble(sampleLines[7]);

        String key = amplicon ? join(" ", gene, chr, Long.toString(start), Long.toString(end), Long.toString(len)).intern() : gene;
        return new Sample(key, sampleName, chr, start, end, gene, len, depth);
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
