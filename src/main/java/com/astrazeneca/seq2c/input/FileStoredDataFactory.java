package com.astrazeneca.seq2c.input;

public abstract class FileStoredDataFactory {

    public static FileStoredDataFactory getFactory(String type, boolean amplicon) {
        switch (type) {
            case "sample":
                return new SampleFactory(amplicon);
            case "statistics":
                return new StatisticsFactory();
            default:
                throw new IllegalArgumentException("Unsupported class " + type);
        }
    }

    public abstract FileStoredData createObjectFromLine(String line);
}
