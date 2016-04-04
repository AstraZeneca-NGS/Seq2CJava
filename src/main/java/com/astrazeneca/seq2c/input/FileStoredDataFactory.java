package com.astrazeneca.seq2c.input;

public abstract class FileStoredDataFactory<T extends FileStoredData> {

    public abstract T createObjectFromLine(String line);
}
