package com.astrazeneca.seq2c.input;

import htsjdk.samtools.util.CloseableIterator;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;


/**
 * Represent a file abstraction and has iterator interface.
 * It will skip a string from file if it meets specified conditions.
 * */
public class FileDataIterator<T extends FileStoredData> implements CloseableIterator<T> {

    private BufferedReader bufferReader;
    private T next;
    private T current;
    private FileStoredDataFactory<T> factory;

    public FileDataIterator(String fileName, FileStoredDataFactory<T> factory) {
        try {
            bufferReader = new BufferedReader(new FileReader(fileName));
            this.factory = factory;
            current = loadNext();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    //for tests
    FileDataIterator(BufferedReader bufferedReader, FileStoredDataFactory<T> factory) {
        bufferReader = bufferedReader;
        this.factory = factory;
        current = loadNext();
    }

    @Override
    public boolean hasNext() {
        return current != null;
    }

    @Override
    public T next() {
        T tmp = current;
        current = loadNext();
        return tmp;
    }

    private T loadNext() {
        T data;
        if (next == null) {
            data = readNextObject();
        } else {
            data = next;
        }
        next = readNextObject();
        while (next != null && data.equals(next)) {
            data.addNextObject(next);
            next = readNextObject();
        }
        return data;
    }

    private T readNextObject() {
        String line = getNextLine();
        if (line == null) return null;
        next = factory.createObjectFromLine(line);
        return next;
    }

    //perl version: 48 str
    private String getNextLine() {
        String line = null;
        try {
            while (true) {
                line = bufferReader.readLine();
                if (line == null || line.isEmpty()) {
                    return null;
                }
                if (line.contains("Sample") || line.contains("Whole") || line.contains("Control") ||
                        line.contains("Undertermined")) {
                    continue;
                }
                return line;
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return line;
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("Cannot remove data from a file");
    }

    public void close() {
        try {
            bufferReader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
