package com.astrazeneca.seq2c.input;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Iterator;

public class FileDataIterator<T extends FileStoredData> implements Iterator<FileStoredData> {

    private BufferedReader bufferReader;
    private T next;
    private T current;
    private FileStoredDataFactory<T> factory;

    public FileDataIterator(String fileName, FileStoredDataFactory factory) {
        try {
            bufferReader = new BufferedReader(new FileReader(fileName));
            this.factory = factory;
            current = loadNext();
        } catch (IOException e) {
            e.printStackTrace();
        }
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
        while (true) {
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
    }

    private T readNextObject() {
        String line = getNextLine();
        if (line == null) return null;
        next = factory.createObjectFromLine(line);
        return next;
    }

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
