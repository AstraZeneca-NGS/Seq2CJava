package com.astrazeneca.seq2c.input;

/**
 * Interface of entity which can be adjusted to other FileStoredData entities.
 */
public interface FileStoredData {
    void addNextObject(FileStoredData next);
}
