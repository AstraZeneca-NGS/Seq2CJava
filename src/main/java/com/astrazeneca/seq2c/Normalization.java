package com.astrazeneca.seq2c;

import java.util.List;
import java.util.Map;

public class Normalization {

    String sample;
    Map<String, List<Sample>> genes;

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Normalization that = (Normalization) o;

        return sample != null ? sample.equals(that.sample) : that.sample == null;

    }

    @Override
    public int hashCode() {
        return sample != null ? sample.hashCode() : 0;
    }
}

