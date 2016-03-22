package com.astrazeneca.seq2c;

import java.util.Map;

public class Normalization {

    String gene;
    Map<String, Double> samples;

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Normalization that = (Normalization) o;

        return gene != null ? gene.equals(that.gene) : that.gene == null;

    }

    @Override
    public int hashCode() {
        return gene != null ? gene.hashCode() : 0;
    }
}

