package com.astrazeneca.seq2c;

import com.astrazeneca.seq2c.input.FileStoredData;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class SampleStatistics implements FileStoredData {

    private String sample;
    private Map<String, List<Sample>> genes;

    public SampleStatistics(String sample) {
        this.sample = sample;
    }

    public String getSample() {
        return sample;
    }

    public void setGenes(Map<String, List<Sample>> genes) {
        this.genes = genes;
    }

    public Map<String, List<Sample>> getGenes() {
        return genes;
    }

    @Override

    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        SampleStatistics that = (SampleStatistics) o;

        return sample != null ? sample.equals(that.sample) : that.sample == null;

    }

    @Override
    public int hashCode() {
        return sample != null ? sample.hashCode() : 0;
    }

    @Override
    public void addNextObject(FileStoredData next) {
        SampleStatistics statistics = (SampleStatistics)next;
        for(Map.Entry<String, List<Sample>> entry : statistics.getGenes().entrySet()) {
            String gene = entry.getKey();
            List<Sample> list;
            if (genes.containsKey(gene)) {
                list = genes.get(gene);
            } else {
                list = new ArrayList<>();
                genes.put(gene, list);
            }
            list.addAll(entry.getValue());
        }
    }
}

