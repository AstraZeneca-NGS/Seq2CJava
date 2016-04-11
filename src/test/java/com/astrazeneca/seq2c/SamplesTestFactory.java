package com.astrazeneca.seq2c;

import com.astrazeneca.seq2c.input.Sample;

public class SamplesTestFactory {
    private String name;
    private String gene;
    private String chr;

    public SamplesTestFactory(String sample, String gene, String chr) {
        this.name = sample;
        this.gene = gene;
        this.chr = chr;
    }

    public Sample getSample(int start, int end, double norm3) {
        Sample sample = new Sample(name, name, chr, start, end, gene, 0, 0);
        sample.setNorm3(norm3);
        return sample;
    }
}
