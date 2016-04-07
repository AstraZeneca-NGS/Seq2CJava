package com.astrazeneca.seq2c;

import com.astrazeneca.seq2c.input.Sample;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertEquals;

/**
 * Created by Mariia_Zueva on 4/5/2016.
 */
public class Lr2geneTest {
    @org.junit.Test
    public void checkBPTestCOMMD1() throws Exception {
        String sampleName = "HG00097.mapped.ILLUMINA.bwa.GBR.exome.20130415";
        String gene = "COMMD1";
        List<Sample> samples = new ArrayList<>();
        Sample sample = new Sample("sampleName", "sampleName", "2", 62132652, 62133086, gene, 435, 19.33);
        sample.setNorm3(-0.03);
        samples.add(sample);

        Sample sample1 = new Sample("sampleName", "sampleName", "2", 62227686, 62228277, gene, 435, 19.33);
        sample1.setNorm3(-0.45);
        samples.add(sample1);

        Sample sample2 = new Sample("sampleName", "sampleName", "2", 62362822, 62363432, gene, 435, 19.33);
        sample2.setNorm3(-0.77);
        samples.add(sample2);

        Sample sample3 = new Sample("sampleName", "sampleName", "2", 62373383, 62374550, gene, 435, 19.33);
        sample3.setNorm3(-4.23);
        samples.add(sample3);

        Sig sig = new Lr2gene().checkBP(samples);
        assertEquals(12, sig.getSig(), 0.001);
        assertEquals("BP", sig.getBp());
        assertEquals("Del", sig.getCn());
        assertEquals(1.0, sig.getBpi(), 0.001);
        assertEquals(4, sig.getTotal());
        assertEquals(-4.23, sig.getSiglr(), 0.001);
        assertEquals(3.813, sig.getSigdiff(), 0.001);
        assertEquals("4", sig.getSigseg());

        System.out.println(sig.getName());
    }

    @org.junit.Test
    public void checkBPTestPKN2() throws Exception {

        String sampleName = "HG00097.mapped.ILLUMINA.bwa.GBR.exome.20130415";
        String gene = "PKN2";
        String chr = "1";
        List<Sample> samples = new ArrayList<>();
        Sample sample = new Sample(sampleName, sampleName, chr, 89149744, 89150413, gene, 670, 19.33);
        sample.setNorm3(-2.15);
        samples.add(sample);

        sample = new Sample(sampleName, sampleName, chr, 89206525, 89207087, gene, 563, 19.33);
        sample.setNorm3(0.43);
        samples.add(sample);

        sample = new Sample(sampleName, sampleName, chr, 89225761, 89226129, gene, 369, 19.33);
        sample.setNorm3(0.40);
        samples.add(sample);

        sample = new Sample(sampleName, sampleName, chr, 89235912, 89236276, gene, 365, 19.33);
        sample.setNorm3(-0.07);
        samples.add(sample);

        sample = new Sample(sampleName, sampleName, chr, 89237023, 89237636, gene, 614, 19.33);
        sample.setNorm3(-0.35);
        samples.add(sample);

        sample = new Sample(sampleName, sampleName, chr, 89250208, 89250579, gene, 372, 19.33);
        sample.setNorm3(0.36);
        samples.add(sample);

        sample = new Sample(sampleName, sampleName, chr, 89251641, 89252008, gene, 368, 19.33);
        sample.setNorm3(0.15);
        samples.add(sample);

        sample = new Sample(sampleName, sampleName, chr, 89270011, 89272037, gene, 368, 19.33);
        sample.setNorm3(-0.42);
        samples.add(sample);

        sample = new Sample(sampleName, sampleName, chr, 89272852, 89273566, gene, 715, 19.33);
        sample.setNorm3(0.01);
        samples.add(sample);

        sample = new Sample(sampleName, sampleName, chr, 89279185, 89279534, gene, 350, 19.33);
        sample.setNorm3(1.05);
        samples.add(sample);

        sample = new Sample(sampleName, sampleName, chr, 89287463, 89287813, gene, 351, 19.33);
        sample.setNorm3(0.08);
        samples.add(sample);

        sample = new Sample(sampleName, sampleName, chr, 89289842, 89290216, gene, 375, 19.33);
        sample.setNorm3(-0.09);
        samples.add(sample);

        sample = new Sample(sampleName, sampleName, chr, 89294053, 89294415, gene, 363, 19.33);
        sample.setNorm3(0.91);
        samples.add(sample);

        sample = new Sample(sampleName, sampleName, chr, 89294531, 89294900, gene, 370, 19.33);
        sample.setNorm3(-1.24);
        samples.add(sample);

        sample = new Sample(sampleName, sampleName, chr, 89298286, 89301988, gene, 3703, 19.33);
        sample.setNorm3(1.36);
        samples.add(sample);

        Sig sig = new Lr2gene().checkBP(samples);
        System.out.println(sig.getName());
    }

}