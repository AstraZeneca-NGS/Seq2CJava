package com.astrazeneca.seq2c.input;

import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import static org.testng.Assert.*;

public class SampleStatisticsTest {

    @DataProvider(name = "lines")
    Object[][] lineData() {
        return new Object[][] {
                {"seq2c-test", "TPTE:11057796-11058160:BAGE", "chr21", 10906903, 10907042, 140, 0.0, 0.1},
                {"seq2c-test", "chr21:19641434-19642441:PRSS7", "chr21", 19641434, 19642441, 1008, 0.2, 0.0},
                {"seq2c-test", "chr21:43334707-43334830:C2CD2", "chr21", 43334707, 43334830, 124, 0.3, 1.2},
        };
    }

    @Test(dataProvider = "lines")
    public void factoryShouldCreateProperSampleStatistics(
            String sample, String gene, String chr, int start, int end, int len, double norm3, double norm3s) {
        String line = createLine(sample, gene, chr, start, end, len, norm3, norm3s);

        StatisticsFactory factory = new StatisticsFactory();
        SampleStatistics statistics = factory.createObjectFromLine(line);

        assertEquals(statistics.getSample(), sample);
        assertEquals(statistics.getGenes().size(), 1);
        assertTrue(statistics.getGenes().containsKey(gene));
        assertEquals(statistics.getGenes().get(gene).size(), 1);

        final Sample smpl = statistics.getGenes().get(gene).get(0);
        assertEquals(smpl.getSample(), sample);
        assertEquals(smpl.getChr(), chr);
        assertEquals(smpl.getGene(), gene);
        assertEquals(smpl.getStart(), start);
        assertEquals(smpl.getEnd(), end);
        assertEquals(smpl.getLen(), len);
        assertEquals(smpl.getNorm3(), norm3);
        assertEquals(smpl.getNorm3s(), norm3s);
    }

    @Test
    public void addNewStatisticsOnTheSameGeneTest() {
        StatisticsFactory factory = new StatisticsFactory();

        final String gene = "chr21:43640008-43640156:ABCG1";
        SampleStatistics stat1 =
                factory.createObjectFromLine(createLine(
                        "seq2c-test",
                        gene,
                        "chr21",
                        43640008,
                        43640156,
                        149,
                        0.00,
                        0.00));

        SampleStatistics stat2 =
                factory.createObjectFromLine(createLine(
                        "seq2c-test",
                        gene,
                        "chr21",
                        43645781,
                        43646024,
                        244,
                        0.00,
                        0.00));

        stat1.addNextObject(stat2);

        assertEquals(stat1.getGenes().size(), 1);
        assertEquals(stat1.getGenes().get(gene).size(), 2);
    }

    @Test
    public void addNewStatisticsOnDifferentGenesTest() {
        StatisticsFactory factory = new StatisticsFactory();

        final String gene1 = "chr21:43562206-43563105:UMODL1";
        SampleStatistics stat1 =
                factory.createObjectFromLine(createLine(
                        "seq2c-test",
                        gene1,
                        "chr21",
                        43640008,
                        43640156,
                        149,
                        0.00,
                        0.00));

        final String gene2 = "chr21:43640008-43640156:ABCG1";
        SampleStatistics stat2 =
                factory.createObjectFromLine(createLine(
                        "seq2c-test",
                        gene2,
                        "chr21",
                        43645781,
                        43646024,
                        244,
                        0.00,
                        0.00));

        stat1.addNextObject(stat2);

        assertEquals(stat1.getGenes().size(), 2);
        assertEquals(stat1.getGenes().get(gene1).size(), 1);
        assertEquals(stat1.getGenes().get(gene2).size(), 1);

    }

    @Test
    public void equalsShouldWorkOnlyWithSampleNames() {
        StatisticsFactory factory = new StatisticsFactory();
        final String sample1 = "seq2c-test";
        SampleStatistics stat1 = factory.createObjectFromLine(createLine(
                        sample1,
                        "chr21:43562206-43563105:UMODL1","chr21",43640008,43640156,
                        149, 0.00, 0.00));

        final String sample2 = "seq2c-test2";
        SampleStatistics stat2 = factory.createObjectFromLine(createLine(
                        sample2,
                        "chr21:43640008-43640156:ABCG1","chr21",43645781, 43646024,
                        244, 0.00, 0.00));

        assertTrue(stat1.equals(stat1));
        assertTrue(stat2.equals(stat2));
        assertFalse(stat2.equals(stat1));
        assertFalse(stat1.equals(stat2));
    }

    private static String createLine(
            String sample, String gene, String chr, int start, int end, int len, double norm3, double norm3s) {
        return sample + '\t' + gene + '\t' + chr + '\t' + start + '\t' + end + '\t' + len + '\t' + norm3 + '\t' + norm3s;
    }
}