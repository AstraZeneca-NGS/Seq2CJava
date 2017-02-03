package com.astrazeneca.seq2c;

import com.astrazeneca.seq2c.input.Sample;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import org.apache.commons.math3.util.Precision;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.*;
import java.util.*;

import static org.testng.Assert.assertEquals;

public class Cov2lrTest {

    private final static Sample badSample = new Sample("gene_name1 chr21 9907194 9908432 1239",
            "sample_name1", "chr21", 9907194, 9908432,
            "gene_name1", 1239, 2.0);
    private final static Sample goodSample = new Sample("gene_name2 chr21 9909047 9909277 231",
            "sample_name1", "chr21", 9909047, 9909277,
            "gene_name2", 231, 100.0);
    private final static Sample goodSample_2 = new Sample("gene_name2 chr21 9909047 9909277 231",
            "sample_name2", "chr21", 9909047, 9909277,
            "gene_name2", 231, 80.0);
    private final static Sample goodSample_3 = new Sample("gene_name2 chr21 9909045 9909279 231",
            "sample_name2", "chr21", 9909045, 9909279,
            "gene_name2", 235, 80.0);
    private final static Map<String, Long> mappingReads = new HashMap<String, Long>() {{
        put("sample_name1", 3L);
        put("sample_name2", 1L);
    }};

    @Test
    public void testReadCoverageAndGetMedDepth() {
        Cov2lr cov2lr = new Cov2lrForTestReadCoverage();
        double actualResult = cov2lr.readCoverageAndGetMedDepth(new HashSet<String>());
        assertEquals(Precision.round(actualResult, 3), 113.335);
    }

    @DataProvider(name = "dataForFilteredDataStreamTest")
    public Object[][] getDataForFilteredDataStreamTest() {
        return new Object[][]{
                {
                    new HashSet<String>() {{
                        add("sample_name1");
                    }}, new HashMap<String, Double>() {{
                        put("sample_name1", 2.0);
                        put("sample_name1", 2.67);
                    }}
                },
        };
    }

    @Test(dataProvider = "dataForFilteredDataStreamTest")
    public void testFilterDataStream(HashSet<String> sampleSet, HashMap<String, Double> map) {
        Cov2lr cov2lr = getCov2lr();
        double actualPercentile = cov2lr.filterDataStream(sampleSet, map, new ArrayList<Double>());
        assertEquals(actualPercentile, 2.67);
    }

    @DataProvider(name = "geneStatisticsMapData")
    public Object[][] getGeneStatisticsMapData() {
        return new Object[][]{
                {
                        new HashMap<String, Map<String, Double>>() {{
                            put("gene_name2 chr21 9909047 9909277 231", new HashMap<String, Double>() {{
                                put("sample_name1", 2.7);
                            }});
                        }}
                }
        };
    }

    @Test(dataProvider = "geneStatisticsMapData")
    public void testShouldFillIntoGeneStatisticsMap(HashMap<String, Map<String, Double>> expectedMap) {
        Cov2lr cov2lr = getCov2lr();
        cov2lr.putInMap(goodSample, 2.7);
        assertEquals(cov2lr.getGeneStatistics(), expectedMap);
    }

    @DataProvider(name = "gooddepth")
    public Object[][] getSplitQualitySamplesData() {
        return new Object[][]{
                {
                    2.67, new LinkedList<Double>() {{
                        add(2.0);
                        add(100.0);
                    }}
                },
                {
                    50.0, new LinkedList<Double>() {{
                        add(100.0);
                    }}
                }
        };
    }

    @Test(dataProvider = "gooddepth")
    public void testSplitQualitySamples(double medDepth, List<Double> expected) {
        Cov2lr cov2lr = getCov2lr();
        cov2lr.putInMap(badSample, 2.0);
        cov2lr.putInMap(goodSample, 100.0);
        List<Double> gooddepth = cov2lr.splitQualitySamples(
                new HashSet<String>() {{
                    add("sample_name1");
                    add("sample_name2");
                }},
                medDepth,
                new HashSet<String>());
        assertEquals(expected, gooddepth);
    }

    @DataProvider(name = "factor2")
    public Object[][] getFactor2Data() {
        return new Object[][]{
                {new HashSet<String>() {{
                    add("gene_name1 chr21 9907194 9908432 1239");}},
                        new HashMap<String, Double>() {{
                            put("gene_name2 chr21 9909047 9909277 231", 0.5);
                        }}
                },
                {new HashSet<String>(),
                        new HashMap<String, Double>() {{
                            put("gene_name1 chr21 9907194 9908432 1239", 25.0);
                            put("gene_name2 chr21 9909047 9909277 231", 0.5);
                        }}
                }
        };
    }

    @Test(dataProvider = "factor2")
    public void testGetFactor2(Set<String> badGenes, Map<String, Double> expected) {
        Cov2lr cov2lr = getCov2lr();
        cov2lr.putInMap(badSample, 2.0);
        cov2lr.putInMap(goodSample, 100.0);

        Map<String, Double> factor2 = cov2lr.getFactor2(50, badGenes);
        assertEquals(factor2, expected);
    }

    @DataProvider(name = "sampleMedian")
    public Object[][] getSampleMedianData() {
        return new Object[][]{
                {new HashSet<String>() {{
                    add("gene_name1 chr21 9907194 9908432 1239");}},
                        new LinkedHashMap<String, Double>(){{
                            put("sample_name1", 100.0);}}
                },
                {new HashSet<String>(),
                        new LinkedHashMap<String, Double>(){{
                            put("sample_name1", 51.0);}}
                }
        };
    }

    @Test(dataProvider = "sampleMedian")
    public void testGetSampleMedian(Set<String> badGenes, Map<String,Double> expected) {
        Cov2lr cov2lr = getCov2lr();
        cov2lr.putInMap(badSample, 2.0);
        cov2lr.putInMap(goodSample, 100.0);

        Map<String, Double> sampleMedian = cov2lr.getSampleMedian(
                new HashSet<String>() {{
                    add("sample_name1");
                }},
                badGenes);
        assertEquals(sampleMedian, expected);
    }

    @DataProvider(name = "controlSamples")
    public Object[][] getDataForControlSamplesTest() {
        return new Object[][] {
                {null, 0.0},
                {"sample_name1", 50.1},
                {"sample_name1:sample_name2", 46.77}
        };
    }

    @Test(dataProvider = "controlSamples")
    public void testPrepareControlSamples(String controlSamples, double expectedValue) {
        Cov2lr cov2lr = getCov2lr(controlSamples);
        cov2lr.putInMap(badSample, 2.0);
        cov2lr.putInMap(goodSample, 100.0);
        cov2lr.putInMap(goodSample_2, 80.0);
        Map<String, Double> factor = new HashMap<String, Double>() {{
            put("gene_name1 chr21 9907194 9908432 1239", 25.0);
            put("gene_name2 chr21 9909047 9909277 231", 0.5);
        }};
        double controlSampleMean = cov2lr.prepareControlSamples(factor);
        assertEquals(Precision.round(controlSampleMean, 2), expectedValue);
    }

    @DataProvider(name = "setNormalization")
    public Object[][] getSetNormalization() {
        return new Object[][] {
                {new HashMap<String, Double>() {{
                    put("gene_name1 chr21 9907194 9908432 1239", 25.0);
                    put("gene_name2 chr21 9909047 9909277 231", 0.5);}},
                new LinkedHashMap<String, Double>(){{
                    put("sample_name1", 51.0);
                    put("sample_name2", 51.0);
                    }},
                        46.77,
                new HashSet<String>() {{
                    add("gene_name1 chr21 9907194 9908432 1239");}},
                        null
                },
                {new HashMap<String, Double>() {{
                    put("gene_name2 chr21 9909047 9909277 231", 0.5);}},
                new LinkedHashMap<String, Double>(){{
                    put("sample_name1", 100.0);
                    put("sample_name2", 100.0);
                    }},
                        50.1,
                new HashSet<String>(),
                        null
                },
                {new HashMap<String, Double>() {{
                    put("gene_name2 chr21 9909047 9909277 231", 0.5);}},
                new LinkedHashMap<String, Double>(){{
                    put("sample_name1", 100.0);
                    put("sample_name2", 100.0);
                    }},
                        50.1,
                new HashSet<String>(),
                    "sample_name1:sample_name2"
                }
        };
    }

    @Test(dataProvider = "setNormalization")
    public void testSetNormalization(Map<String, Double> factor2,
                                     Map<String, Double> sampleMedian,
                                     double controlSamplesMean,
                                     Set<String> badGenes,
                                     String useControlSamples) throws IOException {
        Cov2lr cov2lr = new Cov2lrForTestSetNormalization(useControlSamples);
        cov2lr.putInMap(goodSample, 2.0);
        cov2lr.putInMap(goodSample_2, 100.0);
        List<Sample> sampleResult = cov2lr.setNormalization(50, factor2, sampleMedian, controlSamplesMean, badGenes);
        assertEquals(sampleResult, new ArrayList<Sample>(){{
            add(goodSample_2);
            add(goodSample);
        }});
    }

    @DataProvider(name = "dataToAddLocus")
    public Object[][] getDataForAddLocusTest() {
        return new Object[][] {
                {goodSample, goodSample, "chr21\t9909047\t9909277\t462"},
                {goodSample, badSample, "chr21\t9909047\t9909277\t231"}, // test should do not add this sample to locus
                {goodSample, goodSample_3, "chr21\t9909045\t9909279\t466"},
        };
    }

    @Test(dataProvider = "dataToAddLocus")
    public void testAddLocus(Sample sample1, Sample sample2, String expected) {
        Cov2lr cov2lr = getCov2lr();
        cov2lr.addLocus(sample1);
        cov2lr.addLocus(sample2);
        assertEquals(cov2lr.getLocusMap().get("gene_name2").getName().toString(), expected);
    }

    private Cov2lr getCov2lr(String controlSample) {
        return new Cov2lr(
                true,
                mappingReads,
                "path_to_coverage_file",
                controlSample,
                "path_to_statistic_file");
    }

    private Cov2lr getCov2lr() {
        String controlSamples = null;
        return getCov2lr(controlSamples);
    }

    static class CoverageIteratorMock extends AbstractIterator<Sample> implements CloseableIterator<Sample> {
        Deque<Sample> samples = new ArrayDeque<>();
        {
            samples.push(badSample);
            samples.push(goodSample);
            samples.push(goodSample_2);
            samples.push(goodSample_3);
        }
        @Override protected Sample advance() {
            if (samples.isEmpty())
                return null;
            return samples.pop();
        }

        @Override public void close() {}
    }

    private static class Cov2lrForTestSetNormalization extends Cov2lr {
        public Cov2lrForTestSetNormalization(String useControlSamples) {
            super(true, Cov2lrTest.mappingReads, "path_to_coverage_file", useControlSamples, "path_to_statistic_file");
        }

        @Override
        PrintWriter getPrintWriter() {
            return new PrintWriter(new StringWriter());
        }

        @Override
        CloseableIterator<Sample> getFileDataIterator() {
            return new CoverageIteratorMock();
        }
    }

    private static class Cov2lrForTestReadCoverage extends Cov2lr {
        public Cov2lrForTestReadCoverage() {
            super(true, Cov2lrTest.mappingReads, "path_to_coverage_file", null, "path_to_statistic_file");
        }

        @Override
        CloseableIterator<Sample> getFileDataIterator() {
            return new CoverageIteratorMock();
        }
    }
}
