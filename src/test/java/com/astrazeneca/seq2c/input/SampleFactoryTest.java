package com.astrazeneca.seq2c.input;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class SampleFactoryTest {

    @DataProvider(name = "createObjectFromLineDataProvider")
    public Object[][] createObjectFromLineDataProvider() {
        SampleFactory sampleFactoryAmplicon = new SampleFactory(true);
        SampleFactory sampleFactoryNonAmplicon = new SampleFactory(false);
        String sample1 = "78004-normal chr21:10942921-10943020:TPTE chr21 10942921 10943020 Whole-Gene 100 201.26";
        String sample2 = "78004-normal chr21:14987444-14987558:POTED chr21 14987444 14987558 Amplicon 115 11.97";
        return new Object[][] {
                {
                    sampleFactoryNonAmplicon,
                    sample1,
                    new Sample(
                        "chr21:10942921-10943020:TPTE",
                        "78004-normal",
                        "chr21",
                        10942921,
                        10943020,
                        "chr21:10942921-10943020:TPTE",
                        100,
                        201.26
                    )
                },
                {
                    sampleFactoryNonAmplicon,
                    sample2,
                    new Sample(
                        "chr21:14987444-14987558:POTED",
                        "78004-normal",
                        "chr21",
                        14987444,
                        14987558,
                        "chr21:14987444-14987558:POTED",
                        115,
                        11.97
                    )
                },
                {
                        sampleFactoryAmplicon,
                        sample1,
                        new Sample(
                                "chr21:10942921-10943020:TPTE chr21 10942921 10943020 100",
                                "78004-normal",
                                "chr21",
                                10942921,
                                10943020,
                                "chr21:10942921-10943020:TPTE",
                                100,
                                201.26
                        )
                },
                {
                        sampleFactoryAmplicon,
                        sample2,
                        new Sample(
                                "chr21:14987444-14987558:POTED chr21 14987444 14987558 115",
                                "78004-normal",
                                "chr21",
                                14987444,
                                14987558,
                                "chr21:14987444-14987558:POTED",
                                115,
                                11.97
                        )
                }
        };
    }

    @Test(dataProvider = "createObjectFromLineDataProvider")
    public void testNonAmpliconCreateObjectFromLine(SampleFactory sampleFactory, String line, Sample expected) throws Exception {
        Sample sampleFromLine = sampleFactory.createObjectFromLine(line);
        Assert.assertEquals(sampleFromLine, expected);

        Assert.assertEquals(sampleFromLine.getName(), expected.getName());
        Assert.assertEquals(sampleFromLine.getSample(), expected.getSample());
        Assert.assertEquals(sampleFromLine.getGene(), expected.getGene());
        Assert.assertEquals(sampleFromLine.getChr(), expected.getChr());
        Assert.assertEquals(sampleFromLine.getStart(), expected.getStart());
        Assert.assertEquals(sampleFromLine.getEnd(), expected.getEnd());
        Assert.assertEquals(sampleFromLine.getCov(), expected.getCov());
        Assert.assertEquals(sampleFromLine.getLen(), expected.getLen());
        Assert.assertEquals(sampleFromLine.getNorm1b(), expected.getNorm1b());
        Assert.assertEquals(sampleFromLine.getNorm3(), expected.getNorm3());

    }

}