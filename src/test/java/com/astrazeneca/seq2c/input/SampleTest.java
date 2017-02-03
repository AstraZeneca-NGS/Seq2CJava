package com.astrazeneca.seq2c.input;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public class SampleTest {

    @DataProvider(name = "getTitleDataProvider")
    public Object[][] getTitleDataProvider() {
        String sample1 = "chr21:10942921-10943020:TPTE\tchr21\t10942921\t10943020\t100\t";
        String sample2 = "chr21:14987444-14987558:POTED\tchr21\t14987444\t14987558\t115\t";
        return new Object[][] {
                {
                        new Sample(
                                "chr21:10942921-10943020:TPTE",
                                "78004-normal",
                                "chr21",
                                10942921,
                                10943020,
                                "chr21:10942921-10943020:TPTE",
                                100,
                                201.26
                        ),
                        sample1
                },
                {
                        new Sample(
                                "chr21:14987444-14987558:POTED",
                                "78004-normal",
                                "chr21",
                                14987444,
                                14987558,
                                "chr21:14987444-14987558:POTED",
                                115,
                                11.97
                        ),
                        sample2
                },
                {
                        new Sample(
                                "chr21:10942921-10943020:TPTE chr21 10942921 10943020 100",
                                "78004-normal",
                                "chr21",
                                10942921,
                                10943020,
                                "chr21:10942921-10943020:TPTE",
                                100,
                                201.26
                        ),
                        sample1
                },
                {
                        new Sample(
                                "chr21:14987444-14987558:POTED chr21 14987444 14987558 115",
                                "78004-normal",
                                "chr21",
                                14987444,
                                14987558,
                                "chr21:14987444-14987558:POTED",
                                115,
                                11.97
                        ),
                        sample2
                }
        };
    }

    @DataProvider(name = "getResultStringDataProvider")
    public Object[][] getResultStringProvider() {
        String sample1 = "78004-normal\tchr21:10942921-10943020:TPTE\tchr21\t10942921\t10943020\t100\t\t0.00\t0.00\t";
        String sample2 = "78004-normal\tchr21:14987444-14987558:POTED\tchr21\t14987444\t14987558\t115\t\t0.00\t0.00\t";
        return new Object[][] {
                {
                        new Sample(
                                "chr21:10942921-10943020:TPTE",
                                "78004-normal",
                                "chr21",
                                10942921,
                                10943020,
                                "chr21:10942921-10943020:TPTE",
                                100,
                                201.26
                        ),
                        sample1
                },
                {
                        new Sample(
                                "chr21:14987444-14987558:POTED",
                                "78004-normal",
                                "chr21",
                                14987444,
                                14987558,
                                "chr21:14987444-14987558:POTED",
                                115,
                                11.97
                        ),
                        sample2
                },
                {
                        new Sample(
                                "chr21:10942921-10943020:TPTE chr21 10942921 10943020 100",
                                "78004-normal",
                                "chr21",
                                10942921,
                                10943020,
                                "chr21:10942921-10943020:TPTE",
                                100,
                                201.26
                        ),
                        sample1
                },
                {
                        new Sample(
                                "chr21:14987444-14987558:POTED chr21 14987444 14987558 115",
                                "78004-normal",
                                "chr21",
                                14987444,
                                14987558,
                                "chr21:14987444-14987558:POTED",
                                115,
                                11.97
                        ),
                        sample2
                }
        };
    }

    @Test(dataProvider = "getResultStringDataProvider")
    public void testGetResultString(Sample sample, String expectedResultString) throws Exception {
        Assert.assertEquals(sample.getResultString(), expectedResultString);
    }

    @Test(dataProvider = "getTitleDataProvider")
    public void testGetTitle(Sample sample, String expectedTitle) throws Exception {
        Assert.assertEquals(sample.getTitle(), expectedTitle);
    }

}