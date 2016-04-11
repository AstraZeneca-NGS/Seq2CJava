package com.astrazeneca.seq2c;

import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class TTestStatisticsTest {
    @Test
    public void getT() throws Exception {
        double[] x = {-0.19, 0.47, 1.03, 1.54, 1.58, 1.64, 1.72, 1.84};
        double[] y = {-5.12, -1.76, -1.4, -1.29, -1.13, -0.91, -0.78, -0.19, -0.13};
        double result = TTestStatistics.getT(x, y);
        assertEquals(0.0004177, result, 0.000001);
    }

    @Test
    public void getALDH9A() throws Exception {
        double[] x = {0.93, 0.45, 1.32};
        double[] y = {0.01, -0.12, 0.03, -0.09, -0.31, 0.10, -0.19};

        double result = TTestStatistics.getT(x, y);
        assertEquals(0.062319, result, 0.000001);
    }
}