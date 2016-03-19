package com.astrazeneca.seq2c;

import org.junit.Test;


public class SampleIteratorTest {

    @Test
    public void testSampleIterator() {
        SampleIterator iterator = new SampleIterator("cov_java.txt", false);

        Sample sample = iterator.nextSample();
        while (sample != null) {
            System.out.println(sample);
            sample = iterator.nextSample();
        }
        iterator.close();

    }
}