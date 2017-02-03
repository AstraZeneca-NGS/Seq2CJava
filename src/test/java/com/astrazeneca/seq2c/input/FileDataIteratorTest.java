package com.astrazeneca.seq2c.input;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.BufferedReader;
import java.io.CharArrayReader;
import java.io.IOException;

public class FileDataIteratorTest {
    @DataProvider (name = "fileDataIteratorDataProvider")
    public Object[][] fileDataIteratorDataProvider() {
        return new Object[][] {
                {
                        new FileDataIterator<>(new TestBufferedReader(new String[]{
                                "78004-normal chr17:14987449-14987555:POTED chr17 14987448 14987558 Amplicon 115 11.97",
                                "78004-normal chr18:14987444-14987558:POTED chr18 14987444 14987558 Amplicon 115 11.97",
                                "78004-normal chr19:14987445-14987558:POTED chr19 14987445 14987558 Amplicon 112 11.97",
                                "78004-normal chr20:14987442-14987558:POTED chr20 14987446 14987558 Amplicon 115 11.97",
                                "78004-normal chr21:14987447-14987558:POTED chr21 14987447 14987558 Amplicon 114 11.97"
                        }), new SampleFactory(false)),
                        5
                },
                {
                    // 4 because we have one Whole-Gene sample, and it will be filtered
                        new FileDataIterator<>(new TestBufferedReader(new String[]{
                                "78004-normal chr18:14987444-14987558:POTED chr18 14987444 14987558 Amplicon 115 11.97",
                                "78004-normal chr19:14987445-14987558:POTED chr19 14987445 14987558 Amplicon 112 11.97",
                                "78004-normal chr17:14987449-14987555:POTED chr17 14987448 14987558 Whole-Gene 115 11.97",
                                "78004-normal chr20:14987442-14987558:POTED chr20 14987446 14987558 Amplicon 115 11.97",
                                "78004-normal chr21:14987447-14987558:POTED chr21 14987447 14987558 Amplicon 114 11.97"
                        }), new SampleFactory(false)),
                        4
                },
        };
    }

    @Test(dataProvider = "fileDataIteratorDataProvider")
    public void hasNextShouldReturnFalseWhenAllDataWasRead(FileDataIterator<Sample> dataIterator, int sampleCount) throws Exception {
        Assert.assertTrue(dataIterator.hasNext());
        for (int i = 0; i < sampleCount; i++) {
            Assert.assertNotNull(dataIterator.next());
        }
        Assert.assertFalse(dataIterator.hasNext());
    }

    @Test(expectedExceptions = UnsupportedOperationException.class)
    public void removeShouldThrowUnsupportedOperationException() throws Exception {
        FileDataIterator<Sample> dataIterator = getSampleFileDataIterator();
        dataIterator.remove();
    }

    private FileDataIterator<Sample> getSampleFileDataIterator() {
        return new FileDataIterator<>(new BufferedReader(new CharArrayReader(new char[0])), new SampleFactory(false));
    }

    private class TestBufferedReader extends BufferedReader {

        int count = 0;
        String[] data;

        TestBufferedReader(String[] data) {
            super(new CharArrayReader(new char[0]));
            this.data = data;
        }


        @Override
        public String readLine() throws IOException {
            String tmp;
            if (count >= data.length) {
                tmp = "";
            } else {
                tmp = data[count];
            }
            count++;
            return tmp;
        }
    }
}