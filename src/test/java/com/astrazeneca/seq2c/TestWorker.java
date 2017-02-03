package com.astrazeneca.seq2c;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.concurrent.LinkedBlockingQueue;

import static org.mockito.Matchers.anyString;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;
import static org.testng.Assert.assertEquals;

public class TestWorker {

    @DataProvider(name = "testData")
    public Object[][] testData() {
        SAMRecord supplementaryAlignment = mock(SAMRecord.class);
        when(supplementaryAlignment.getStringAttribute(anyString())).thenReturn("test");
        when(supplementaryAlignment.getFlags()).thenReturn(0x800);


        SAMRecord record_1 = TestSAMReader.createRecord(20,1, new ArrayList<CigarElement>(){{
            add(new CigarElement(2, CigarOperator.M));
            add(new CigarElement(3, CigarOperator.D));
        }});
        SAMRecord record_2 = TestSAMReader.createRecord(15, 1, new ArrayList<CigarElement>(){{
            add(new CigarElement(2, CigarOperator.S));
            add(new CigarElement(2, CigarOperator.M));
            add(new CigarElement(1, CigarOperator.D));
        }});
        SAMRecord record_3 = TestSAMReader.createRecord(15, 255, new ArrayList<CigarElement>(){{
            add(new CigarElement(2, CigarOperator.M));
            add(new CigarElement(1, CigarOperator.D));
            add(new CigarElement(2, CigarOperator.S));
        }});

        return new Object[][]{
                {supplementaryAlignment, 10, 20, 0, false},
                {record_1, 10, 20, 1, false},
                {record_1, 10, 20, 0, true},
                {record_2, 10, 20, 1, true},
                {record_3, 10, 20, 1, true}
        };
    }

    @Test(dataProvider = "testData")
    public void testProcessRecord(Object record, int start, int end, int exonCoverage, boolean PCRamplbc) {
        Seq2cov.Worker testWorker = getTestWorker(PCRamplbc);
        testWorker.processRecord(start, end, (SAMRecord) record);
        assertEquals(testWorker.exoncov, exonCoverage);
    }

    private Seq2cov.Worker getTestWorker(boolean PCRamplbc) {
        Seq2cov.Worker testWorker = new Seq2cov.Worker(new TestSAMReader());
        testWorker.ctx = new Seq2cov.GeneCtx(new Seq2cov.Seq2covGene("","", new int[]{}),
                Dispatcher.getService(1), new LinkedBlockingQueue<Seq2cov.Worker>(),
                new ArrayList<Gene>(), "", PCRamplbc, true);
        return testWorker;
    }
}
