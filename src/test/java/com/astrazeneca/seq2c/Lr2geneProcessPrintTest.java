package com.astrazeneca.seq2c;

import com.astrazeneca.seq2c.input.SampleStatistics;
import com.astrazeneca.seq2c.input.StatisticsFactory;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import org.mockito.Mockito;
import org.powermock.api.mockito.PowerMockito;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.FileDescriptor;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.io.StringReader;
import java.util.ArrayDeque;
import java.util.Deque;
import java.util.Map;

import static org.testng.Assert.assertEquals;

public class Lr2geneProcessPrintTest {

    private ByteArrayOutputStream outBytes;

    @BeforeClass
    public void setUp() {
        outBytes = new ByteArrayOutputStream();
        System.setOut(new PrintStream(outBytes));
    }

    @AfterClass
    public void tearDown(){
        //reset System.out
        System.setOut(new PrintStream(new FileOutputStream(FileDescriptor.out)));
    }


    @DataProvider
    public Object[][] lr2geneDataProvider() {
        return new Object[][] {
                {
                    "HG0004\tnot_a_gene_4\tmockchr\t0\t0\t0\t-5.335\t3.2e+02\tBP\tDel\t3\t4\t-5.38\t3.2\t1,2,3(62132652-62363432)\t1\t0.2",
                    "HG0003\tnot_a_gene_3\tmockchr\t0\t0\t0\t0.25\t\t\t\t\t1",
                    "HG0002\tnot_a_gene_2\tmockchr\t0\t0\t0\t2.25\t0.00\tWhole\tAmp\t1\t1\t2.25\t2.3\tALL\t0\t0.0",
                    "HG0001\tnot_a_gene_1\tmockchr\t0\t0\t0\t-5.49\t0.00\tWhole\tDel\t1\t1\t-5.49\t-5.5\tALL\t0\t0.0",
                    "HG0008\tnot_a_gene_5\tmockchr\t0\t0\t0\t5.235\t0.00\tWhole\tAmp\t2\t2\t5.235\t5.2\tALL\t0\t0.0",
                    null
                }
        };
    }

    @Test(dataProvider = "lr2geneDataProvider")
    public void lr2geneProcessAcceptanceTest(String[] expected) throws IOException {

        Lr2gene lr2gene = createLr2geneWithMockLocusMap();
        lr2gene.process();

        try (BufferedReader outData = new BufferedReader(new StringReader(new String(outBytes.toByteArray())))) {
            /*skip header*/ outData.readLine();
            for (String expectedLine : expected) {
                assertEquals(outData.readLine(), expectedLine);
            }
        }
    }

    private Lr2gene createLr2geneWithMockLocusMap() {
        Map locusMapMock = Mockito.mock(Map.class);
        PowerMockito
                .when(locusMapMock.get(Mockito.any()))
                .thenReturn(new Locus("mockname", 0, 0, "mockchr"));

        return new Lr2gene(locusMapMock, "mockfile") {
            @Override
            CloseableIterator<SampleStatistics> getFileDataIterator() {
                return new StatisticIteratorMock();
            }
        };
    }

    static class StatisticIteratorMock extends AbstractIterator<SampleStatistics> implements CloseableIterator<SampleStatistics> {
        Deque<SampleStatistics> statistics = new ArrayDeque<>();
        {
            statistics.push(createStats("1", "1", 71963329, 71963668, -5.49));
            statistics.push(createStats("2", "2", 71967724, 71968138, 2.25));
            statistics.push(createStats("3", "3", 71986683, 71987097, 0.25));

            final SampleStatistics statistic = createStats("4", "4", 62132652, 62133086, -5.47);
            statistic.addNextObject(createStats("5", "4", 62227686, 62228277, -5.46));
            statistic.addNextObject(createStats("6", "4", 62362822, 62363432, -5.21));
            statistic.addNextObject(createStats("7", "4", 62373383, 62374550, -2.23));
            statistics.push(statistic);

            final SampleStatistics statisticSeveral = (createStats("8", "5", 71963802, 71964197, 5.24));
            statisticSeveral.addNextObject(createStats("9", "5", 71967163, 71967526, 5.23));
            statistics.push(statisticSeveral);
        }
        @Override protected SampleStatistics advance() {
            if (statistics.isEmpty())
                return null;
            return statistics.pop();
        }

        @Override public void close() {}
    }

    private static SampleStatistics createStats(String sample, String gene, int start, int end, double norm3) {
        sample = "HG000" + sample;
        gene = "not_a_gene_" + gene;
        String line = sample + '\t' + gene + '\t' + "16" + '\t' + start + '\t' + end + '\t' + 149 + '\t' + norm3 + '\t' + 0.25;
        return new StatisticsFactory().createObjectFromLine(line);
    }


}
