package com.astrazeneca.seq2c;


import htsjdk.samtools.*;
import org.powermock.api.mockito.PowerMockito;
import org.powermock.core.classloader.annotations.PrepareForTest;
import org.powermock.modules.testng.PowerMockObjectFactory;
import org.testng.IObjectFactory;
import org.testng.annotations.DataProvider;
import org.testng.annotations.ObjectFactory;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.ExecutionException;

import static com.astrazeneca.seq2c.Seq2cov.*;
import static org.mockito.Matchers.any;
import static org.mockito.Mockito.mock;
import static org.powermock.api.mockito.PowerMockito.mockStatic;
import static org.powermock.api.mockito.PowerMockito.when;
import static org.testng.AssertJUnit.assertEquals;

@PrepareForTest(SamReaderFactory.class)
public class Seq2covTest {

    Seq2cov seq2cov = new Seq2cov("", "", "");

    @ObjectFactory
    public IObjectFactory setObjectFactory() {
        return new PowerMockObjectFactory();
    }

    @DataProvider(name = "processBedTestData")
    public Object[][] testData() {
        Map<String, Seq2covGene> result = getProcessBedTestResult();
        return new Object[][]{
                {"src/test/resources/testBed.bed", result, false},
                {"src/test/resources/testAmpBed.bed", result, true}
        };
    }

    @Test(dataProvider = "processBedTestData")
    public void testProcessBed(String filePath, Map<String, Seq2covGene> expected, boolean PCRamplbc) throws IOException {
        Seq2cov seq2cov = new Seq2cov(filePath, "", "");
        Map<String, Seq2covGene> actual = seq2cov.processBed();
        assertEquals(actual, expected);
        assertEquals(seq2cov.PCRamplbc, PCRamplbc);
    }

    private Map<String, Seq2covGene> getProcessBedTestResult() {
        Map<String, Seq2covGene> result = new LinkedHashMap<>();
        Seq2covGene doubledCDC = new Seq2covGene("chr21",
                "chr21:9907194-9908432:LOC100132288", new int[]{9907194, 9908432});
        doubledCDC.addCDS(new int[]{9907194, 9908432});
        result.put("chr21:9907194-9908432:LOC100132288", doubledCDC);
        result.put("chr21:9909047-9909277:LOC100132288", new Seq2covGene("chr21",
                "chr21:9909047-9909277:LOC100132288", new int[]{9909047, 9909277}));
        result.put("chr21:9966322-9966380:LOC100132288", new Seq2covGene("chr21",
                "chr21:9966322-9966380:LOC100132288", new int[]{9966322, 9966380}));
        result.put("chr21:9968516-9968585:LOC100132288", new Seq2covGene("chr21",
                "chr21:9968516-9968585:LOC100132288", new int[]{9968516, 9968585}));
        result.put("chr21:10906743-10907040:TPTE", new Seq2covGene("chr21",
                "chr21:10906743-10907040:TPTE", new int[]{10906743, 10907040}));
        return result;
    }

    @DataProvider(name = "processBamTestData")
    public Object[][] processBamNonAmpPCRData() {
        Map<String, Seq2covGene> resultBed = new LinkedHashMap<String, Seq2covGene>() {{
            put("gene_name", new Seq2covGene("chr1", "gene_name", new int[]{10, 20}));
        }};

        Collection<Gene> expected_nonAmp = new ArrayList<>();
        Collection<Gene> expected_Amp = new ArrayList<>();

        double mdepth = (double) 7 / (20 - 10 + 1);

        expected_nonAmp.add(new GeneX("", "gene_name", "chr1",
                10, 20, "Amplicon", 11, mdepth, 7));
        expected_nonAmp.add(new Gene("", "gene_name", "chr1",
                10, 20, "Whole-Gene", 11, mdepth));

        expected_Amp.add(new GeneX("", "gene_name", "chr1",
                10, 20, "Amplicon", 11, 2.0, 2));
        expected_Amp.add(new Gene("", "gene_name", "chr1",
                10, 20, "Whole-Gene", 11, 2.0));

        return new Object[][]{
                {resultBed, expected_nonAmp, false},
                {resultBed, expected_Amp, true}
        };
    }


    @Test(dataProvider = "processBamTestData")
    public void testProcessBam(Map<String, Seq2covGene> resultBed, Collection<Gene> expected, boolean PCRamplbc)
            throws InterruptedException, ExecutionException, IOException {

        mockSamReaderFactory();

        seq2cov.PCRamplbc = PCRamplbc;
        Collection<Gene> actual = seq2cov.processBam(resultBed);
        assertEquals(expected, actual);
    }

    private void mockSamReaderFactory() {
        final SamReaderFactory factory = PowerMockito.mock(SamReaderFactory.class);
        when(factory.open(any(File.class))).thenReturn(new TestSAMReader());
        when(factory.validationStringency(ValidationStringency.SILENT)).thenReturn(factory);
        mockStatic(SamReaderFactory.class);
        when(SamReaderFactory.makeDefault()).thenReturn(factory);
    }

    @DataProvider(name = "isHGenomeTestData")
    public Object[][] isHGenomeTestData() {
        List<SAMSequenceRecord> records1 = new ArrayList<>();
        List<SAMSequenceRecord> records2 = new ArrayList<>();
        List<SAMSequenceRecord> records3 = new ArrayList<>();

        SAMSequenceRecord record1 = new SAMSequenceRecord("chr1", 2);
        SAMSequenceRecord record2 = new SAMSequenceRecord("grch", 2);

        records1.add(record1);
        records2.add(record2);
        records2.add(record1);
        records3.add(record2);

        return new Object[][]{
                {records1, true},
                {records2, true},
                {records3, false}
        };
    }

    @Test(dataProvider = "isHGenomeTestData")
    public void testIsHGenome(List<SAMSequenceRecord> sequences, boolean expected) {
        assertEquals(expected, seq2cov.isHGenome(sequences));
    }

    @DataProvider(name = "getMrnmTestData")
    public Object[][] getMrnmTestData() {
        SAMRecord equal = new SAMRecord(new SAMFileHeader());
        equal.setReferenceName("equalNames");
        equal.setMateReferenceName("equalNames");

        SAMRecord differ = new SAMRecord(new SAMFileHeader());
        differ.setReferenceName("refName");
        differ.setMateReferenceName("mateName");

        SAMRecord nullMate = mock(SAMRecord.class);
        when(nullMate.getMateReferenceName()).thenReturn(null);
        nullMate.setMateReferenceName("");

        return new Object[][]{
                {equal, "="},
                {differ, "mateName"},
                {nullMate, "*"}
        };
    }

    @Test(dataProvider = "getMrnmTestData")
    public void testGetMrnm(Object record, String expected) {
        assertEquals(expected, getMrnm((SAMRecord) record));
    }

    @DataProvider(name = "getMDOperationLengthTestData")
    public Object[][] getMDOperationLengthTestData() {
        ArrayList<CigarElement> cigarElements = new ArrayList<>();
        cigarElements.add(new CigarElement(1, CigarOperator.M));
        cigarElements.add(new CigarElement(1, CigarOperator.D));
        cigarElements.add(new CigarElement(1, CigarOperator.S));
        Cigar cigar1 = new Cigar(cigarElements);

        ArrayList<CigarElement> cigarElements2 = new ArrayList<>();
        cigarElements.add(new CigarElement(1, CigarOperator.S));
        cigarElements.add(new CigarElement(1, CigarOperator.I));
        cigarElements.add(new CigarElement(1, CigarOperator.H));
        Cigar cigar2 = new Cigar(cigarElements2);
        return new Object[][]{
                {cigar1, 2},
                {cigar2, 0}
        };
    }

    @Test(dataProvider = "getMDOperationLengthTestData")
    public void testGetMDOperationLength(Object cigar, int expected) {
        assertEquals(expected, getMDOperationLength((Cigar) cigar));
    }
}
