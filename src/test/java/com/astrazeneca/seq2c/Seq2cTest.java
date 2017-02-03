package com.astrazeneca.seq2c;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.mockito.Mockito;
import org.powermock.api.mockito.PowerMockito;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.nio.file.Files;
import java.util.*;

public class Seq2cTest {
    @DataProvider(name = "printStatDataProvider")
    public Object[][] printStatDataProvider() {
        Map<String, Long> sam2bam = new LinkedHashMap<String, Long>(){{
            put("P-E0001004-ff", 125L);
            put("P-E0001004-ffpe", 1250L);
            put("P-E0001004-normal", 124L);
        }};

        ArrayList<String> expected = new ArrayList<String>(){{
            add("P-E0001004-ff\t125");
            add("P-E0001004-ffpe\t1250");
            add("P-E0001004-normal\t124");
        }};
        return new Object[][]{
                {sam2bam, expected}
        };
    }

    @DataProvider(name = "printCoverageProvider")
    public Object[][] printCoverageProvider() {
        List<Gene> sam2bam = new ArrayList<Gene>(){{
            add(
                new Gene("P-E0001004-normal", "BAGE5:11057796-11058160:BAGE",
                    "chr21", 11098633, 11098768,
                    "Amplicon", 136, 133.77
                )
            );

            add(new Gene("P-E0001004-normal", "chr21:10910307-10910399:TPTE",
                        "chr21", 10910307, 10910399,
                        "Whole-Gene", 93, 688.68
                        )
            );

            add(
                new Gene("P-E0001004-ffpe", "chr21:10916370-10916475:TPTE",
                        "chr21", 10916370, 10916475,
                        "Whole-Gene", 106, 127.25
                )
            );
        }};

        ArrayList<String> expected = new ArrayList<String>(){{
            add("Sample\tGene\tChr\tStart\tEnd\tTag\tLength\tMeanDepth");
            add("P-E0001004-normal\tBAGE5:11057796-11058160:BAGE\tchr21\t11098633\t11098768\tAmplicon\t136\t133.77");
            add("P-E0001004-normal\tchr21:10910307-10910399:TPTE\tchr21\t10910307\t10910399\tWhole-Gene\t93\t688.68");
            add("P-E0001004-ffpe\tchr21:10916370-10916475:TPTE\tchr21\t10916370\t10916475\tWhole-Gene\t106\t127.25");
        }};
        return new Object[][]{
                {sam2bam, expected}
        };
    }

    @DataProvider(name = "configDataProvider")
    public Object[][] configDataProvider() throws ParseException {
        CommandLine cmdAvailableProcessorThreadCount = getCmdWithParameter("i", null);
        CommandLine cmd1Thread = getCmdWithParameter("i", 1);

        return new Object[][]{
                {cmdAvailableProcessorThreadCount, Runtime.getRuntime().availableProcessors()},
                {cmd1Thread, 1}
        };
    }

    @DataProvider(name = "runPartDataProvider")
    public Object[][] runPartDataProvider() throws ParseException {
        CommandLine cmd0RunPart = getCmdWithParameter("r", null);
        CommandLine cmd1RunPart = getCmdWithParameter("r", 1);
        CommandLine cmd2RunPart = getCmdWithParameter("r", 2);

        return new Object[][]{
                {cmd0RunPart, 0},
                {cmd1RunPart, 1},
                {cmd2RunPart, 2}
        };
    }

    @Test(dataProvider = "printStatDataProvider")
    public void testPrintStat(Map<String, Long> sam2bam, List<String> expected) throws Exception {
        File tempFile = File.createTempFile("readstat", "txt");
        Seq2c.printStat(sam2bam, tempFile.getAbsolutePath());
        try (BufferedReader reader = new BufferedReader(new FileReader(tempFile.getAbsolutePath()))) {
            for (String expectedString : expected)
                Assert.assertEquals(reader.readLine(), expectedString);
        } finally {
            Files.deleteIfExists(tempFile.toPath());
        }
    }

    @Test(dataProvider = "printCoverageProvider")
    public void testPrintCoverage(List<Gene> genes, List<String> expected) throws Exception {
        File tempFile = File.createTempFile("readstat", "txt");
        Seq2c.printCoverage(genes, tempFile.getAbsolutePath(), true);
        try (BufferedReader reader = new BufferedReader(new FileReader(tempFile.getAbsolutePath()))) {
            for (String expectedString : expected)
                Assert.assertEquals(reader.readLine(), expectedString);
        } finally {
            Files.deleteIfExists(tempFile.toPath());
        }
    }

    @Test(dataProvider = "configDataProvider")
    public void testGetThreadsCount(CommandLine cmd, int expected) throws Exception {
        Assert.assertEquals(Seq2c.getThreadsCount(cmd),expected);
    }

    @Test(dataProvider = "runPartDataProvider")
    public void testGetRunPart(CommandLine cmd, int expected) throws Exception {
        Assert.assertEquals(Seq2c.getRunPart(cmd),expected);
    }

    @Test
    public void testBuildCmdOptions() throws Exception {
        Options options = Seq2c.buildCmdOptions();
        Assert.assertEquals(options.getOption("i").getArgName(), "number of threads");
        Assert.assertEquals(options.getOption("i").getType(), Number.class);
        Assert.assertEquals(options.getOption("i").isRequired(), false);

        Assert.assertEquals(options.getOption("r").getArgName(), "run part");
        Assert.assertEquals(options.getOption("r").getType(), Number.class);
        Assert.assertEquals(options.getOption("r").isRequired(), false);

        Assert.assertEquals(options.getOption("h").getArgName(), "help");
        Assert.assertEquals(options.getOption("h").getType(), Number.class);
        Assert.assertEquals(options.getOption("h").isRequired(), false);
    }

    private CommandLine getCmdWithParameter(String param, Object value) throws ParseException {
        CommandLine cmd1Thread = Mockito.mock(CommandLine.class);
        PowerMockito
                .when(cmd1Thread.hasOption(param))
                .thenReturn(true);

        PowerMockito.when(cmd1Thread.getParsedOptionValue(param)).thenReturn(value);
        return cmd1Thread;
    }
}
