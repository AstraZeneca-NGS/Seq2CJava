package com.astrazeneca.seq2c;

import java.io.*;
import java.util.LinkedHashMap;
import java.util.Map;

import htsjdk.samtools.AbstractBAMFileIndex;
import htsjdk.samtools.BAMIndexMetaData;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

//perl version: bam2reads
public class Bam2Reads {

    /**
     * Reads a map of sample name/bam filename from fileIn. For each bam file reads its indices, using samtools and gets a number of
     * aligned records from index information. Number of aligned records for each file is recorded into output file, Format of
     * output file: one line per bam file each line: "sample name" \t "number of aligned records"
     *
     * @param fileIn
     *            contains list of sample names and bam files to process, each line in file represent one file
     * @throws IOException
     *
     */

    //perl version: 12-21 str
    public static Map<String, Long> printStatsToFile(String fileIn) throws IOException {
        Map<String, String> files = parseFile(fileIn);
        Map<String, Long> result = new LinkedHashMap<>();

        for (Map.Entry<String, String> entry : files.entrySet()) {
            try (SamReader sam = SamReaderFactory.makeDefault().open(new File(entry.getValue()))) {

                SamReader.Indexing ind = sam.indexing();
                AbstractBAMFileIndex index = (AbstractBAMFileIndex)ind.getIndex();
                long count = 0;
                for (int i = 0; i < index.getNumberOfReferences(); i++) {
                    BAMIndexMetaData meta = index.getMetaData(i);
                    count += meta.getAlignedRecordCount();
                }
                result.put(entry.getKey(), count);
            }
        }
        return result;
    }

    /**
     * Reads input file line by line and creates a map containing sample name/bam file name
     *
     * @param fileName
     *            contains list of sample names and bam files to process, each line in file represent one file
     * @return map with sample names/bam file names
     * @throws IOException
     */
    //perl version: 14 str
    public static Map<String, String> parseFile(String fileName) throws IOException {
        Map<String, String> map = new LinkedHashMap<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(fileName))) {
            String currentLine;
            while ((currentLine = reader.readLine()) != null) {
                String[] samples = currentLine.split("\\s+");
                if (samples.length == 2) {
                    map.put(samples[0], samples[1]);
                }
            }
        }
        return map;
    }
}
