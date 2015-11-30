/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package seq2c;
import htsjdk.samtools.AbstractBAMFileIndex;
import htsjdk.samtools.BAMIndexMetaData;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.*;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 *
 * @author Petr_Rastegaev
 */
public class Bam2Reads {
    
    /**
     * Reads a map of sample name/bam filename from fileIn. For each bam file reads its indices,
     * using samtools and gets a number of aligned records from index information.
     * Number of aligned records for each file is recorded into output file,
     * Format of output file:
     * one line per bam file
     * each line: "sample name" \t "number of aligned records"
     *
     * @param fileIn  contains list of sample names and bam files to process,
     *                each line in file represent one file
     * 
     */

    public static Map<String, Long> printStatsToFile(String fileIn) {
        Map<String, String> files = parseFile(fileIn);
        Map<String, Long> result = new LinkedHashMap<>();
        
            for (Map.Entry<String, String> entry : files.entrySet()) {
                //writer.write(entry.getKey() + "\t");
                //System.out.println(entry.getKey() + "\t" + entry.getValue());
                try (SamReader sam = SamReaderFactory.makeDefault().open(new File(entry.getValue()))) {

                    SamReader.Indexing ind = sam.indexing();
                    AbstractBAMFileIndex index = (AbstractBAMFileIndex) ind.getIndex();
                    long count = 0;
                    for (int i = 0; i < index.getNumberOfReferences(); i++) {
                        BAMIndexMetaData meta = index.getMetaData(i);
                        count += meta.getAlignedRecordCount();
                    }
                    System.out.println(entry.getKey()+" "+ count);
                    result.put(entry.getKey(), count);
                    //writer.write(count + "\n");
                    //writer.flush();
                } catch (IOException e) {
                    System.err.println("Cannot read file " + entry.getValue());
                }
            }
        return result;
    }

    /**
     * Reads input file line by line and creates a map containing
     * sample name/bam file name
     *
     * @param fileName contains list of sample names and bam files to process,
     *                 each line in file represent one file
     * @return map with sample names/bam file names
     */
    private static Map<String, String> parseFile(String fileName) {
        Map<String, String> map = new LinkedHashMap<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(new File(fileName)))) {
            String currentLine;
            while ((currentLine = reader.readLine()) != null) {
                String[] samples = currentLine.split("\\s+");
                if (samples.length == 2) {
                    map.put(samples[0], samples[1]);
                }
            }
        } catch (IOException e) {
            System.err.println("Cannot open file " + fileName);
        }
        return map;
    }
}
