package com.astrazeneca.seq2c.input;

import com.astrazeneca.seq2c.SampleStatistics;
import com.astrazeneca.seq2c.Sample;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * Created by Mariia_Zueva on 4/4/2016.
 */
public class StatisticsFactory extends FileStoredDataFactory {

    @Override
    public FileStoredData createObjectFromLine(String line) {
        if (line == null) return null;

        String[] sampleLines = line.split("\\t");
        if (sampleLines.length < 8) return null;

        String sample = sampleLines[0];
        String gene = sampleLines[1];
        Sample sampleObj = readSample(line);

        SampleStatistics statistics = new SampleStatistics(sample);
        HashMap<String, List<Sample>> genes = new HashMap<>();
        List<Sample> list = new ArrayList<>();
        list.add(sampleObj);
        genes.put(gene, list);
        statistics.setGenes(genes);
        return statistics;
    }

    private Sample readSample(String currentLine) {

        if (currentLine == null) return null;

        String[] sampleLines = currentLine.split("\\s+");
        if (sampleLines.length < 8) return null;

        String sampleName = sampleLines[0];
        String gene = sampleLines[1];
        String chr = sampleLines[2];
        int start = Integer.parseInt(sampleLines[3]);
        int end = Integer.parseInt(sampleLines[4]);
        int len = Integer.parseInt(sampleLines[5]);

        String key = gene;
        Sample sample = new Sample(key, sampleName, chr, start, end, gene, len, 0.0);
        sample.setNorm3(Double.parseDouble(sampleLines[6]));
        sample.setNorm3s(Double.parseDouble(sampleLines[7]));
        return sample;
    }
}
