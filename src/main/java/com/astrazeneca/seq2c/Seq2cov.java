package com.astrazeneca.seq2c;

import com.astrazeneca.seq2c.Dispatcher.Service;
import htsjdk.samtools.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentLinkedDeque;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.atomic.AtomicInteger;

import static java.lang.Math.max;

/**
 * Class for calculating candidate variance for a given region(s) in an indexed BAM file.
 * */
//perl version: seq2cov.pl
public class Seq2cov {

    private final String bedfileName;
    private final String bamfileName;
    private final String sampleName;
    private final int workersCount;

    boolean PCRamplbc = false;

    /**
     * Creates an instance of Seq2cov class.
     * @param bename bed file name
     * @param bamname bam file name
     * @param spname sample name
     * @param workers number of worker for parallel data processing
     * */
    public Seq2cov(String bename, String bamname, String spname, int workers) {
        bedfileName = bename;
        bamfileName = bamname;
        sampleName = spname;
        this.workersCount = workers;
    }

    /**
     * Creates an instance of Seq2cov that will work in single threaded mode.
     * */
    public Seq2cov(String bename, String bamname, String spname) {
        this(bename, bamname, spname, max(Dispatcher.getThreadsCount(), 1));
    }

    public Collection<Gene> process() throws IOException {
        try {
            return processBam(processBed());
        } catch (InterruptedException | ExecutionException e) {
            throw new IOException(e);
        }
    }

    /**
     * Reads an input bed file and creates a map of gene names and related to them data structures.
     * */
    //perl version: 36 str
    Map<String, Seq2covGene> processBed() throws IOException {
        Map<String, Seq2covGene> resultBed = new LinkedHashMap<>();
        try (BufferedReader bufferReader = new BufferedReader(new FileReader(bedfileName))) {
            String line;
            while ((line = bufferReader.readLine()) != null) {
                String[] parts = line.split("\t");
                String chr = parts[0];
                String gene = parts[3];
                int start = Integer.parseInt(parts[1]);
                int end = Integer.parseInt(parts[2]);
                if (parts.length == 8) {
                    PCRamplbc = true;
                }
                start += 1;
                int[] CDS = { start, end };
                Seq2covGene gn;
                if (resultBed.containsKey(gene)) {
                    gn = resultBed.get(gene);
                    gn.addCDS(CDS);
                } else {
                    gn = new Seq2covGene(chr, gene, CDS);
                }
                resultBed.put(gene, gn);

            }
            return resultBed;
        }
    }

    /**
     * Analyse a given bam file for each region from given map.
     * In parallel mode performs asynchronously by region.
     * @param resultBed map of gene names and related to them data structures from bed file.
     * @return collection of genes with exon coverage
     * */
    //perl version: 66 str
    Collection<Gene> processBam(Map<String, Seq2covGene> resultBed) throws IOException, InterruptedException, ExecutionException {

        Service service = Dispatcher.getService(workersCount);
        Collection<Gene> resultBam = workersCount > 1 ? new ConcurrentLinkedDeque<Gene>() : new ArrayList<Gene>();

        BlockingQueue<Worker> workers = new LinkedBlockingQueue<>();

        workers.add(new Worker(SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bamfileName))));
        boolean genome = isHGenome(workers.peek().sam.getFileHeader().getSequenceDictionary().getSequences());

        //perl version: 68 str
        for (Seq2covGene gn : resultBed.values()) {
            GeneCtx ctx = new GeneCtx(gn, service, workers, resultBam, sampleName, PCRamplbc, genome);
            ctx.process(bamfileName);
        }

        service.await();

        for (Worker worker : workers) {
            worker.sam.close();
        }

        return resultBam;
    }

    boolean isHGenome(List<SAMSequenceRecord> sequences) {
        for (SAMSequenceRecord record : sequences) {
            if (record.getSequenceName().startsWith("chr")) {
                return true; // "hg";
            }
        }
        return false; // "grch";
    }

    /**
     * Data structure for information coming from a bed file.
     * */
    static class Seq2covGene {

        public final List<int[]> CDS = new ArrayList<>();
        public final String chr;
        public final String gene;

        /**
         * @param chr chromosome name
         * @param gene gene name
         * @param CDS array of 2 elements: CDS[0] - start position; CDS[1] - end position
         * */
        public Seq2covGene(String chr, String gene, int[] CDS) {
            this.chr = chr;
            this.gene = gene;
            this.CDS.add(CDS);
        }

        public void addCDS(int[] CDS) {
            this.CDS.add(CDS);
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            Seq2covGene that = (Seq2covGene) o;

            if (CDS.size() != that.CDS.size()) return false;
            for (int i = 0; i < CDS.size(); i++) {
                if (!Arrays.equals(CDS.get(i), that.CDS.get(i))) return false;
            }

            if (!chr.equals(that.chr)) return false;
            return gene.equals(that.gene);
        }

        @Override
        public int hashCode() {
            int result = CDS.hashCode();
            result = 31 * result + chr.hashCode();
            result = 31 * result + gene.hashCode();
            return result;
        }
    }

    /**
     * Class that provides an information for asynchronous workers and collects results of their work.
     * */
    static class GeneCtx {

        final ArrayList<GeneX> gns;
        final AtomicInteger done = new AtomicInteger(0);
        final Seq2covGene gene;
        final Collection<Gene> genes;
        final String sampleName;
        final boolean PCRamplbc;
        final String chr;
        final Service es;
        final BlockingQueue<Worker> workers;

        int start = Integer.MAX_VALUE;
        int end = Integer.MIN_VALUE;
        int length = 0;

        /**
         * @param gene analysing gene
         * @param workers queue of workers
         * @param genes structure for collecting genes with exon coverage
         * @param sampleName sample name
         * @param PCRamplbc indicates that it's PCR amplicon based calling
         * @param isHGenome indicates that bam file presents a human genome
         * */
        public GeneCtx(Seq2covGene gene, Service service, BlockingQueue<Worker> workers, Collection<Gene> genes, String sampleName, boolean PCRamplbc, boolean isHGenome) {
            this.gene = gene;
            this.es = service;
            this.workers = workers;
            this.genes = genes;
            this.sampleName = sampleName;
            this.PCRamplbc = PCRamplbc;
            this.gns = new ArrayList<>(gene.CDS.size());
            for (int i = 0; i < gene.CDS.size(); i++) {
                gns.add(null);
            }
            this.chr = getChrName(isHGenome, gene);
        }

        //perl version: 81 str
        public void process(String bamfileName) throws InterruptedException {
            int workerNumber = 0;
            for (int[] region : gene.CDS) {
                final int rStart = region[0];
                final int rEnd = region[1];
                length += rEnd - rStart + 1;
                if (rStart < start) {
                    start = rStart;
                }
                if (rEnd > end) {
                    end = rEnd;
                }
                Worker worker = workers.poll();
                if (worker == null) {
                    worker = new Worker(SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bamfileName)));
                }
                worker.ctx = this;
                worker.numer = workerNumber++;
                worker.region = region;

                es.submit(worker);
            }

        }

        static String getChrName(boolean genome, Seq2covGene gene) {
            String tchr = gene.chr;

            if (genome && !tchr.startsWith("chr")) {
                tchr += "chr";
            } else if (!genome && tchr.startsWith("chr")) {
                tchr = tchr.replace("chr", "");
            }
            return tchr;
        }

        public void done(Worker worker) throws InterruptedException {
            gns.set(worker.numer,
                    new GeneX(sampleName, gene.gene, chr, worker.region[0], worker.region[1], "Amplicon", worker.length(), getMDepth(worker.exoncov, worker.length()), worker.exoncov));

            int cnt = done.incrementAndGet();
            if (cnt == gene.CDS.size()) {
                int total = 0;
                for (GeneX exoncov : gns) {
                    total += exoncov.exoncov;
                }
                genes.addAll(gns);
                genes.add(new Gene(sampleName, gene.gene, chr, start, end, "Whole-Gene", length, getMDepth(total, length)));
            }
            workers.put(worker);
        }

        double getMDepth(int total, int len) {
            return PCRamplbc ? total : ((double)total / len);
        }

    }

    /**
     * Data structure for genes with exon coverage.
     * */
    static class GeneX extends Gene {

        private int exoncov;

        public GeneX(String sample, String name, String chr, long start, long end, String tag, long len, double mdepth, int exoncov) {
            super(sample, name, chr, start, end, tag, len, mdepth);
            this.exoncov = exoncov;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            if (!super.equals(o)) return false;

            GeneX geneX = (GeneX) o;

            return exoncov == geneX.exoncov;
        }

        @Override
        public int hashCode() {
            int result = super.hashCode();
            result = 31 * result + exoncov;
            return result;
        }
    }

    /**
     * Runnable task for bam file analysis of exon coverage on specific region.
     * */
    static class Worker implements Runnable {

        final SamReader sam;
        int numer;
        GeneCtx ctx;
        int[] region;
        int exoncov = 0;

        public Worker(SamReader sam) {
            this.sam = sam;
        }

        @Override
        public void run() {
            try {
                exoncov = 0;
                int cdsStart = region[0];
                int cdsEnd = region[1];

                //perl version: 88 str
                try (SAMRecordIterator iterator = sam.query(ctx.chr, cdsStart, cdsEnd, false)) {
                    while (iterator.hasNext()) {
                        processRecord(cdsStart, cdsEnd, iterator.next());
                    }
                }
                ctx.done(this);
            } catch (Exception e) {
                try {
                    ctx.done(this);
                } catch (InterruptedException e1) {
                    e1.printStackTrace();
                }
            }
        }

        /**
         * Process sam records and increments exon coverage
         * */
        void processRecord(int cdsStart, int cdsEnd, SAMRecord rec) {
            final Cigar cigar = rec.getCigar();
            final int flags = rec.getFlags();

            if (rec.getStringAttribute(SAMTag.SA.name()) != null && (flags & 0x800) != 0) { // isSupplementaryAlignment
                return;
            }

            final int start = rec.getAlignmentStart();

            final int end = start - 1 + getMDOperationLength(cigar);

            if (ctx.PCRamplbc) { // $opt_a
                double dis = 10;
                double ovlp = 0.95;
                final int ts1;
                final int te1;
                int seqStart = start;
                int seqEnd = end;
                if (!cigar.isEmpty() && cigar.getCigarElement(0).getOperator() == CigarOperator.S && (flags & 0x10) == 0) {
                    final int add = cigar.getCigarElement(0).getLength();
                    if (seqStart - add > cdsStart) {
                        ts1 = seqStart - add;
                    } else {
                        ts1 = cdsStart;
                    }
                    if (seqEnd < cdsEnd) {
                        te1 = seqEnd;
                    } else {
                        te1 = cdsEnd;
                    }

                    if (!(Math.abs((double)ts1 - (double)te1) / ((double)seqEnd - (double)seqStart + add) > ovlp)) {
                        return;
                    }

                } else if (!cigar.isEmpty() && cigar.getCigarElement(cigar.numCigarElements() - 1).getOperator() == CigarOperator.S && (flags & 0x10) != 0) {
                    int add = cigar.getCigarElement(cigar.numCigarElements() - 1).getLength();
                    if (seqStart > cdsStart) {
                        ts1 = seqStart;
                    } else {
                        ts1 = cdsStart;
                    }
                    if (seqEnd + add < cdsEnd) {
                        te1 = seqEnd + add;
                    } else {
                        te1 = cdsEnd;
                    }
                    if (!(Math.abs((double)te1 - (double)ts1) / ((double)seqEnd + (double)add - seqStart) > ovlp)) {
                        return;
                    }
                } else {
                    final String mrnm = getMrnm(rec);

                    if ("=".equals(mrnm) && rec.getInferredInsertSize() != 0) {
                        int pNext = rec.getInferredInsertSize();
                        int rNext = rec.getMateAlignmentStart();
                        if (pNext > 0) {
                            seqEnd = seqStart + pNext - 1;
                        } else {
                            seqStart = rNext;
                            seqEnd = rNext - pNext - 1;
                        }
                    }
                    if (seqStart > cdsStart) {
                        ts1 = seqStart;
                    } else {
                        ts1 = cdsStart;
                    }
                    if (seqEnd < cdsEnd) {
                        te1 = seqEnd;
                    } else {
                        te1 = cdsEnd;
                    }
                    if (!((Math.abs((double)seqStart - (double)cdsStart) <= dis && Math.abs((double)seqEnd - (double)cdsEnd) <= dis) && Math.abs(((double)ts1 - (double)te1) / ((double)seqEnd - (double)seqStart)) > ovlp)) {
                        return;
                    }

                }
                exoncov++;
            } else {
                final double alEnd;
                if (cdsEnd > end) {
                    alEnd = end;
                } else {
                    alEnd = cdsEnd;
                }

                final double alStart;
                if (cdsStart > start) {
                    alStart = cdsStart;
                } else {
                    alStart = start;
                }

                final double alignen = (alEnd - alStart) + 1;
                exoncov += alignen;
            }
        }

        int length() {
            return region[1] - region[0] + 1;
        }

    }

    static String getMrnm(SAMRecord record) {
        if (record.getMateReferenceName() == null) {
            return "*";
        }

        if (record.getReferenceName().equals(record.getMateReferenceName())) {
            return "=";
        }
        return record.getMateReferenceName();
    }

    //perl version: 94 str
    static int getMDOperationLength(Cigar cigar) {
        int length = 0;
        for (CigarElement ce : cigar.getCigarElements()) {
            if (ce.getOperator() == CigarOperator.M || ce.getOperator() == CigarOperator.D) {
                length += ce.getLength();
            }
        }
        return length;
    }
}
