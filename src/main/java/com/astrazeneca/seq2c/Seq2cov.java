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
import static java.lang.Math.min;

public class Seq2cov {

    private final String bedfileName;
    private final String bamfileName;
    private final String sampleName;
    private final int workersCount;

    private boolean PCRamplbc = false;

    public Seq2cov(String bename, String bamname, String spname, int workers) {
        bedfileName = bename;
        bamfileName = bamname;
        sampleName = spname;
        this.workersCount = workers;
    }

    public Seq2cov(String bename, String bamname, String spname) {
        this(bename, bamname, spname, max(min(Dispatcher.getThreadsCount(), 8), 1));
    }

    public Collection<Gene> process() throws IOException {
        try {
            return processBam(processBed());
        } catch (InterruptedException | ExecutionException e) {
            throw new IOException(e);
        }
    }

    private Map<String, Seq2covGene> processBed() throws IOException {
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

    private Collection<Gene> processBam(Map<String, Seq2covGene> resultBed) throws IOException, InterruptedException, ExecutionException {

        Service service = Dispatcher.getService(workersCount);
        Collection<Gene> resultBam = workersCount > 1 ? new ConcurrentLinkedDeque<Gene>() : new ArrayList<Gene>();

        BlockingQueue<Worker> workers = new LinkedBlockingQueue<>();
        workers.add(new Worker(SamReaderFactory.makeDefault().open(new File(bamfileName))));
        boolean genome = isHGenome(workers.peek().sam.getFileHeader().getSequenceDictionary().getSequences());

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

    private boolean isHGenome(List<SAMSequenceRecord> sequences) {
        for (SAMSequenceRecord record : sequences) {
            if (record.getSequenceName().startsWith("chr")) {
                return true; // "hg";
            }
        }
        return false; // "grch";
    }

    private static class Seq2covGene {

        public final List<int[]> CDS = new ArrayList<>();
        public final String chr;
        public final String gene;

        public Seq2covGene(String chr, String gene, int[] CDS) {
            this.chr = chr;
            this.gene = gene;
            this.CDS.add(CDS);
        }

        public void addCDS(int[] CDS) {
            this.CDS.add(CDS);
        }
    }

    private static class GeneCtx {

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
                    worker = new Worker(SamReaderFactory.makeDefault().open(new File(bamfileName)));
                }
                worker.ctx = this;
                worker.numer = workerNumber++;
                worker.region = region;

                es.submit(worker);
            }

        }

        private static String getChrName(boolean genome, Seq2covGene gene) {
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

        private double getMDepth(int total, int len) {
            return PCRamplbc ? total : ((double)total / len);
        }

    }

    private static class GeneX extends Gene {

        private int exoncov;

        public GeneX(String sample, String name, String chr, long start, long end, String tag, long len, double mdepth, int exoncov) {
            super(sample, name, chr, start, end, tag, len, mdepth);
            this.exoncov = exoncov;
        }

    }

    private static class Worker implements Runnable {

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

                try (SAMRecordIterator iterator = sam.query(ctx.chr, cdsStart, cdsEnd, false)) {
                    while (iterator.hasNext()) {
                        final SAMRecord rec = iterator.next();
                        final Cigar cigar = rec.getCigar();
                        final int flags = rec.getFlags();

                        if (rec.getStringAttribute(SAMTag.SA.name()) != null && (flags & 0x800) != 0) { // isSupplementaryAlignment
                            continue;
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
                                    continue;
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
                                    continue;
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
                                    continue;
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

        int length() {
            return region[1] - region[0] + 1;
        }

    }

    private static String getMrnm(SAMRecord record) {
        if (record.getMateReferenceName() == null) {
            return "*";
        }

        if (record.getReferenceName().equals(record.getMateReferenceName())) {
            return "=";
        }
        return record.getMateReferenceName();
    }

    private static int getMDOperationLength(Cigar cigar) {
        int length = 0;
        for (CigarElement ce : cigar.getCigarElements()) {
            if (ce.getOperator() == CigarOperator.M || ce.getOperator() == CigarOperator.D) {
                length += ce.getLength();
            }
        }
        return length;
    }
}
