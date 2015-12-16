package seq2c;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

import htsjdk.samtools.*;

/**
 *
 * @author Petr_Rastegaev
 */
public class Seq2cov {

    //public static void main(String[] args)
    //{
       //processBed(fileName);
     //}

    boolean PCRamplbc = false;
    private String bedfileName;
    private String bamfileName;
    private String sampleName;
    private Map<String,Seq2cov_gene> result_bed = new LinkedHashMap<>();
    private List<Gene> result_bam = new ArrayList<>();

    public Seq2cov(String bename, String bamname, String spname){
         bedfileName = bename;
         bamfileName = bamname;
         sampleName = spname;
    }
    public List<Gene> process() throws IOException{
        processBed();
        processBam();
        return result_bam;
    }

    private Map<String, Seq2cov_gene> processBed() throws IOException {
        // System.out.println("Reading File from Java code");
        // Name of the file
        // Create object of FileReader
        // Variable to hold the one line data
        try (BufferedReader bufferReader = new BufferedReader(new FileReader(bedfileName))) {
            // Variable to hold the one line data
            String line;
            // Read file line by line and print on the console
            while ((line = bufferReader.readLine()) != null) {
                String[] parts = line.split("\t");
                String chr = parts[0];
                String gene = parts[3];
                int start = Integer.parseInt(parts[1]);
                int end = Integer.parseInt(parts[2]);
                if (parts.length == 8) {
                    PCRamplbc = true;
                    // System.out.println("Its opt_a!");
                }
                start += 1;
                Integer[] CDS = { start, end };
                Seq2cov_gene gn;
                if (result_bed.containsKey(gene)) {
                    gn = result_bed.get(gene);
                    gn.addCDS(CDS);
                } else {
                    gn = new Seq2cov_gene(chr, gene, CDS);
                }
                result_bed.put(gene, gn);

            }
            return result_bed;
        }
    }

    private List<Gene> processBam() throws IOException {
        long time = 0;

        try (SamReader sam = SamReaderFactory.makeDefault().open(new File(bamfileName))) {
            List<SAMSequenceRecord> sequences = sam.getFileHeader().getSequenceDictionary().getSequences();
            boolean genome;
            if (!sequences.isEmpty() &&  sequences.get(0).getSequenceName().startsWith("chr")) {
                genome = true; // "hg";
            } else {
                genome = false; // "grch";
            }

            for (Seq2cov_gene gn : result_bed.values()) {
                int exoncov;
                int total = 0;
                int geneLength = 0;
                int geneStart = 500000000;
                int geneEnd = 0;
                String tchr = gn.chr;
                if (genome && !tchr.startsWith("chr")) {
                    tchr += "chr";
                } else if (!genome && tchr.startsWith("chr")) {
                    tchr = tchr.replace("chr", "");
                }
                for (Integer[] CDS_arr : gn.CDS) {
                    exoncov = 0;
                    final int cdsStart = CDS_arr[0];
                    final int cdsEnd = CDS_arr[1];
                    geneLength += cdsEnd - cdsStart + 1;
                    if (cdsStart < geneStart) {
                        geneStart = cdsStart;
                    }
                    if (cdsEnd > geneEnd) {
                        geneEnd = cdsEnd;
                    }
                    SAMRecordIterator iterator = sam.query(tchr, cdsStart, cdsEnd, false);
                    while (iterator.hasNext()) {
                        long t = System.currentTimeMillis();
                        final SAMRecord rec = iterator.next();
                        Flags flags = new Flags(rec.getFlags());
                        Cigar cigar = rec.getCigar();

                        if (rec.getStringAttribute(SAMTag.SA.name()) != null && flags.isSupplementaryAlignment()) {
                            continue;
                        }

                        int start = rec.getAlignmentStart();

                        int end = start - 1 + getMDOperationLength(cigar);

                        // System.out.println(gn.gene + " " + start + " " + end );
                        if (PCRamplbc) { // $opt_a
                            double dis = 10;
                            double ovlp = 0.95;
                            int ts1;
                            int te1;
                            int seqStart = start;
                            int seqEnd = end;
                            if (!cigar.isEmpty() && cigar.getCigarElement(0).getOperator() == CigarOperator.S && !flags.isReverseStrand()) {
                                int add = cigar.getCigarElement(0).getLength();
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
                                // System.out.println("entered if " +cdsStart+" "+cdsEnd + tchr);
                                if (!(Math.abs((double)ts1 - (double)te1) / ((double)seqEnd - (double)seqStart + add) > ovlp)) {
                                    continue;
                                }

                            } else if (!cigar.isEmpty() && cigar.getCigarElement(cigar.numCigarElements() -1 ).getOperator() == CigarOperator.S && flags.isReverseStrand()) {
                                int add = cigar.getCigarElement(cigar.numCigarElements() -1).getLength();
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
                                // System.out.println("entered elseif " +cdsStart+" "+cdsEnd + tchr);
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
                                // System.out.println("entered else " +cdsStart+" "+cdsEnd + tchr);
                                if (!((Math.abs((double)seqStart - (double)cdsStart) <= dis && Math.abs((double)seqEnd - (double)cdsEnd) <= dis) && Math.abs(((double)ts1 - (double)te1) / ((double)seqEnd - (double)seqStart)) > ovlp)) {
                                    continue;
                                }

                            }
                            exoncov++;
                            total++;
                        } else {
                            double alEnd;
                            double alStart;
                            if (cdsEnd > end) {
                                alEnd = end;
                            } else {
                                alEnd = cdsEnd;
                            }

                            if (cdsStart > start) {
                                alStart = cdsStart;
                            } else {
                                alStart = start;
                            }

                            double alignen = (alEnd - alStart) + 1;
                            exoncov += alignen;
                            total += alignen;
                        }
                        time += (System.currentTimeMillis() - t);

                    } // iterator
                    iterator.close();
                    // System.out.println(exoncov + " " + total);
                    Gene sqr;
                    if (PCRamplbc) {
                        sqr = new Gene(sampleName, gn.gene, tchr, cdsStart, cdsEnd, "Amplicon", cdsEnd - cdsStart + 1, exoncov);
                    } else {
                        sqr = new Gene(sampleName, gn.gene, tchr, cdsStart, cdsEnd, "Amplicon", cdsEnd - cdsStart + 1, (double)exoncov / (double)(cdsEnd - cdsStart + 1));
                    }
                    result_bam.add(sqr);

                } // for loop2
                Gene sqr;
                if (PCRamplbc) {
                    sqr = new Gene(sampleName, gn.gene, tchr, geneStart, geneEnd, "Whole-Gene", geneLength, total);
                } else {
                    sqr = new Gene(sampleName, gn.gene, tchr, geneStart, geneEnd, "Whole-Gene", geneLength, (double)total / geneLength);
                }
                result_bam.add(sqr);
            } // for loop1
            System.err.println("TIME: " + time);
            return result_bam;
        }
        //
    }

    private void processRecords(Collection<SAMRecord> records, int cdsStart, int cdsEnd) {
        int exoncov = 0;
        int total = 0;

        for (SAMRecord rec : records) {

            Flags flags = new Flags(rec.getFlags());
            Cigar cigar = rec.getCigar();

            if (rec.getStringAttribute(SAMTag.SA.name()) != null && flags.isSupplementaryAlignment()) {
                continue;
            }

            int start = rec.getAlignmentStart();

            int end = start - 1 + getMDOperationLength(cigar);

            // System.out.println(gn.gene + " " + start + " " + end );
            if (PCRamplbc) { // $opt_a
                double dis = 10;
                double ovlp = 0.95;
                int ts1;
                int te1;
                int seqStart = start;
                int seqEnd = end;
                if (!cigar.isEmpty() && cigar.getCigarElement(0).getOperator() == CigarOperator.S && !flags.isReverseStrand()) {
                    int add = cigar.getCigarElement(0).getLength();
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
                    // System.out.println("entered if " +cdsStart+" "+cdsEnd + tchr);
                    if (!(Math.abs((double)ts1 - (double)te1) / ((double)seqEnd - (double)seqStart + add) > ovlp)) {
                        continue;
                    }

                } else if (!cigar.isEmpty() && cigar.getCigarElement(cigar.numCigarElements() - 1).getOperator() == CigarOperator.S && flags.isReverseStrand()) {
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
                    // System.out.println("entered elseif " +cdsStart+" "+cdsEnd + tchr);
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
                    // System.out.println("entered else " +cdsStart+" "+cdsEnd + tchr);
                    if (!((Math.abs((double)seqStart - (double)cdsStart) <= dis && Math.abs((double)seqEnd - (double)cdsEnd) <= dis) && Math.abs(((double)ts1 - (double)te1) / ((double)seqEnd - (double)seqStart)) > ovlp)) {
                        continue;
                    }

                }
                exoncov++;
                total++;
            } else {
                double alEnd;
                double alStart;
                if (cdsEnd > end) {
                    alEnd = end;
                } else {
                    alEnd = cdsEnd;
                }

                if (cdsStart > start) {
                    alStart = cdsStart;
                } else {
                    alStart = start;
                }

                double alignen = (alEnd - alStart) + 1;
                exoncov += alignen;
                total += alignen;
            }
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

    private int getMDOperationLength(Cigar cigar) {
        int length = 0;
        for (CigarElement ce : cigar.getCigarElements()) {
            if (ce.getOperator() == CigarOperator.M || ce.getOperator() == CigarOperator.D) {
                length += ce.getLength();
            }
        }
        return length;
    }

    private static class BatchBuffer<T> {

        private final int bufferSize;
        private final BatchConsumer<T> consumer;
        private Collection<T> batch;

        public BatchBuffer(int bufferSize, BatchConsumer<T> consumer) {
            this.bufferSize = bufferSize;
            this.consumer = consumer;
            batch = buildBatch();
        }

        public void add(T e) {
            batch.add(e);
            if (batch.size() == bufferSize) {
                flush();
            }
        }

        public void flush() {
            consumer.consume(batch);
            batch = buildBatch();
        }

        protected Collection<T> buildBatch() {
            return new ArrayList<T>(bufferSize);
        }

    }

    private static interface BatchConsumer<T> {

        void consume(Collection<T> batch);
    }

    private static class Flags {
        private int flag;

        public Flags(int flag) {
            this.flag = flag;
        }

        /**
         * @return indicates that the corresponding alignment line is part of a chimeric alignment.
         */
        public boolean isSupplementaryAlignment() {
            return (flag & 0x800) != 0;
        }

        /**
         *
         * @return next segment in the template is unmapped
         */
        public boolean isUnmappedMate() {
            return (flag & 0x8) != 0;
        }

        /**
         *
         * @return sequence is reverse complemented
         */
        public boolean isReverseStrand() {
            return (flag & 0x10) != 0;
        }

        /**
         * Bit 0x100 marks the alignment not to be used in certain analyses when the tools in use are aware of this bit.
         * It is typically used to flag alternative mappings when multiple mappings are presented in a SAM.
         * @return Is secondary alignment
         */
        public boolean isNotPrimaryAlignment() {
            return (flag & 0x100) != 0;
        }

    }
}
