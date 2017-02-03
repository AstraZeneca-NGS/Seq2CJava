package com.astrazeneca.seq2c;

import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.concurrent.LinkedBlockingQueue;

import static com.astrazeneca.seq2c.Seq2cov.*;
import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertTrue;

public class TestGeneCtx {

    @DataProvider(name = "getChrNameTestData")
    public Object[][] getChrNameTestData() {
        return new Object[][]{
                {new Seq2covGene("chr1", "", new int[]{}), true, "chr1"},
                {new Seq2covGene("1", "", new int[]{}), true, "1chr"},
                {new Seq2covGene("chr1", "", new int[]{}), false, "1"}
        };
    }

    @Test(dataProvider = "getChrNameTestData")
    public void testGetChrName(Seq2covGene gene, boolean genome, String expected) {
        assertEquals(expected, GeneCtx.getChrName(genome, gene));
    }

    @DataProvider(name = "doneTestData")
    public Object[][] doneTestData() {
        Worker testWorker = new Worker(new TestSAMReader());
        testWorker.region = new int[]{0, 10};
        testWorker.numer = 0;
        testWorker.exoncov = 2;
        Seq2covGene gene1 = new Seq2covGene("chr1", "testGene", testWorker.region);
        Seq2covGene gene2 = new Seq2covGene("chr1", "testGene", testWorker.region);
        gene2.addCDS(new int[]{15, 20});
        GeneX geneX = new GeneX("testSample", gene1.gene, gene1.chr,
                testWorker.region[0], testWorker.region[1],
                "Amplicon", 11, (double) 2 / 11, testWorker.exoncov);
        Gene resultGene = new Gene("testSample", gene1.gene, gene1.chr,
                Integer.MAX_VALUE, Integer.MIN_VALUE, "Whole-Gene", 0, Double.POSITIVE_INFINITY);
        return new Object[][]{
                {testWorker, gene1, geneX, resultGene, 0},
                {testWorker, gene2, geneX, null, 0}
        };
    }

    @Test(dataProvider = "doneTestData")
    public void testDone(Worker worker,
                         Seq2covGene gene,
                         GeneX expectedGeneX,
                         Gene expectedGene,
                         int index) throws InterruptedException {

        GeneCtx ctx = new GeneCtx(gene, Dispatcher.getService(1),
                new LinkedBlockingQueue<Seq2cov.Worker>(), new ArrayList<Gene>(),
                "testSample", false, true);
        ctx.done(worker);

        assertEquals(ctx.gns.get(index), expectedGeneX);
        if (expectedGene != null) {
            ctx.genes.contains(expectedGene);
            assertTrue(ctx.genes.contains(expectedGene));
            assertTrue(ctx.genes.contains(expectedGeneX));
        }
        assertTrue(ctx.workers.contains(worker));
    }

    @DataProvider(name = "getMdepthTestData")
    public Object[][] getMdepthTestData() {
        return new Object[][]{
                {false, 10, 5, 2.0},
                {true, 10, 5, 10.0}
        };
    }

    @Test(dataProvider = "getMdepthTestData")
    public void testGetMdepth(boolean PCRamplbc, int total, int length, double result) {
        GeneCtx  ctx = new GeneCtx(new Seq2covGene("", "", new int[]{}),
                Dispatcher.getService(1),
                new LinkedBlockingQueue<Seq2cov.Worker>(),
                new ArrayList<Gene>(),
                "testSample",
                PCRamplbc,
                true);
        assertEquals(ctx.getMDepth(total, length), result);
    }

}
