package com.astrazeneca.seq2c;

import com.astrazeneca.seq2c.input.Sample;
import org.junit.Ignore;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertEquals;

public class Lr2geneTest {

    @Test
    public void checkBPTestCOMMD1() throws Exception {
        List<Sample> samples = getSampleCOMMD1();
        Sig sig = new Lr2gene().checkBP(samples);
        assertEquals("12", sig.getSigStr());
        assertEquals("BP", sig.getBp());
        assertEquals("Del", sig.getCn());
        assertEquals(1.0, sig.getBpi(), 0.001);
        assertEquals(4, sig.getTotal());
        assertEquals(-4.23, sig.getSiglr(), 0.001);
        assertEquals(3.813, sig.getSigdiff(), 0.001);
        assertEquals("4", sig.getSigseg());
    }

    @Test
    public void checkBPTestPKN2() throws Exception {

        List<Sample> samples = getSamplePKN2();
        Sig sig = new Lr2gene().checkBP(samples);
        assertEquals(-1.0, sig.getSig(), 0.001);
        assertEquals(15, sig.getTotal());
        assertEquals(0.0, sig.getSiglr(), 0.001);
        assertEquals(0.0, sig.getSigdiff(), 0.001);
        assertEquals("", sig.getCn());
    }

    @Test
    public void checkBPTestNotAGene4865() throws Exception {

        List<Sample> samples = getSampleNotAGene4865();
        Sig sig = new Lr2gene().checkBP(samples);
        assertEquals(-1.0, sig.getSig(), 0.001);
        assertEquals(28, sig.getTotal());
        assertEquals(0.0, sig.getSiglr(), 0.001);
        assertEquals(0.0, sig.getSigdiff(), 0.001);
        assertEquals("", sig.getCn());
    }

    @Test
    public void checkBPTestNotAGene4998() throws Exception {

        List<Sample> samples = getSampleNotAGene4998();
        Sig sig = new Lr2gene().checkBP(samples);
        assertEquals(0.0004177, sig.getSig(), 0.000001);
        assertEquals(17, sig.getTotal());
        assertEquals(10, sig.getBpi());
        assertEquals(2.61597, sig.getSigdiff(), 0.0001);
        assertEquals(-1.41222, sig.getSiglr(), 0.0001);
        assertEquals("Del", sig.getCn());
        assertEquals("9,10,11,12,13,14,15,16,17", sig.getSigseg());
    }

    @Test
    public void checkBPTestNotAGene6082() throws Exception {

        List<Sample> samples = getSampleNotAGene6082();
        Sig sig = new Lr2gene().checkBP(samples);
        assertEquals("45", sig.getSigStr());
        assertEquals("BP", sig.getBp());
        assertEquals("Del", sig.getCn());
        assertEquals(3, sig.getBpi());
        assertEquals(5, sig.getTotal());
        assertEquals(-0.757, sig.getSiglr(), 0.001);
        assertEquals(1.251666, sig.getSigdiff(), 0.001);
        assertEquals("2,3,4", sig.getSigseg());
    }

    @Test
    public void checkBPTestAACS() throws Exception {

        List<Sample> samples = getSampleAACS();
        Sig sig = new Lr2gene().checkBP(samples);
        assertEquals(-1, sig.getSig(), 0.001);
        assertEquals("", sig.getBp());
        assertEquals("", sig.getCn());
        assertEquals(0, sig.getBpi());
        assertEquals(16, sig.getTotal());
        assertEquals(0.0, sig.getSiglr(), 0.001);
        assertEquals(0.0, sig.getSigdiff(), 0.001);
        assertEquals("", sig.getSigseg());
    }

    @Test
    public void checkBPTestAASDHPPT() throws Exception {

        List<Sample> samples = getSampleAASDHPPT();
        Sig sig = new Lr2gene().checkBP(samples);
        assertEquals("6.6", sig.getSigStr());
        assertEquals("BP", sig.getBp());
        assertEquals("Amp", sig.getCn());
        assertEquals(2, sig.getBpi());
        assertEquals(6, sig.getTotal());
        assertEquals(0.995, sig.getSiglr(), 0.001);
        assertEquals(1.655, sig.getSigdiff(), 0.001);
        assertEquals("3,4", sig.getSigseg());
    }

    @Test
    public void checkBPTestACOX1() throws Exception {

        List<Sample> samples = getSampleACOX1();
        Sig sig = new Lr2gene().checkBP(samples);
        assertEquals("-1.0", sig.getSigStr());
        assertEquals("", sig.getBp());
        assertEquals("", sig.getCn());
        assertEquals(0, sig.getBpi());
        assertEquals(11, sig.getTotal());
        assertEquals(0.0, sig.getSiglr(), 0.001);
        assertEquals(0.0, sig.getSigdiff(), 0.001);
        assertEquals("", sig.getSigseg());
    }

    @Test
    public void checkBPTestALDH9A() throws Exception {

        List<Sample> samples = getSampleALDH9A1();
        Sig sig = new Lr2gene().checkBP(samples);
        assertEquals("-1.0", sig.getSigStr());
        assertEquals("", sig.getBp());
        assertEquals("", sig.getCn());
        assertEquals(0, sig.getBpi());
        assertEquals(10, sig.getTotal());
        assertEquals(0.0, sig.getSiglr(), 0.001);
        assertEquals(0.0, sig.getSigdiff(), 0.001);
        assertEquals("", sig.getSigseg());
    }

    @Test
    public void checkBPTestANKRD44() throws Exception {

        List<Sample> samples = getSampleANKRD44();
        Sig sig = new Lr2gene().checkBP(samples);
//        assertEquals("1.62e-12", sig.getSigStr());
        assertEquals("BP", sig.getBp());
        assertEquals("Del", sig.getCn());
        assertEquals(3, sig.getBpi());
        assertEquals(30, sig.getTotal());
        assertEquals(-4.99, sig.getSiglr(), 0.01);
        assertEquals(4.6318, sig.getSigdiff(), 0.001);
        assertEquals("1,2,3", sig.getSigseg());
    }

    @Test
    public void checkBPTestAPEH() throws Exception {

        List<Sample> samples = getSampleAPEH();
        Sig sig = new Lr2gene().checkBP(samples);
        assertEquals("10", sig.getSigStr());
        assertEquals("BP", sig.getBp());
        assertEquals("Del", sig.getCn());
        assertEquals(1, sig.getBpi());
        assertEquals(4, sig.getTotal());
        assertEquals(-1.98, sig.getSiglr(), 0.01);
        assertEquals(1.26, sig.getSigdiff(), 0.001);
        assertEquals("3", sig.getSigseg());
    }

    @Test
    public void checkBPTestCASC4() throws Exception {

        List<Sample> samples = getSampleCASC4();
        Sig sig = new Lr2gene().checkBP(samples);
        assertEquals("0.000716", sig.getSigStr());
        assertEquals("BP", sig.getBp());
        assertEquals("Amp", sig.getCn());
        assertEquals(3, sig.getBpi());
        assertEquals(10, sig.getTotal());
        assertEquals(1.21, sig.getSiglr(), 0.01);
        assertEquals(1.481, sig.getSigdiff(), 0.001);
        assertEquals("6,7,8", sig.getSigseg());
    }

    @Ignore
    @Test
    public void checkBPTestCPSF4L() throws Exception {

        List<Sample> samples = getSampleCPSF4L();
        Sig sig = new Lr2gene().checkBP(samples);
        assertEquals("-1.0", sig.getSigStr());
        assertEquals("", sig.getBp());
        assertEquals("", sig.getCn());
        assertEquals(0, sig.getBpi());
        assertEquals(7, sig.getTotal());
        assertEquals(0.0, sig.getSiglr(), 0.01);
        assertEquals(0.0, sig.getSigdiff(), 0.001);
        assertEquals("", sig.getSigseg());
    }

    @Test
    public void checkBPTestSLC38A11() throws Exception {

        List<Sample> samples = getSampleSLC38A11();
        Sig sig = new Lr2gene().checkBP(samples);
        assertEquals("-1.0", sig.getSigStr());
        assertEquals("", sig.getBp());
        assertEquals("", sig.getCn());
        assertEquals(0, sig.getBpi());
        assertEquals(12, sig.getTotal());
        assertEquals(0.0, sig.getSiglr(), 0.01);
        assertEquals(0.0, sig.getSigdiff(), 0.001);
        assertEquals("", sig.getSigseg());
    }

    @Test
    public void checkBPTestFIG4() throws Exception {

        List<Sample> samples = getSampleFIG4();
        Sig sig = new Lr2gene().checkBP(samples);
        assertEquals("0.000205", sig.getSigStr());
        assertEquals("BP", sig.getBp());
        assertEquals("Amp", sig.getCn());
        assertEquals(4, sig.getBpi());
        assertEquals(23, sig.getTotal());
        assertEquals(1.21, sig.getSiglr(), 0.01);
        assertEquals(1.3075, sig.getSigdiff(), 0.001);
        assertEquals("20,21,22,23", sig.getSigseg());
    }


    private List<Sample> getSampleNotAGene4865() {

        String sampleName = "HG00096";
        String gene = "not_a_gene_4865";
        String chr = "16";

        List<Sample> samples = new ArrayList<>();
        SamplesTestFactory factory = new SamplesTestFactory(sampleName, gene, chr);
        samples.add(factory.getSample(71963329, 71963668, -6.32));
        samples.add(factory.getSample(71963802, 71964197, -4.88));
        samples.add(factory.getSample(71967163, 71967526, -5.49));
        samples.add(factory.getSample(71967724, 71968138, -9.06));
        samples.add(factory.getSample(71969130, 71969520, -9.06));
        samples.add(factory.getSample(71971064, 71971437, -7.16));
        samples.add(factory.getSample(71976293, 71977077, -7.58));
        samples.add(factory.getSample(71977823, 71978210, -7.97));
        samples.add(factory.getSample(71981245, 71981651, -7.42));
        samples.add(factory.getSample(71983527, 71984300, -7.58));
        samples.add(factory.getSample(71985381, 71985753, -6.43));
        samples.add(factory.getSample(71986683, 71987097, -6.84));
        samples.add(factory.getSample(71987892, 71988309, -7.94));
        samples.add(factory.getSample(71996980, 71997383, -5.54));
        samples.add(factory.getSample(72000908, 72001346, -5.49));
        samples.add(factory.getSample(72001644, 72001966, -9.06));
        samples.add(factory.getSample(72003668, 72004123, -7.43));
        samples.add(factory.getSample(72007098, 72007544, -7.09));
        samples.add(factory.getSample(72011016, 72011481, -8.65));
        samples.add(factory.getSample(72012050, 72012408, -7.41));
        samples.add(factory.getSample(72013529, 72014157, -7.40));
        samples.add(factory.getSample(72015849, 72016218, -5.45));
        samples.add(factory.getSample(72017799, 72018180, -4.15));
        samples.add(factory.getSample(72019997, 72020448, -7.91));
        samples.add(factory.getSample(72024033, 72024335, -5.08));
        samples.add(factory.getSample(72026993, 72027378, -7.21));
        samples.add(factory.getSample(72032010, 72032402, -6.96));
        samples.add(factory.getSample(72033468, 72033951, -9.06));
        return samples;
    }

    private List<Sample> getSamplePKN2() {

        String sampleName = "HG00096";
        String gene = "PKN2";
        String chr = "1";
        List<Sample> samples = new ArrayList<>();
        SamplesTestFactory factory = new SamplesTestFactory(sampleName, gene, chr);

        samples.add(factory.getSample(89149744, 89150413, -2.15));
        samples.add(factory.getSample(89206525, 89207087, 0.43));
        samples.add(factory.getSample(89225761, 89226129, 0.40));
        samples.add(factory.getSample(89235912, 89236276, -0.07));
        samples.add(factory.getSample(89237023, 89237636, -0.35));
        samples.add(factory.getSample(89250208, 89250579, 0.36));
        samples.add(factory.getSample(89251641, 89252008, 0.15));
        samples.add(factory.getSample(89270011, 89272037, -0.42));
        samples.add(factory.getSample(89272852, 89273566, 0.01));
        samples.add(factory.getSample(89279185, 89279534, 1.05));
        samples.add(factory.getSample(89287463, 89287813, 0.08));
        samples.add(factory.getSample(89289842, 89290216, -0.09));
        samples.add(factory.getSample(89294053, 89294415, 0.91));
        samples.add(factory.getSample(89294531, 89294900, -1.24));
        samples.add(factory.getSample(89298286, 89301988, 1.36));
        return samples;
    }

    private List<Sample> getSampleCOMMD1() {
        String sampleName = "HG00097.mapped.ILLUMINA.bwa.GBR.exome.20130415";
        String gene = "COMMD1";
        String chr = "1";
        List<Sample> samples = new ArrayList<>();
        SamplesTestFactory factory = new SamplesTestFactory(sampleName, gene, chr);

        samples.add(factory.getSample(62132652, 62133086, -0.03));
        samples.add(factory.getSample(62227686, 62228277, -0.45));
        samples.add(factory.getSample(62362822, 62363432, -0.77));
        samples.add(factory.getSample(62373383, 62374550, -4.23));
        return samples;
    }

    private List<Sample> getSampleNotAGene4998() {
        String sampleName = "HG00096";
        String gene = "not_a_gene_4998";
        String chr = "17";
        SamplesTestFactory factory = new SamplesTestFactory(sampleName, gene, chr);
        List<Sample> samples = new ArrayList<>();
        samples.add(factory.getSample(28935361, 28935661, 1.54));
        samples.add(factory.getSample(28937857, 28938010, 1.84));
        samples.add(factory.getSample(28942831, 28943091, 1.03));
        samples.add(factory.getSample(28944673, 28944947, 1.58));
        samples.add(factory.getSample(28951587, 28951908, -0.19));
        samples.add(factory.getSample(28951969, 28952605, 0.47));
        samples.add(factory.getSample(28956535, 28956895, 1.72));
        samples.add(factory.getSample(28958625, 28958625, 1.64));
        samples.add(factory.getSample(28959984, 28960491, -0.78));
        samples.add(factory.getSample(28960885, 28964534, -1.13));
        samples.add(factory.getSample(28998930, 28999243, -0.13));
        samples.add(factory.getSample(29058674, 29059244, -5.12));
        samples.add(factory.getSample(29061781, 29062209, -1.40));
        samples.add(factory.getSample(29070195, 29070503, -1.76));
        samples.add(factory.getSample(29086382, 29086666, -1.29));
        samples.add(factory.getSample(29093456, 29093709, -0.91));
        samples.add(factory.getSample(29095756, 29096186, -0.19));
        return samples;
    }

    private List<Sample> getSampleNotAGene6082() {
        String sampleName = "HG000100";
        String gene = "not_a_gene_6082";
        String chr = "X";
        SamplesTestFactory factory = new SamplesTestFactory(sampleName, gene, chr);
        List<Sample> samples = new ArrayList<>();

        samples.add(factory.getSample(101033151, 101033783, 0.99));
        samples.add(factory.getSample(101033802, 101034165, -0.90));
        samples.add(factory.getSample(101034176, 101034562, -0.87));
        samples.add(factory.getSample(101034577, 101034935, -0.50));
        samples.add(factory.getSample(101034957, 101035657, 0.00));

        return samples;
    }

    private List<Sample> getSampleAACS() {
        String sampleName = "HG00096";
        String gene = "AACS";
        String chr = "12";
        List<Sample> samples = new ArrayList<>();
        SamplesTestFactory factory = new SamplesTestFactory(sampleName, gene, chr);

        samples.add(factory.getSample(125549783, 125550385, -5.29));
        samples.add(factory.getSample(125558270, 125558647, -0.66));
        samples.add(factory.getSample(125560884, 125561297, 0.31));
        samples.add(factory.getSample(125570732, 125571074, -0.67));
        samples.add(factory.getSample(125575823, 125576141, -1.98));

        samples.add(factory.getSample(125585382, 125587733, 0.61));
        samples.add(factory.getSample(125591525, 125591919, -0.29));
        samples.add(factory.getSample(125598844, 125599204, -0.13));
        samples.add(factory.getSample(125602056, 125603418, -1.47));
        samples.add(factory.getSample(125604354, 125609714,	-0.16));

        samples.add(factory.getSample(125610221, 125612893, -0.61));
        samples.add(factory.getSample(125613013, 125614168, 0.44));
        samples.add(factory.getSample(125618383, 125618680, -1.39));
        samples.add(factory.getSample(125619165, 125619498, -0.87));
        samples.add(factory.getSample(125621079, 125621488, -0.84));
        samples.add(factory.getSample(125626476, 125627943, 0.37));

        return samples;
    }

    private List<Sample> getSampleAASDHPPT() {
        String sampleName = "HG00096";
        String gene = "AASDHPPT";
        String chr = "11";
        List<Sample> samples = new ArrayList<>();
        SamplesTestFactory factory = new SamplesTestFactory(sampleName, gene, chr);

        samples.add(factory.getSample(105947989, 105948719, -1.52));
        samples.add(factory.getSample(105950077, 105950550, -0.24));
        samples.add(factory.getSample(105961208, 105961538, 1.14));
        samples.add(factory.getSample(105961953, 105962321, 0.85));
        samples.add(factory.getSample(105965188, 105965517, -0.66));
        samples.add(factory.getSample(105967332, 105969492, -0.22));
        return samples;
    }

    private List<Sample> getSampleACOX1() {
        String sampleName = "HG00096";
        String gene = "ACOX1";
        String chr = "17";
        List<Sample> samples = new ArrayList<>();
        SamplesTestFactory factory = new SamplesTestFactory(sampleName, gene, chr);

        samples.add(factory.getSample(73937538, 73942992, -0.67));
        samples.add(factory.getSample(73944228, 73944620, 0.65));
        samples.add(factory.getSample(73945179, 73946063, 0.17));
        samples.add(factory.getSample(73946754, 73947744, -0.06));
        samples.add(factory.getSample(73949280, 73949280, -0.48));

        samples.add(factory.getSample(73951524, 73952115, -0.33));
        samples.add(factory.getSample(73952193, 73952550, 2.86));
        samples.add(factory.getSample(73953402, 73953746, -0.04));
        samples.add(factory.getSample(73956144, 73956600, -0.31));
        samples.add(factory.getSample(73968953, 73969989, -0.08));
        samples.add(factory.getSample(73974470, 73975642, -0.16));

        return samples;
    }

    private List<Sample> getSampleALDH9A1() {
        String sampleName = "HG00096";
        String gene = "ALDH9A1";
        String chr = "1";
        List<Sample> samples = new ArrayList<>();
        SamplesTestFactory factory = new SamplesTestFactory(sampleName, gene, chr);

        samples.add(factory.getSample(165631386, 165632501, 0.01));
        samples.add(factory.getSample(165634123, 165634508, -0.12));
        samples.add(factory.getSample(165636427, 165636763, 0.03));
        samples.add(factory.getSample(165638016, 165638016, -0.09));
        samples.add(factory.getSample(165648539, 165648946, -0.31));

        samples.add(factory.getSample(165649661, 165650019, 0.10));
        samples.add(factory.getSample(165651189, 165651587, 0.93));
        samples.add(factory.getSample(165652089, 165652490, 0.45));
        samples.add(factory.getSample(165664342, 165664751, 1.32));
        samples.add(factory.getSample(165667437, 165668182, -0.19));

        return samples;
    }

    private List<Sample> getSampleANKRD44() {
        String sampleName = "HG00096";
        String gene = "ANKRD44";
        String chr = "2";
        List<Sample> samples = new ArrayList<>();
        SamplesTestFactory factory = new SamplesTestFactory(sampleName, gene, chr);

        samples.add(factory.getSample(197831691,	197832348,	-4.81));
        samples.add(factory.getSample(197851336,	197855701,	-5.14));
        samples.add(factory.getSample(197858152,	197858544,	-5.03));
        samples.add(factory.getSample(197859154,	197860303,	-0.48));
        samples.add(factory.getSample(197862900,	197863283,	-1.12));

        samples.add(factory.getSample(197863511,	197863907,	-0.62));
        samples.add(factory.getSample(197865005,	197865346,	-0.40));
        samples.add(factory.getSample(197866327,	197866625,	1.35));
        samples.add(factory.getSample(197870346,	197870772,	-0.28));
        samples.add(factory.getSample(197872410,	197872776,	0.60));

        samples.add(factory.getSample(197873517,	197873874,	1.47));
        samples.add(factory.getSample(197878094,	197878491,	-0.14));
        samples.add(factory.getSample(197889782,	197890110,	-0.31));
        samples.add(factory.getSample(197942468,	197943615,	-0.84));
        samples.add(factory.getSample(197946198,	197946535,	0.72));

        samples.add(factory.getSample(197947981,	197948351,	-0.01));
        samples.add(factory.getSample(197951228,	197951594,	-0.45));
        samples.add(factory.getSample(197953302,	197953670,	1.01));
        samples.add(factory.getSample(197954538,	197954903,	-0.43));
        samples.add(factory.getSample(197964091,	197964792,	-0.82));

        samples.add(factory.getSample(197975325,	197975704,	0.09));
        samples.add(factory.getSample(197985946,	197986323,	-1.50));
        samples.add(factory.getSample(197987225,	197987615,	-0.09));
        samples.add(factory.getSample(197989959,	197990349,	-0.27));
        samples.add(factory.getSample(197990440,	197990885,	-0.57));

        samples.add(factory.getSample(198001149,	198001459,	0.06));
        samples.add(factory.getSample(198011622,	198011982,	-0.98));
        samples.add(factory.getSample(198051563,	198051940,	1.49));
        samples.add(factory.getSample(198175135,	198175571,	-3.22));
        samples.add(factory.getSample(198175799,	198175947,	-4.02));

        return samples;
    }

    private List<Sample> getSampleCASC4() {
        String sampleName = "HG00096";
        String gene = "CASC4";
        String chr = "15";
        List<Sample> samples = new ArrayList<>();
        SamplesTestFactory factory = new SamplesTestFactory(sampleName, gene, chr);

        samples.add(factory.getSample(44580830,	44581687,	-0.00));
        samples.add(factory.getSample(44614998,	44615390,	-0.49));
        samples.add(factory.getSample(44620753,	44621128,	0.14));
        samples.add(factory.getSample(44624059,	44630225,	0.20));
        samples.add(factory.getSample(44629823,	44630225,	-0.43));

        samples.add(factory.getSample(44630283,	44630675,	1.05));
        samples.add(factory.getSample(44671785,	44672131,	1.70));
        samples.add(factory.getSample(44672889,	44673292,	0.88));
        samples.add(factory.getSample(44694937,	44695417,	-0.45));
        samples.add(factory.getSample(44705351,	44708057,	-0.87));

        return samples;
    }

    private List<Sample> getSampleAPEH() {
        String sampleName = "HG00096";
        String gene = "APEH";
        String chr = "3";
        List<Sample> samples = new ArrayList<>();
        SamplesTestFactory factory = new SamplesTestFactory(sampleName, gene, chr);

        samples.add(factory.getSample(49711385,	49711385,	-0.93));
        samples.add(factory.getSample(49716130,	49716475,	-0.20));
        samples.add(factory.getSample(49716827,	49716827,	-1.98));
        samples.add(factory.getSample(49717793,	49727305,	-1.03));

        return samples;
    }

    private List<Sample> getSampleCPSF4L() {
        String sampleName = "HG00096";
        String gene = "CPSF4L";
        String chr = "17";
        List<Sample> samples = new ArrayList<>();
        SamplesTestFactory factory = new SamplesTestFactory(sampleName, gene, chr);

        samples.add(factory.getSample(71245882,	71246557,	-7.74));
        samples.add(factory.getSample(71246959,	71247400,	-3.71));
        samples.add(factory.getSample(71248645,	71249041,	-5.81));
        samples.add(factory.getSample(71249928,	71250336,	-4.06));
        samples.add(factory.getSample(71253695,	71254104,	-5.63));
        samples.add(factory.getSample(71256905,	71257300,	-3.86));
        samples.add(factory.getSample(71257758,	71258112,	-5.91));

        return samples;
    }

    private List<Sample> getSampleSLC38A11() {
        String sampleName = "HG00099";
        String gene = "SLC38A11";
        String chr = "2";
        List<Sample> samples = new ArrayList<>();
        SamplesTestFactory factory = new SamplesTestFactory(sampleName, gene, chr);

        samples.add(factory.getSample(165752601,	165755363,	-0.27));
        samples.add(factory.getSample(165765051,	165765416,	-0.52));
        samples.add(factory.getSample(165767977,	165768381,	-0.05));
        samples.add(factory.getSample(165771522,	165771862,	-0.29));
        samples.add(factory.getSample(165772265,	165772562,	-1.86));

        samples.add(factory.getSample(165793721,	165794040,	-0.77));
        samples.add(factory.getSample(165795833,	165796127,	-0.06));
        samples.add(factory.getSample(165800923,	165801336,	0.54));
        samples.add(factory.getSample(165801958,	165802320,	-0.36));
        samples.add(factory.getSample(165809123,	165809423,	-0.43));

        samples.add(factory.getSample(165811021,	165811316,	0.53));
        samples.add(factory.getSample(165811526,	165812104,	0.74));

        return samples;
    }

    private List<Sample> getSampleFIG4() {
        String sampleName = "HG000100";
        String gene = "FIG4";
        String chr = "6";
        List<Sample> samples = new ArrayList<>();
        SamplesTestFactory factory = new SamplesTestFactory(sampleName, gene, chr);

        samples.add(factory.getSample(110012300,	110012815,	-1.18));
        samples.add(factory.getSample(110022759,	110023069,	-4.34));
        samples.add(factory.getSample(110036172,	110036513,	-0.64));
        samples.add(factory.getSample(110037498,	110037892,	0.42));
        samples.add(factory.getSample(110048168,	110048547,	0.21));

        samples.add(factory.getSample(110052779,	110054059,	1.04));
        samples.add(factory.getSample(110056239,	110056614,	0.75));
        samples.add(factory.getSample(110059380,	110059738,	0.73));
        samples.add(factory.getSample(110062519,	110062972,	0.94));
        samples.add(factory.getSample(110064216,	110065081,	-0.51));

        samples.add(factory.getSample(110081298,	110081718,	-0.10));
        samples.add(factory.getSample(110083135,	110083521,	0.04));
        samples.add(factory.getSample(110085002,	110085309,	-0.15));
        samples.add(factory.getSample(110086078,	110086433,	-0.52));
        samples.add(factory.getSample(110087813,	110088197,	0.86));

        samples.add(factory.getSample(110098016,	110098337,	0.61));
        samples.add(factory.getSample(110106041,	110106374,	0.54));
        samples.add(factory.getSample(110106746,	110107785,	-0.28));
        samples.add(factory.getSample(110110667,	110110976,	-0.32));
        samples.add(factory.getSample(110112450,	110112839,	1.32));

        samples.add(factory.getSample(110113639,	110113945,	1.40));
        samples.add(factory.getSample(110117839,	110118151,	1.04));
        samples.add(factory.getSample(110146162,	110146720,	1.07));

        return samples;
    }
 }