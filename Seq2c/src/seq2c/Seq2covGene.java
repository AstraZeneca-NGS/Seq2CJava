package seq2c;

import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author Petr_Rastegaev
 */
public class Seq2covGene {

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
