/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package seq2c;
import java.util.*;
import java.util.Collections;
import org.apache.commons.math3.stat.*;
import org.apache.commons.math3.stat.inference.TTest;
/**
 *
 * @author Petr_Rastegaev
 */
public class Lr2gene {
   private ArrayList<Sample> inputGenes;
   public Lr2gene(Seq2cov_return sq){
       Seq2cov_return inputGenes = sq;
   }

private void Lr2gene_mainloop() {
    HashMap<String,HashMap<String,ArrayList<Sample>>> g2amp = new HashMap();
    boolean opt_c = false;
    for(Sample sqr : inputGenes){
        HashMap<String,ArrayList<Sample>> gq;
        ArrayList<Sample> sq2amparr;
        if(g2amp.containsKey(sqr.getName())){
            gq = g2amp.get(sqr.getName());
            if(gq.containsKey(sqr.getGene())){
                sq2amparr = gq.get(sqr.getGene());
                sq2amparr.add(sqr);
                
            }else{
                 
                 sq2amparr = new ArrayList();
                 sq2amparr.add(sqr);
            }
            
        }else{
            gq = new HashMap<>();
            sq2amparr = new ArrayList();
            sq2amparr.add(sqr);
            
        }
        gq.put(sqr.getGene(),sq2amparr);
        g2amp.put(sqr.getName(), gq);
    }
    
    for(Map.Entry<String,HashMap<String,ArrayList<Sample>>> entry: g2amp.entrySet()){
        String Sample = entry.getKey();
        HashMap<String,ArrayList<Sample>> gq = entry.getValue();
        for(Map.Entry<String,ArrayList<Sample>> entry2: gq.entrySet()){
        String gene = entry2.getKey();
        ArrayList<Sample> sq2amparr = entry2.getValue();
        
        Collections.sort(sq2amparr,new Comparator<Sample>() {
                public int compare(Sample ints, Sample otherInts) {
                    //return ints[3].compareTo(otherInts[3]);
                    return Double.compare(ints.getStart(),otherInts.getStart());
                }
            });
        int i;
        if(opt_c){
            i = 11;
        } else {
            i = 10;
        }
        ArrayList<Double> lr = new ArrayList();
        for(Sample sqarr  : sq2amparr){
            lr.add(sqarr.getNorm3());
        }
        double lr_med;
        if(lr.size()> 1){
           lr_med = 1.0;//median
        } else {
            lr_med = lr.get(0);
        }
        Sig sig = checkBP(sq2amparr);
        if(sig == null){
            if(lr_med > AMP){

            }else if(lr_med <= DEL){

            }
        }
    }    
} 
}
private Sig checkBP(ArrayList<Sample> segs){
    double[][] arr = new double[segs.size()][3];
    double[] lr = new double[segs.size()];
    int i = 0;
    for(Sample sqarr  : segs){ 
        arr[i][0] = sqarr.getStart();
        arr[i][1] = sqarr.getNorm3();
        arr[i][2] = i;
        i++;
    }
    i = 0;
    for(double[] l : arr){
        lr[i] = l[1];
        i++;
    }
    double max = StatUtils.max(lr);
    double min = StatUtils.min(lr);
    double mid = (max + min)/2;
    double[] bps = getBPS(lr);
    double minbp = 1;
    double maxmd = 1;
    for(double bp : bps){
        
        ArrayList<double[]> bm = new ArrayList<>();
        ArrayList<double[]> up = new ArrayList<>();
        ArrayList<Double> lrup = new ArrayList<>();
        ArrayList<Double> lrbm = new ArrayList<>();
        ArrayList<Integer> upsegArr = new ArrayList<>();
        ArrayList<Integer> bmsegArr = new ArrayList<>();
        i =0;
        for(double[]intarr : arr){
            if(intarr[1] > bp) {
                    up.add(intarr);
                    lrup.add(intarr[1]);
                    upsegArr.add(i+1);
                } else{ 
                    bm.add(intarr);
                    lrbm.add(intarr[1]);
                    bmsegArr.add(i+1);
                }
            i++;
        }
        double lrupm = StatUtils.percentile(convertDoubles(lrup),50);
        double lrbmm = StatUtils.percentile(convertDoubles(lrbm),50);
        String cn = lrbmm < -0.35 ? "Del" : (lrupm > 0.35 && Math.abs(lrbmm) < Math.abs(lrupm) ? "Amp" : "NA");
        if("NA".equals(cn)) continue;
        double[] bmiscArr = isConsecutive(bm);
        double[] upiscArr = isConsecutive(up);
        if(bmiscArr[0] != 0){
            
        } else if(upiscArr[0] != 0){
            
        }
    }
    return sig;    
} 

private double[] isConsecutive(ArrayList<double[]> ref) {
    
    double skip = 0;
    double si = -1;
    double sii = -1;
    double[] result = new double[3];
    for(int i = 0; i < ref.size(); i++){
        skip += ref.get(i)[2] - ref.get(i - 1)[2] - 1;
        if(ref.get(i)[2] - ref.get(i - 1)[2] == 2){
            si = ref.get(i)[2] -1;
            sii = i;
        }    
    }
    if((skip == 0)||(skip == 1 && ref.size() > 10)) {
        result[0] = 1;  
    } else {
        result[0] = 0;
    }
    result[1] = si;
    result[2] = sii;
    return result;
}

private double[] getBPS(double[] lr){
    ArrayList<Double> bpsArr = new ArrayList<>();
    double[][] dis = new double[lr.length][3];
    java.util.Arrays.sort(dis);
    
    for(int i =0; i < lr.length; i++){
        dis[i][0] = lr[i] - lr[i-1];
        dis[i][1] = lr[i];
        dis[i][2] = lr[i -1];
    }
    
    java.util.Arrays.sort(dis, new java.util.Comparator<double[]>() {
        public int compare(double[] a, double[] b) {
        return Double.compare(a[0], b[0]);
        }});
    
    for(double[] bp : dis){
        if(bp[0]< 0.1) break;
        bpsArr.add((bp[1]+bp[2])/2);
        
    }
    return convertDoubles(bpsArr);
}

private Sig findBP(double[] lr){
    if (lr.length < 15) return null;
    double minp = 1;
    double bpi = 0;
    double siglr = 0;
    String cn = "NA";
    double mindiff = 0;
    String sigseg = "";
    
    
}

private void getCalls(ArrayList<double[]> bm, ArrayList<double[]> up){
  double[] tlr1 = new double[bm.size()];
  double[] ti1 = new double[bm.size()];
  double[] tlr2 = new double[up.size()];
  double[] ti2 = new double[up.size()];
  double[] result = new double[2];
  int i = 0;
  for(double[] b: bm){
    tlr1[i] = b[1];
    ti1[i] = b[2];
    i++;
  }
  i = 0;
  for(double[] u: up){
    tlr2[i] = u[1];
    ti2[i] = u[2];
    i++;
  }
  i = 0;    
    
}

private double[] isSig(ArrayList<double[]> bm, ArrayList<double[]> up ){
  double[] bm_x = new double[bm.size()];
  double[] up_y = new double[up.size()];
  double[] result = new double[2];
  int i = 0;
  for(double[] b: bm){
    bm_x[i] = b[1];
    i++;
  }
  i = 0;
  for(double[] u: up){
    up_y[i] = u[1];
    i++;
  }
  i = 0;
  if(bm_x.length >=3 && up_y.length >= 3){
      TTest ttest = new TTest();
      double p = ttest.tTest(bm_x,up_y);
      double diff = StatUtils.mean(bm_x) - StatUtils.mean(up_y);
      if((p < PVALUE && Math.abs(diff) >= MINDIFF)||(p < 0.001 && Math.abs(diff) >= MINDIFF &&(Math.abs(StatUtils.mean(bm_x)) > 0.8 ||Math.abs(StatUtils.mean(up_y))>0.8))){
        result[0] = p;
        result[1] = Math.abs(diff);
      } 
  }else if(bm_x.length >= MINSEGS && up_y.length >= 3){
        double med = StatUtils.percentile(up_y, 50);
        double[] d = new double[up_y.length];
        i =0;
        for(double m: up_y){
            d[i] = Math.abs(m - med);  
            i++;
        }
        double mad = StatUtils.percentile(d, 50);
        if(mad ==0) mad += 0.1;
        double[] t = new double[bm_x.length];
        i=0;
        for(double m: bm_x){
            t[i] = (m - med)/mad;  
            i++;
        }
        double mean = StatUtils.mean(t);
        double sum = StatUtils.sum(t);
        double diff = Math.abs(StatUtils.mean(bm_x) - StatUtils.mean(up_y));
        if(Math.abs(sum) > MINMAD && diff > EXONDIFF) {
            result[0] = Math.abs(mean);
            result[1] = diff;
        }
        
      } else if(up_y.length >= MINSEGS && bm_x.length >= 3){
          double med = StatUtils.percentile(up_y, 50);
        double[] d = new double[bm_x.length];
        i =0;
        for(double m: bm_x){
            d[i] = Math.abs(m - med);  
            i++;
        }
        double mad = StatUtils.percentile(d, 50);
        if(mad ==0) mad += 0.1;
        double[] t = new double[up_y.length];
        i=0;
        for(double m: up_y){
            t[i] = (m - med)/mad;  
            i++;
        }
        double mean = StatUtils.mean(t);
        double sum = StatUtils.sum(t);
        double diff = Math.abs(StatUtils.mean(bm_x) - StatUtils.mean(up_y));
        if(Math.abs(sum) > MINMAD && diff > EXONDIFF) {
            result[0] = Math.abs(mean);
            result[1] = diff;
        }
        
      }else{
          result[0] = -1;
          result[1] = 0;
      }
    return result;
}

private double[] convertDoubles(List<Double> doubles)
{
    double[] ret = new double[doubles.size()];
    Iterator<Double> iterator = doubles.iterator();
    int i=0;
    while(iterator.hasNext())
    {
        ret[i] = iterator.next().doubleValue();
        i++;
    }
    return ret;
}

}