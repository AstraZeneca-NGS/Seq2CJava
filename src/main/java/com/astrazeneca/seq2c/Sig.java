package com.astrazeneca.seq2c;

import org.apache.commons.math3.util.Precision;

public class Sig {

    private double sig;
    private String bp;
    private double minp;
    private int bpi;
    private double siglr;
    private String cn;
    private double mindiff;
    private String sigseg;
    private double sigdiff;
    private int total;

    public Sig(double sig, double minp, int bpi, String bp, double siglr, String cn, int total, double mindiff, String sigseg, double sigdiff) {
        this.sig = sig;
        this.minp = minp;
        this.bpi = bpi;
        this.bp = bp;
        this.siglr = siglr;
        this.cn = cn;
        this.mindiff = mindiff;
        this.sigseg = sigseg;
        this.sigdiff = sigdiff;
        this.total = total;
    }

    public Sig() {
        this.sig = -1;
    }

    public void addSig(double sig) {
        this.setSig(sig);
    }

    public void addSdiff(double sigdiff) {
        this.setSigdiff(sigdiff);
    }

    /**
     * @return the bp
     */
    public String getBp() {
        return bp;
    }

    /**
     * @param bp
     *            the bp to set
     */
    public void setBp(String bp) {
        this.bp = bp;
    }

    /**
     * @return the minp
     */
    public double getMinp() {
        return minp;
    }

    /**
     * @param minp
     *            the minp to set
     */
    public void setMinp(double minp) {
        this.minp = minp;
    }

    /**
     * @return the bpi
     */
    public int getBpi() {
        return bpi;
    }

    /**
     * @param bpi
     *            the bpi to set
     */
    public void setBpi(int bpi) {
        this.bpi = bpi;
    }

    /**
     * @return the siglr
     */
    public double getSiglr() {
        return siglr;
    }

    /**
     * @param siglr
     *            the siglr to set
     */
    public void setSiglr(double siglr) {
        this.siglr = siglr;
    }

    /**
     * @return the cn
     */
    public String getCn() {
        return cn;
    }

    /**
     * @param cn
     *            the cn to set
     */
    public void setCn(String cn) {
        this.cn = cn;
    }

    /**
     * @return the mindiff
     */
    public double getMindiff() {
        return mindiff;
    }

    /**
     * @param mindiff
     *            the mindiff to set
     */
    public void setMindiff(double mindiff) {
        this.mindiff = mindiff;
    }

    /**
     * @return the sigseg
     */
    public String getSigseg() {
        return sigseg;
    }

    /**
     * @param sigseg
     *            the sigseg to set
     */
    public void setSigseg(String sigseg) {
        this.sigseg = sigseg;
    }

    /**
     * @return the sigdiff
     */
    public double getSigdiff() {
        return sigdiff;
    }

    /**
     * @param sigdiff
     *            the sigdiff to set
     */
    public void setSigdiff(double sigdiff) {
        this.sigdiff = sigdiff;
    }

    /**
     * @return the total
     */
    public int getTotal() {
        return total;
    }

    /**
     * @param total
     *            the total to set
     */
    public void setTotal(int total) {
        this.total = total;
    }

    /**
     * @return the sig
     */
    public double getSig() {
        return sig;
    }

    /**
     * @param sig
     *            the sig to set
     */
    public void setSig(double sig) {
        this.sig = sig;
    }

    public StringBuilder getName() {
        StringBuilder result = new StringBuilder();
        result.append(Precision.round(getSig(), 3)).append("\t");
        result.append(getBp()).append("\t");
        result.append(getCn()).append("\t");
        result.append(getBpi()).append("\t");
        result.append(getTotal()).append("\t");
        result.append(Precision.round(getSiglr(), 3)).append("\t");
        result.append(Precision.round(getSigdiff(), 1)).append("\t");
        result.append(getSigseg());
        return result;
    }

    public void update(Sig tmpSig) {
        this.bpi = tmpSig.getBpi();
        this.bp = tmpSig.getBp();
        this.siglr = tmpSig.getSiglr();
        this.cn = tmpSig.getCn();
        this.total = tmpSig.getTotal();
        this.sigseg = tmpSig.getSigseg();
    }
}
