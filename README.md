# Seq2C Java
## Introduction

Seq2CJava is a CNV caller written in Java. It is a Java port of [Seq2C Perl Script](https://github.com/AstraZeneca-NGS/Seq2C).

The original Perl Seq2C is a CNV caller, a tool for Copy-number Vaiation (CNV) detection.

Seq2C normalizes the coverage from targeted sequencing to CNV log2 ratio, and detects CNVs in regions that have abnormally high or low read depths compared to the rest. 

As input, Seq2CJava takes list of BAM files to process in txt format, a BED file with target regions and, optional, list of control samples names.

The output columns contain a description and statistical info for the input samples. See section Output Columns for list of columns in the output.

## Quick start
The following is an example command to run Seq2c:

'Seq2C sample2bam.txt regions.bed coverage.txt [sample_name1[:sample_name2]]'

where:

    1. sample2bam.txt   -   file of samples and bam files, each sample at new line
    2. regions.bed      -   bed file with regions of interest, with at least 4 columns
    3. coverage.txt     -   file for coverage output
    4. sample_name      -   optional control sample names. For multiple controls, separate them using :

## Program Options

    *  `-i`
        Manages the number of threads that do the work
        If this parameter is missing, then the mode is one-thread
        If you add the `-i` parameter, the number of threads equals to the number of processor cores
        The parameter '-i threads' sets the number of threads explicitly
    *   `-r`
        Manages the separate launch of two parts of code
        For launching the first  part of code option '–r 1' (analogue to seq2cov.pl script in Perl) is used
        For launching the second part of code option '–r 2' (analogue to bam2reads.pl, cov2lr.pl and lr2gene.pl scripts in Perl)
        Without option '–r' all code will run from start to end
    *   '-M' float
        When considering partial deletions less than 3 exons/amplicons, the minimum MAD value, in addition to -d,
        before considering it to be amplified or deleted
        Default: 10
    *   '-c'
        Indidate that control sample is used for normalization
    *   '-d' float
        When considering >=3 exons deletion/amplification within a gene, the minimum differences between the log2 of two segments
        Default: 0.7
    *   '-p' float (0-1)
        The p-value for t-test when the breakpoint is in the middle with min exons/amplicons >= [-e]
        Default: 0.000001
    *   '-A' float
        Minimum log2 ratio for a whole gene to be considered amplified
        Default: 1.50
    *   '-D' float
        Minimum log2 ratio for a whole gene to be considered deleted
        Default: -2.00
    *   '-E' float
        Minimum mean log2 ratio difference for <3 exon deletion/amplification to be called
        Default: 1.25
    *   '-R' float (0-1)
        If a breakpoint has been detected more than "float" fraction of samples, it's considered false positive and removed
        Default: 0.1, or 10%.  Use in combination with -N
    *   '-N'
        If a breakpoint has been detected more than "float" fraction of samples, it's considered false positive and removed
        Default: 0.1, or 10%.  Use in combination with -N
    *   '-t' float
        When considering breakpoint in the middle of a gene, the minimum differences between the log2 of two segments
        Default: 0.7
    *   '-P' float (0-1)
        The p-value for t-test when the breakpoint is in the middle with min exons/amplicons >= [-e]
        Default: 0.000001
    *   '-e' float
        When considering breakpoint in the middle of a gene, the minimum number of exons
        Default: 8



## Output columns
    1.Sample - sample name 
    2.Gene - gene name
    3.Chr - chromosome name
    4.Start - start position of a gene
    5.End - end position of a gene
    6.Length - length of a gene
    7.Log2ratio - CNV log2 ratio of a sample
    8.Sig - significance
    9.BP_Whole - "BP" or "Whole"
    10.Amp_Del - type of CNV: Amplification("Amp") or Deletion("Del")
    11.Ab_seg - affected segments
    12.Total_seg - total number of segments
    13.Ab_log2ratio - CNV log2 ratio
    14.Log2r_Diff - difference between CNV log2 ratio
    15.Ab_Seg_Loc - segment location
    16.Ab_Samples - control samples
    17.Ab_Samples_Pcnt - control samples percent
 

