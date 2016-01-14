# Seq2C Java
## Introduction

Seq2CJava is a CNV caller written in Java. It is a Java port of [Seq2C Perl Script](https://github.com/AstraZeneca-NGS/Seq2C).

The original Perl Seq2C is a CNV caller, a tool for Copy-number Vaiation (CNV) detection.

Seq2C normalizes the coverage from targeted sequencing to CNV log2 ratio, and detects CNVs in regions that have abnormally high or low read depths compared to the rest. 

As input, Seq2CJava takes list of BAM files to process in txt format, a BED file with target regions and, optional, list of control samples names.

The output columns contain a description and statistical info for the input samples. See section Output Columns for list of columns in the output.
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
 
