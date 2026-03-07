#!/bin/bash 

## Running CPAT on output by ISAR 

annotations=/home/neuro/Documents/NMD_analysis/Analysis/Genome/cpat
fasta=/home/neuro/Documents/NMD_analysis/Analysis/NMD-analysis/output/ISAR

cpat.py -x ${annotations}/Mouse_Hexamer.tsv -d ${annotations}/Mm.39.logit.RData -g ${fasta}/isoformSwitchAnalyzeR_isoform_nt.fasta -o ./