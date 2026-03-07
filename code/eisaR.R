## Performing exon-intron split analysis
## This analysis allows separately quantifying exonic and intronic alignments in RNA-seq data, in order to measure changes
## in mature RNA and pre-mRNA reads


## First step is to
## 1) Preare the annotations (exonic and gene body coordinates)
## 2) quantifying RNA-seq alignments in exons and introns

source("libraries.R")
library(eisaR)
library(AnnotationDbi)
library(EnsDb.Mmusculus.v79)
library(org.Mm.eg.db)
library(QuasR)
library(doParallel)
library(parallel)


parallel::makeCluster(16, setup_strategy = "sequential")

## Preparing the annotation file
### extract filtered exonic and gene body regions

regS <- getRegionsFromTxDb(txdb = EnsDb.Mmusculus.v79, strandedData = TRUE)

##
regU <- getRegionsFromTxDb(txdb = EnsDb.Mmusculus.v79, strandedData = FALSE)


## fastq files

fq = list.files("/home/neuro/Documents/NMD_analysis/Analysis/Results/Trimming", pattern = "1.fq.gz$", full.names = TRUE)
fq2 = list.files("/home/neuro/Documents/NMD_analysis/Analysis/Results/Trimming", pattern = "2.fq.gz$", full.names = TRUE)

fq = data.frame("FileName1" = fq,
                "FileName2" = fq2,
                "SampleName" = c(212:229))
write_delim(fq, "/home/neuro/Documents/NMD_analysis/Analysis/NMD-analysis/data/samples.txt", delim = "\t")

sampleFile <- "/home/neuro/Documents/NMD_analysis/Analysis/NMD-analysis/data/samples.txt"


proj <- qAlign(sampleFile = sampleFile,
               genome = "BSgenome.Mmusculus.UCSC.mm39",
               aligner = "Rhisat2", splicedAlignment = TRUE,
               paired="fr", clObj = cluster)


save(prof, file="/home/neuro/Documents/NMD_analysis/Analysis/NMD-analysis/output/eisaR.Rda")
