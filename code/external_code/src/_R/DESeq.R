library(DESeq2)
library(gplots)
library(tibble)


jobid="NMD_analysis"
metadata="/home/neuro/Documents/NMD_analysis/Analysis/NMD-analysis/data/data2.txt"
refdir="/home/neuro/Documents/NMD_analysis/Analysis/Results/EISA"
stringency=0.99
fitMode="linear"

args <- commandArgs(trailingOnly = TRUE)

jobID <- jobid
tablePath <- metadata
inputFolder <- refdir

path= getwd()

outputFolder <- file.path(path, "tmp")

sampleTable <- read.csv(tablePath, sep=",")


dir.create(outputFolder)

for( readType in c("exonic","intronic") )
{
	print( paste("Analyzing ", readType, " reads ...", sep="") )

	for( batch in 1:max(sampleTable$Batch) )
	{
		print( paste("Analyzing batch ", batch, sep="") )

		cds <- DESeqDataSetFromHTSeqCount( # omits special rows from htseq-count
		    sampleTable = sampleTable[sampleTable$Batch==batch & sampleTable$ReadType==readType,],
		    directory = inputFolder,
		    design = ~ 1 # required parameter
		)

		# estimate size factors, dispersions, and perform VST
		vsd <- vst(object = cds, blind = T)

        normalized <- assay(vsd)
		normalized <- normalized[apply(normalized,1,sd)>0,] # remove genes that have zero reads all across the samples

		write.table(
			normalized,
			paste(outputFolder,"/vsd_normalized.",readType,".batch",batch,".mx.txt",sep=""),
			quote=F, sep="\t")


		if( batch == 1 )
			all_normalized <- normalized
		else
			all_normalized <- cbind( all_normalized,
			   normalized[ match( rownames(all_normalized), rownames(normalized) ), ] )

		# transform the normalized values, by first centering the columns (mean=0), and then subtracting row medians
		centered <- t(apply(scale(normalized,scale=F), 1, function(y) y - median(y,na.rm=T) ))
		write.table(
			centered,
			paste(outputFolder,"/vsd_normalized.",readType,".batch",batch,".centered.mx.txt",sep=""),
			quote=F, sep="\t" )

		if( batch == 1 )
			all_centered <- centered
		else
			all_centered <- cbind( all_centered,
			   centered[ match( rownames(all_centered), rownames(centered) ), ] )
	}

	all_normalized <- all_normalized[ apply(all_normalized,1,function(x) (sum(is.na(x))==0) ) , ]
	write.table(
		rownames_to_column(as.data.frame(all_normalized), var = 'GeneID'),
		paste(outputFolder,"/vsd_normalized.",readType,".all.mx.txt",sep=""),
		quote=F, sep="\t", append=FALSE )

	all_centered <- all_centered[ apply(all_centered,1,function(x) (sum(is.na(x))==0) ) , ]
	write.table(
		rownames_to_column(as.data.frame(all_centered), var = 'GeneID'),
		paste(outputFolder,"/vsd_normalized.",readType,".all.centered.mx.txt",sep=""),
		quote=F, sep="\t", append=FALSE )

	mat = 1 - cor( all_centered, use="na.or.complete", method="pearson" )
	width  <- ncol(mat) * 36.9 + 231
	height <- nrow(mat) * 36.9 + 231
	jpeg(file=paste(outputFolder,"/vsd_normalized.",readType,".all.centered.correl.heatmap.jpg",sep=""),
		width=width,height=height)
	heatmap.2(mat, trace="none", margin=c(13, 13),breaks=seq(0,2,length.out=256),col=colorRampPalette(c("red","white","blue"))(255), key.title=NA, key.xlab="Pearson distance", density.info="none", key.ylab=NA)
}



