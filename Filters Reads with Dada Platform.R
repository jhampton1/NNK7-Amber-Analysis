library(dada2)
#create path to unzipped NGS fastq files
path <- "/Users/traehampton/Desktop/NGS FastQ"
NNK7F <- "/Users/traehampton/Desktop/NGS FastQ/1_S1_L001_R1_001.fastq"
NNK7R <- "/Users/traehampton/Desktop/NGS FastQ/1_S1_L001_R2_001.fastq"
plotQualityProfile(NNK7F)
plotQualityProfile(NNK7R)
#creates directory for filtered files
filtFs <- file.path(path, "filtered", paste0(NNK7F, "_F2_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(NNK7R, "_R2_filt.fastq"))

#this filters the reads, set trim length corresponding to poor sequence quality from plots
out <- filterAndTrim(NNK7F, filtFs, NNK7R, filtRs, truncLen=c(120,95),
                     maxN=1, maxEE=c(1,1), rm.phix=FALSE,
                     compress=FALSE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)