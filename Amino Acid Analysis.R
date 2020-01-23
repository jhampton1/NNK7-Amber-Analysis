library(microseq)
library(RColorBrewer)
library(plyr)
library(dplyr)
library(stringr)
library(gplots)
NNK7Ffilt <- readFastq("/Users/traehampton/Downloads/18347Wns_N18209/1_S1_L001_R1_001.fastq")
NNK7Rfilt <- readFastq("/Users/traehampton/Documents/Research/CD44/19086Wns_N19031/Round_04_S4_L001_R2_001.fastq")
#Define the following variables
libraryseq <- "GCC.{36}GCG" #change this regex to match specific library
beginning <- 22 #beginning of library in DNA string
end <- 55 #end of library in DNA string
initialcodon <- 27 #initial position of library with spaces between codons
endcodon <- 76 #terminal position of library with spaces between codons
lib <- 12 #number of codons in the library region
del <- 7 #number of codons before library


#slices out matches that contain start followed by 24 bases to reverse primer
NNK7Ffilt21 <- gregexpr(libraryseq,NNK7Ffilt[,2],extract = TRUE)
NNK7Rrevcomp <- reverseComplement(NNK7Rfilt[,2],reverse = TRUE) #gives reverse complement of reverse reads
NNK7Rcompfilt21 <- gregexpr(libraryseq,NNK7Rrevcomp,extract = TRUE)

#this compares the forward and reverse strands, only allowing for one mismatch in the primers, no mismatches allowed in the library region
n <- length(NNK7Ffilt21)
NNK7Fgood <- vector()
for(i in c(1:n)){
  if(NNK7Ffilt21[[i]][1] == NNK7Rcompfilt21[[i]][1]){
    NNK7Fgood[i] <- NNK7Ffilt21[[i]][1]
  }
  else{
    split <- strsplit(c(NNK7Ffilt21[[i]][1],NNK7Rcompfilt21[[i]][1]), split = "")
    diff <- which(split[[1]] != split[[2]])
    if(length(diff) < 2 && length(diff) > 0){
      for(x in c(1:length(diff))){
        if(diff[[x]] < beginning || diff[[x]] > end){
          NNK7Fgood[i] <- NNK7Ffilt21[[i]][1]
        }
        else{
          NNK7Fgood[i] <- ""
        }
      }
    }
    else{
      NNK7Fgood[i] <- ""
    }
  }
}
NNK7Fgood <- as.data.frame(NNK7Fgood)
NNK7Fgood <- NNK7Fgood[!apply(is.na(NNK7Fgood) | NNK7Fgood == "", 1, all),]

#this separates nucleotides into codons
codons <- gsub("(...)", "\\1 \\2", NNK7Ffilt21)

#this creates dataframe of sequences with reads organized by frequency
seqcount <- as.data.frame(sort(table(codons), decreasing = TRUE))

#this generates a matrix that contains amino acids in library region
l <- length(codons)
AAs <- matrix(0,l,lib)
Phereg <- gregexpr("\\s(TT[TC])",codons,useBytes = FALSE)
  lphe <- length(Phereg)
  for(a in c(1:lphe)){
    lphee <- length(Phereg[[a]])
    for(b in c(1:lphee)){
      value <- Phereg[[a]][b]
      if(value > initialcodon && value < endcodon){
        AAs[a,(value%/%4 - (del-1))] <- "F"
      }  
    }
  }
Leureg <- gregexpr("(\\sTT[AG])|(\\sCT[GACT])",codons,useBytes = FALSE)
  lleu <- length(Leureg)
  for(a in c(1:lleu)){
    lleuu <- length(Leureg[[a]])
    for(b in c(1:lleuu)){
      value <- Leureg[[a]][b]
      if(value > initialcodon && value < endcodon){
        AAs[a,(value%/%4 - (del-1))] <- "L"
      }  
    }
  }
Serreg <- gregexpr("(\\sTC[GCAT])|(\\sAG[TC])",codons,useBytes = FALSE)
  lser <- length(Serreg)
  for(a in c(1:lser)){
    lserr <- length(Serreg[[a]])
    for(b in c(1:lserr)){
      value <- Serreg[[a]][b]
      if(value > initialcodon && value < endcodon){
        AAs[a,(value%/%4 - (del-1))] <- "S"
      }  
    }
  }
Tyrreg <- gregexpr("\\sTA[TC]",codons,useBytes = FALSE)
  ltyr <- length(Tyrreg)
  for(a in c(1:ltyr)){
    ltyrr <- length(Tyrreg[[a]])
    for(b in c(1:ltyrr)){
      value <- Tyrreg[[a]][b]
      if(value > initialcodon && value < endcodon){
        AAs[a,(value%/%4 - (del-1))] <- "Y"
      }  
    }
  }
TAGgreg <- gregexpr("\\sTAG",codons,useBytes = FALSE)
  ltag <- length(TAGgreg)
  for(a in c(1:ltag)){
    ltagg <- length(TAGgreg[[a]])
    for(b in c(1:ltagg)){
      value <- TAGgreg[[a]][b]
      if(value > initialcodon && value < endcodon){
        AAs[a,(value%/%4 - (del-1))] <- "TAG"
      }  
    }
  }
TAAgreg <- gregexpr("\\sTAA",codons,useBytes = FALSE)
  ltaa <- length(TAAgreg)
  for(a in c(1:ltaa)){
    ltaaa <- length(TAAgreg[[a]])
    for(b in c(1:ltaaa)){
      value <- TAAgreg[[a]][b]
      if(value > initialcodon && value < endcodon){
        AAs[a,(value%/%4 - (del-1))] <- NA
      }  
    }
  }
Cysreg <- gregexpr("\\sTG[TC]",codons,useBytes = FALSE)
  lcys <- length(Cysreg)
  for(a in c(1:lcys)){
    lcyss <- length(Cysreg[[a]])
    for(b in c(1:lcyss)){
      value <- Cysreg[[a]][b]
      if(value > initialcodon && value < endcodon){
        AAs[a,(value%/%4 - (del-1))] <- "C"
      }  
    }
  }
TGAgreg <- gregexpr("\\sTGA",codons,useBytes = FALSE)
  ltga <- length(TGAgreg)
  for(a in c(1:ltga)){
    ltgaa <- length(TGAgreg[[a]])
    for(b in c(1:ltgaa)){
      value <- TGAgreg[[a]][b]
      if(value > initialcodon && value < endcodon){
        AAs[a,(value%/%4 - (del-1))] <- NA
      }  
    }
  }
Trpgreg <- gregexpr("\\sTGG",codons,useBytes = FALSE)
  ltrp <- length(Trpgreg)
  for(a in c(1:ltrp)){
    ltrpp <- length(Trpgreg[[a]])
    for(b in c(1:ltrpp)){
      value <- Trpgreg[[a]][b]
      if(value > initialcodon && value < endcodon){
        AAs[a,(value%/%4 - (del-1))] <- "W"
      }  
    }
  }
Progreg <- gregexpr("\\sCC[GCAT]",codons,useBytes = FALSE)
  lpro <- length(Progreg)
  for(a in c(1:lpro)){
    lproo <- length(Progreg[[a]])
    for(b in c(1:lproo)){
      value <- Progreg[[a]][b]
      if(value > initialcodon && value < endcodon){
        AAs[a,(value%/%4 - (del-1))] <- "P"
      }  
    }
  }
Hisgreg <- gregexpr("\\sCA[CT]",codons,useBytes = FALSE)
  lhis <- length(Hisgreg)
  for(a in c(1:lhis)){
    lhiss <- length(Hisgreg[[a]])
    for(b in c(1:lhiss)){
      value <- Hisgreg[[a]][b]
      if(value > initialcodon && value < endcodon){
        AAs[a,(value%/%4 - (del-1))] <- "H"
      }  
    }
  }
Glngreg <- gregexpr("\\sCA[AG]",codons,useBytes = FALSE)
  lgln <- length(Glngreg)
  for(a in c(1:lgln)){
    lglnn <- length(Glngreg[[a]])
    for(b in c(1:lglnn)){
      value <- Glngreg[[a]][b]
      if(value > initialcodon && value < endcodon){
        AAs[a,(value%/%4 - (del-1))] <- "Q"
      }  
    }
  }
Arggreg <- gregexpr("(\\sCG[GCAT])|(\\sAG[GA])",codons,useBytes = FALSE)
  larg <- length(Arggreg)
  for(a in c(1:larg)){
    largg <- length(Arggreg[[a]])
    for(b in c(1:largg)){
      value <- Arggreg[[a]][b]
      if(value > initialcodon && value < endcodon){
        AAs[a,(value%/%4 - (del-1))] <- "R"
      }  
    }
  }
Ilegreg <- gregexpr("\\sAT[CAT]",codons,useBytes = FALSE)
  lile <- length(Ilegreg)
  for(a in c(1:lile)){
    lilee <- length(Ilegreg[[a]])
    for(b in c(1:lilee)){
      value <- Ilegreg[[a]][b]
      if(value > initialcodon && value < endcodon){
        AAs[a,(value%/%4 - (del-1))] <- "I"
      }  
    }
  }
Metgreg <- gregexpr("\\sATG",codons,useBytes = FALSE)
  lmet <- length(Metgreg)
  for(a in c(1:lmet)){
    lmett <- length(Metgreg[[a]])
    for(b in c(1:lmett)){
      value <- Metgreg[[a]][b]
      if(value > initialcodon && value < endcodon){
        AAs[a,(value%/%4 - (del-1))] <- "M"
      }  
    }
  }
Thrgreg <- gregexpr("\\sAC[GCAT]",codons,useBytes = FALSE)
  lthr <- length(Thrgreg)
  for(a in c(1:lthr)){
    lthrr <- length(Thrgreg[[a]])
    for(b in c(1:lthrr)){
      value <- Thrgreg[[a]][b]
      if(value > initialcodon && value < endcodon){
        AAs[a,(value%/%4 - (del-1))] <- "T"
      }  
    }
  }
Asngreg <- gregexpr("\\sAA[CT]",codons,useBytes = FALSE)
  lasn <- length(Asngreg)
  for(a in c(1:lasn)){
    lasnn <- length(Asngreg[[a]])
    for(b in c(1:lasnn)){
      value <- Asngreg[[a]][b]
      if(value > initialcodon && value < endcodon){
        AAs[a,(value%/%4 - (del-1))] <- "N"
      }  
    }
  }
Lysgreg <- gregexpr("\\sAA[AG]",codons,useBytes = FALSE)
  llys <- length(Lysgreg)
  for(a in c(1:llys)){
    llyss <- length(Lysgreg[[a]])
    for(b in c(1:llyss)){
      value <- Lysgreg[[a]][b]
      if(value > initialcodon && value < endcodon){
        AAs[a,(value%/%4 - (del-1))] <- "K"
      }  
    }
  }
Valgreg <- gregexpr("\\sGT[GACT]",codons,useBytes = FALSE)
  lval <- length(Valgreg)
  for(a in c(1:lval)){
    lvall <- length(Valgreg[[a]])
    for(b in c(1:lvall)){
      value <- Valgreg[[a]][b]
      if(value > initialcodon && value < endcodon){
        AAs[a,(value%/%4 - (del-1))] <- "V"
      }  
    }
  }
Alagreg <- gregexpr("\\sGC[GACT]",codons,useBytes = FALSE)
  lala <- length(Alagreg)
  for(a in c(1:lala)){
    lalaa <- length(Alagreg[[a]])
    for(b in c(1:lalaa)){
      value <- Alagreg[[a]][b]
      if(value > initialcodon && value < endcodon){
        AAs[a,(value%/%4 - (del-1))] <- "A"
      }  
    }
  }
Aspgreg <- gregexpr("\\sGA[TC]",codons,useBytes = FALSE)
  lasp <- length(Aspgreg)
  for(a in c(1:lasp)){
    laspp <- length(Aspgreg[[a]])
    for(b in c(1:laspp)){
      value <- Aspgreg[[a]][b]
      if(value > initialcodon && value < endcodon){
        AAs[a,(value%/%4 - (del-1))] <- "D"
      }  
    }
  }
Glugreg <- gregexpr("\\sGA[AG]",codons,useBytes = FALSE)
  lglu <- length(Glugreg)
  for(a in c(1:lglu)){
    lgluu <- length(Glugreg[[a]])
    for(b in c(1:lgluu)){
      value <- Glugreg[[a]][b]
      if(value > initialcodon && value < endcodon){
        AAs[a,(value%/%4 - (del-1))] <- "E"
      }  
    }
  }
Glygreg <- gregexpr("\\sGG[GACT]",codons,useBytes = FALSE)
  lgly <- length(Glygreg)
  for(a in c(1:lgly)){
    lglyy <- length(Glygreg[[a]])
    for(b in c(1:lglyy)){
      value <- Glygreg[[a]][b]
      if(value > initialcodon && value < endcodon){
        AAs[a,(value%/%4 - (del-1))] <- "G"
      }  
    }
  }
AAs <- as.data.frame(AAs)
countAAs <- apply(AAs, 2, table)

#this counts number of unique amino acid sequences and sorts them by frequency
UniqueAAs <- AAs %>% group_by_all() %>% count()
UniqueAAs <- UniqueAAs[order(-UniqueAAs$n),]
UniqueAAs <- UniqueAAs[apply(UniqueAAs,1,function(row) all(row != 0)),]
UniqueAAs <- na.omit(UniqueAAs)

#this counts sequences that have TAG codons, sequences that have more than one are only counted once
TAGreg <- regexpr("\\sTAG",codons)
TAGtable <- table(TAGreg)
percentTAG <- sum(TAGtable[2:length(TAGtable)])/length(codons)*100#percent of sequences containing TAG

#this creates heatmap for amino acid frequency per library position, change scale according to values
AAtable <- apply(AAs, 2, table)
AAtable <- as.matrix(AAtable/length(codons))
colnames(AAtable) <- c(1:lib)
heatmapcolors <- colorRampPalette(brewer.pal(9,"Blues"))(100)
sc <- seq(0.0,0.3,by=0.003)
AAheatmap <- heatmap.2(AAtable, Rowv = NA, Colv = NA, col = heatmapcolors, density.info = "none", scale = "none", trace = "none", breaks = sc,  xlab = "Position in Library", ylab = "Codon", margins = c(3,4), dendrogram = "none")

#this creates projected heatmap based on NNK randomized codons
randomAAs <- matrix(0,21,lib,dimnames = list(rownames(AAtable),c(1:lib)))
randomAAs[c("A","G","P","T","V"),] <- 2/32
randomAAs[c("C","H","Q","N","K","Y","D","E","W","I","M","TAG","F"),] <- 1/32
randomAAs[c("L","S","R"),] <- 3/32
NNKheatmap <- heatmap.2(randomAAs, Rowv = NA, Colv = NA, col = heatmapcolors, density.info = "none", scale = "none", trace = "none", breaks = sc,  xlab = "Position in Library", ylab = "Codon", margins = c(3,4), dendrogram = "none")

#this creates heatmap showing bias from random, change scale with respect to range of values
lscale <- seq(-1,4,by=5/100)
librarybias <- (AAtable - randomAAs)/randomAAs
Biasheatmap <- heatmap.2(librarybias, Rowv = NA, Colv = NA, col = heatmapcolors, density.info = "none", scale = "none", trace = "none", breaks = lscale,  xlab = "Position in Library", ylab = "Codon", margins = c(3,4), dendrogram = "none")