#this gives overall enrichment of each peptide sequence by comparing Round 1 and Round 4 Library
library(prodlim)
path <- "/Users/traehampton/Documents/Research/Sequencing Results/Next Gen Sequencing/21510Wns_N21170"
lib <- 7
match <- row.match(UniqueAAsR1[,1:lib],UniqueAAsR3[,1:lib])
matchseq <- which(is.na(match)==FALSE)
percentenriched <- (UniqueAAsR3[match[matchseq],lib+1]/sum(UniqueAAsR3[,lib+1])-UniqueAAsR1[matchseq,lib+1]/sum(UniqueAAsR1$n))/(UniqueAAsR1[matchseq,lib+1]/sum(UniqueAAsR1$n))
enrichedseq <- UniqueAAsR3[match[matchseq],]
enrichedseq$enrichment <- percentenriched[,1]
enrichedseq <- enrichedseq[order(-enrichedseq$enrichment),]
enrichedseq <- as.data.frame(enrichedseq)
UniqueAAsR1$percent <- UniqueAAsR1$n/sum(UniqueAAsR1$n)*100
UniqueAAsR3$percent <- UniqueAAsR3$n/sum(UniqueAAsR3$n)*100
#Change names to whatever you want
write.csv(enrichedseq, paste(path,"/HDAC1_R3vR1.csv", sep = ""), row.names = F)
write.csv(UniqueAAsR1, paste(path,"/HDAC1_R1.csv", sep = ""), row.names = F)
write.csv(UniqueAAsR3, paste(path,"/HDAC1_R3.csv", sep = ""), row.names = F)
write.csv(TAGtable, paste(path,"/12mer_plus_neg_LibraryBias.csv", sep = ""), row.names = F)
write.csv(UniqueAAsR3, paste(path,"/KetoKR3.csv", sep = ""), row.names = F)