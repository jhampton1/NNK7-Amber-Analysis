#NNK Analysis of Sequences in library,ignoring TAG codons, do this after defining TAG positions in AAs
#Make sure variables are defined from the AA analysis
l <- length(codons)
matrixseq <- matrix(0,l,3*lib)
Ggreg <- gregexpr("G",codons,useBytes = FALSE)
  lg <- length(Ggreg)
  for(a in c(1:lg)){
    lgg <- length(Ggreg[[a]])
    for(b in c(1:lgg)){
      value <- Ggreg[[a]][b]
      if(value > initialcodon && value < endcodon+1){
        if(AAs[a,(value%/%4 - (del-1))] != "TAG" | is.na(AAs[a,(value%/%4 - (del-1))])){
          matrixseq[a,(value-del*3-(value%/%4))] <- "G"
        }
      }
    }
  }
Cgreg <- gregexpr("C",codons,useBytes = FALSE)
  lc <- length(Cgreg)
  for(a in c(1:lc)){
    lcc <- length(Cgreg[[a]])
    for(b in c(1:lcc)){
      value <- Cgreg[[a]][b]
      if(value > initialcodon && value < endcodon+1){
        if(AAs[a,(value%/%4 - (del-1))] != "TAG" | is.na(AAs[a,(value%/%4 - (del-1))])){
          matrixseq[a,(value-del*3-(value%/%4))] <- "C"
        }
      }
    }
  }
Agreg <- gregexpr("A",codons,useBytes = FALSE)
  la <- length(Agreg)
  for(a in c(1:la)){
    laa <- length(Agreg[[a]])
    for(b in c(1:laa)){
      value <- Agreg[[a]][b]
      if(value > initialcodon && value < endcodon+1){
        if(AAs[a,(value%/%4 - (del-1))] != "TAG" | is.na(AAs[a,(value%/%4 - (del-1))])){
          matrixseq[a,(value-del*3-(value%/%4))] <- "A"
        }
      }
    }
  }
Tgreg <- gregexpr("T",codons,useBytes = FALSE)
  lt <- length(Tgreg)
  for(a in c(1:lt)){
    ltt <- length(Tgreg[[a]])
    for(b in c(1:ltt)){
      value <- Tgreg[[a]][b]
      if(value > initialcodon && value < endcodon+1){
        if(AAs[a,(value%/%4 - (del-1))] != "TAG" | is.na(AAs[a,(value%/%4 - (del-1))])){
          matrixseq[a,(value-del*3-(value%/%4))] <- "T"
        }
      }
    }
  }
matrixseq <- as.data.frame(matrixseq)

#counts nucleotides for all of the positions in matrixseq, returns a table with amount of nucleotides at each position
count <- apply(matrixseq, 2, table)
View(count)