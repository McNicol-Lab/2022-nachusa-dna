library(dada2); packageVersion("dada2") 
library(ggplot2)
library(phyloseq)
library(DECIPHER) #for assigning taxonomy 
library(phangorn) #for tree building 
require(tidyverse) #data wrangling 
library(dplyr) #data wrangling 
library(rmarkdown)

setwd("~/Desktop/sarah_proj/") #set working directory 

#load file paths 

forward_path <- file.path("data/forward/")
reverse_path <- file.path("data/reverse/")

forward <- sort(list.files(forward_path, full.names = T))
reverse <- sort(list.files(reverse_path, full.names = T))

#check quality
Forwardi <- sample(length(forward), 4) #the number should be the number of files in the directory with fastq files

for(i in Forwardi) { 
  print(plotQualityProfile(forward[i]))  
} 


Reversei <- sample(length(reverse), 4) #the number should be the number of files in the directory with fastq files

for(i in Reversei) { 
  print(plotQualityProfile(reverse[i]))  
} 


#Filter and trimm primers
all_path <- file.path("all/")
filt_path <- file.path("data/", "filtered")

fns <- sort(list.files(all_path, full.names = TRUE))
fnFs <- fns[grepl("R1", fns)]
fnRs <- fns[grepl("R2", fns)]


if(!file_test("-d", filt_path)) dir.create(filt_path)

filtFs <- file.path(filt_path, basename(fnFs))
filtRs <- file.path(filt_path, basename(fnRs))

for(i in seq_along(fnFs)) {
  fastqPairedFilter(c(fnFs[[i]], fnRs[[i]]),
                    c(filtFs[[i]], filtRs[[i]]),
                    trimLeft=c(19, 20), truncLen=c(150, 150),
                    maxN=0, maxEE=2, truncQ=2,
                    compress=TRUE)}


#check quality after filtering 
  
  filtFi <- sample(length(filtFs), 4) #the number should be the number of files in the directory with fastq files
  
  for(i in filtFi) { 
    print(plotQualityProfile(filtFs[i]))  
  }   
  
  
  
  filtRi <- sample(length(filtRs), 4) #the number should be the number of files in the directory with fastq files
  
  for(i in filtRi) { 
    print(plotQualityProfile(filtRs[i]))  
  }   
  
#Infer Seq variants 
  derepsFs <- derepFastq(filtFs)
  deprepRs <- derepFastq(filtRs)
  derepsRs <- deprepRs #typo 
  sam.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
  names(derepsFs) <- sam.names
  names(derepsRs) <- sam.names

ddF <- dada(derepsFs[1:4], err=NULL, selfConsist=TRUE)
ddR <- dada(derepsRs[1:4], err = NULL, selfConsist = TRUE)

plotErrors(ddF)
plotErrors(ddR)


dadaFs <- dada(derepsFs, err=ddF[[1]]$err_out, pool=TRUE)

dadaRs <- dada(derepsRs, err=ddR[[1]]$err_out, pool=TRUE)


mergers <- mergePairs(dadaFs, derepsFs, dadaRs, derepsRs)

mergers2 <- mergePairs(dadaFs, derepsFs, dadaRs, derepsRs, minOverlap = 6)


seqtab.all <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))]) #make seq tab 
seqtab <- removeBimeraDenovo(seqtab.all) #remove chimeras 

seqtab.all2 <- makeSequenceTable(mergers2[!grepl("Mock", names(mergers2))]) #make seq tab 
seqtab <- removeBimeraDenovo(seqtab.all) #remove chimeras 

