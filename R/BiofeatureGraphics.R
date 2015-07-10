##########################################################################
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software

##########################################################################

##########################################################################
# File: BiofeatureGraphics.R
# Author: Tiphaine Martin
# Email: tiphaine.martin@kcl.ac.uk
# Purpose: coMET allows the display of p-values from association
#           with a correlation heatmap.
# Version : 0.99.9
###########################################################################

#-------------------- CpG pvalue ------------------
cpgPvalue<-function(cprange,data,chr,start,end,typefunction,title){
  DataTrack(range=cprange,start,end,data=data,name="CpG pvalue")
}

#-------------------- Convertion chromosome UCSC to ENSEMBL ------------------
chrUCSC2ENSEMBL<-function(chr){
  gsub("^chr", "", chr)
}
#-------------------- CREATION track all elements but only one line per element ------------------
genesNameENSEMBL<-function(gen,chr,start,end,dataset){
  if(is.null(chr)){
    stop("Invalid in function genesENSEMB :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid in function genesENSEMB :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid in function genesENSEMB :end null:\n")
  }
  if(is.null(gen)){
    stop("Invalid in function genesENSEMB :gen null:\n")
  }
  chrEnsembl=chrUCSC2ENSEMBL(chr)
  
  genTrunk <- gsub("\\..*","",gen)
  
  ens_ENSEMBL <- NULL
  if(!is.na(match(tolower(genTrunk), c("grch37","hg19")))){
    martENSEMBL=useMart(host='grch37.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL',
                        dataset=dataset)
    ens_ENSEMBL <- getBM(c("ensembl_gene_id","external_gene_name"),
                         filters = c("chromosome_name","start","end"),
                         values = list(chrEnsembl, start, end), mart=martENSEMBL) 
  } else {
    martENSEMBL=useMart("ensembl",dataset=dataset)
    ens_ENSEMBL <- getBM(c("ensembl_gene_id","external_gene_name"),
                         filters = c("chromosome_name","start","end"),
                         values = list(chrEnsembl, start, end), mart=martENSEMBL) 
  }
  
  
  
  if(nrow(ens_ENSEMBL) == 0) {
    ens_ENSEMBL <- NULL
  } 
  ens_ENSEMBL
}

#-------------------- CREATION track all elements but only one line per element ------------------
genesENSEMBL<-function(gen,chr,start,end,showId=FALSE){
  if(is.null(chr)){
    stop("Invalid in function genesENSEMBL :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid in function genesENSEMBL :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid in function genesENSEMBL :end null:\n")
  }
  if(is.null(gen)){
    stop("Invalid in function genesENSEMBL :gen null:\n")
  }
  #cat("data",gen,"\t",chr,"\t",start,"\t",end,"\n")
  
  genTrunk <- gsub("\\..*","",gen)
  
  biomTrack=NULL
  if(!is.na(match(tolower(genTrunk), c("grch37","hg19")))){
    martENSEMBL=useMart(host='grch37.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL',
                        dataset='hsapiens_gene_ensembl')
    #  fm <- Gviz:::.getBMFeatureMap()
    #  fm["symbol"] <- "external_gene_id"
    #  biomTrack <- BiomartGeneRegionTrack(genome = gen, featureMap=fm, biomart=martENSEMBL,
    #                                      chromosome = chr, start = start, 
    #                                      end = end,  name = "ENSEMBL",
    #                                      groupAnnotation = "group",
    #                                      just.group = "above",
    #                                     fontcolor="black",showId=showId,size=2)
    biomTrack <- BiomartGeneRegionTrack(genome = genTrunk, biomart=martENSEMBL,
                                        chromosome = chr, start = start, 
                                        end = end,  name = "genes ENSEMBL",
                                        groupAnnotation = "group",
                                        just.group = "above",
                                        fontcolor="black",showId=showId,size=2,
                                        col.line = NULL, col = NULL)
    
    
  } else {
    martENSEMBL=useMart("ensembl",dataset='hsapiens_gene_ensembl')
    biomTrack <- BiomartGeneRegionTrack(genome = genTrunk, biomart=martENSEMBL,
                                        chromosome = chr, start = start, 
                                        end = end,  name = "genes ENSEMBL",
                                        groupAnnotation = "group",
                                        just.group = "above",
                                        fontcolor="black",showId=showId,size=2,
                                        col.line = NULL, col = NULL)
  }
  
  
  #  cat("change feature\n")
  if(length(feature(biomTrack)) > 0) {
    feature(biomTrack) <- "protein_coding"
    #  cat("change elements\n")
    r <- split(ranges(biomTrack), gene(biomTrack))
    # cat("change elements1\n")
    rNew <-  endoapply(r, function(x){
      rx <- reduce(x, with.revmap=TRUE) 
      mcols(rx) <- mcols(x)[sapply(rx$revmap, head, 1),]
      rx$transcript <- rx$gene
      rx
    })
    #cat("change elements2\n")
    ranges(biomTrack) <- unlist(rNew)
  } 
  #cat("change elements3\n")
  
  biomTrack
}

#-------------------- CREATION track all elements but only one line per element with specific color for features of interest ------------------
interestGenesENSEMBL<-function(gen,chr,start,end,interestfeatures,interestcolor,showId=FALSE){
  if(is.null(chr)){
    stop("Invalid in function interestGenesENSEMBL :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid in function interestGenesENSEMBL :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid in function interestGenesENSEMBL :end null:\n")
  }
  if(is.null(gen)){
    stop("Invalid in function interestGenesENSEMBL :gen null:\n")
  }
  if(is.null(interestfeatures)){
    stop("Invalid in function interestGenesENSEMBL: interestfeatures null:\n")
  }
  if(is.null(interestcolor)){
    stop("Invalid in function interestGenesENSEMBL: interestcolor null:\n")
  }
  #cat("data",gen,"\t",chr,"\t",start,"\t",end,"\n")
  
  genTrunk <- gsub("\\..*","",gen)
  
  biomTrack=NULL
  if(!is.na(match(tolower(genTrunk), c("grch37","hg19")))){
    martENSEMBL=useMart(host='grch37.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL',
                        dataset='hsapiens_gene_ensembl')
    #  fm <- Gviz:::.getBMFeatureMap()
    #  fm["symbol"] <- "external_gene_id"
    #  biomTrack <- BiomartGeneRegionTrack(genome = gen, featureMap=fm, biomart=martENSEMBL,
    #                                      chromosome = chr, start = start, 
    #                                      end = end,  name = "ENSEMBL",
    #                                      groupAnnotation = "group",
    #                                      just.group = "above",
    #                                     fontcolor="black",showId=showId,size=2)
    biomTrack <- BiomartGeneRegionTrack(genome = genTrunk, biomart=martENSEMBL,
                                        chromosome = chr, start = start, 
                                        end = end,  name = "genes ENSEMBL of interest",
                                        groupAnnotation = "group",
                                        just.group = "above",
                                        fontcolor="black",showId=showId,size=2,
                                        col.line = NULL, col = NULL)
    
    
  } else {
    martENSEMBL=useMart("ensembl",dataset='hsapiens_gene_ensembl')
    biomTrack <- BiomartGeneRegionTrack(genome = genTrunk, biomart=martENSEMBL,
                                        chromosome = chr, start = start, 
                                        end = end,  name = "genes ENSEMBL of interest",
                                        groupAnnotation = "group",
                                        just.group = "above",
                                        fontcolor="black",showId=showId,size=2,
                                        col.line = NULL, col = NULL)
  }
  
  
  #  cat("change feature\n")
  if(length(feature(biomTrack)) > 0) {
    feature(biomTrack) <- "protein_coding"
    #  cat("change elements\n")
    r <- split(ranges(biomTrack), gene(biomTrack))
    # cat("change elements1\n")
    rNew <-  endoapply(r, function(x){
      rx <- reduce(x, with.revmap=TRUE) 
      mcols(rx) <- mcols(x)[sapply(rx$revmap, head, 1),]
      rx$transcript <- rx$gene
      rx
    })
    #cat("change elements2\n")
    ranges(biomTrack) <- unlist(rNew)
  } 
  #cat("change elements3\n")
  
  #Color of interest features
  myfeatures <- feature(biomTrack)
  for( f in 1:length(myfeatures)){
    # cat("change features ",f,"\n")
    for (i in 1:nrow(as.data.frame(interestfeatures))){
      #    cat("change interest features ",i,"\n")
      if (start(biomTrack)[f] == as.numeric(interestfeatures[i,1]) & end(biomTrack)[f] == as.numeric(interestfeatures[i,2])){
        #    cat("change color",f,":",i,"(",start(biomTrack)[f],",",end(biomTrack)[f],") \n")
        feature(biomTrack)[f] <- as.character(interestfeatures[i,3])
      }
    }
  }
  
  displayPars(biomTrack) <- interestcolor
  
  biomTrack
}


#-------------------- CREATION track for all transcript in the region ------------------
transcriptENSEMBL<-function(gen,chr,start,end,showId=FALSE){
  if(is.null(chr)){
    stop("Invalid in function transcriptENSEMBL :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid in function transcriptENSEMBL :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid in function transcriptENSEMBL :end null:\n")
  }
  if(is.null(gen)){
    stop("Invalid in function transcriptENSEMBL :gen null:\n")
  }
  
  genTrunk <- gsub("\\..*","",gen)
  
  biomTrack=NULL
  if(!is.na(match(tolower(genTrunk), c("grch37","hg19")))){
    martENSEMBL=useMart(host='grch37.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL',
                        dataset='hsapiens_gene_ensembl')
    #  fm <- Gviz:::.getBMFeatureMap()
    #  fm["symbol"] <- "external_gene_id"
    #  biomTrack <- BiomartGeneRegionTrack(genome = gen, featureMap=fm, biomart=martENSEMBL,
    #                                      chromosome = chr, start = start, 
    #                                    end = end,  name = "ENSEMBL",
    #                                      fontcolor="black",groupAnnotation = "group",
    #                                     just.group = "above",showId=showId,size=2)
    biomTrack <- BiomartGeneRegionTrack(genome = genTrunk, biomart=martENSEMBL,
                                        chromosome = chr, start = start, 
                                        end = end,  name = "transcripts ENSEMBL",
                                        fontcolor="black",groupAnnotation = "group",
                                        just.group = "above",showId=showId,size=2,
                                        col.line = NULL, col = NULL)
    
  } else {
    martENSEMBL=useMart("ensembl",dataset='hsapiens_gene_ensembl')
    biomTrack <- BiomartGeneRegionTrack(genome = genTrunk, biomart=martENSEMBL,
                                        chromosome = chr, start = start, 
                                        end = end,  name = "transcripts ENSEMBL",
                                        fontcolor="black", groupAnnotation = "group",
                                        just.group = "above",showId=showId,size=2,
                                        col.line = NULL, col = NULL)
  }
  
  #cat("data",gen,"\t",chr,"\t",start,"\t",end,"\n")
  
  
  #stacking="dense"
  
  biomTrack
}


#-------------------- CREATION track for all transcript in the region with color for features of interest------------------
interestTranscriptENSEMBL<-function(gen,chr,start,end,interestfeatures,interestcolor,showId=FALSE){
  if(is.null(chr)){
    stop("Invalid in function interestTranscriptENSEMBL :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid in function interestTranscriptENSEMBL: start null:\n")
  }
  if(is.null(end)){
    stop("Invalid in function interestTranscriptENSEMBL: end null:\n")
  }
  if(is.null(gen)){
    stop("Invalid in function interestTranscriptENSEMBL: gen null:\n")
  }
  if(is.null(interestfeatures)){
    stop("Invalid in function interestTranscriptENSEMBL: interestfeatures null:\n")
  }
  if(is.null(interestcolor)){
    stop("Invalid in function interestTranscriptENSEMBL: interestcolor null:\n")
  }
  
  genTrunk <- gsub("\\..*","",gen)
  
  biomTrack=NULL
  if(!is.na(match(tolower(genTrunk), c("grch37","hg19")))){
    martENSEMBL=useMart(host='grch37.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL',
                        dataset='hsapiens_gene_ensembl')
    #  fm <- Gviz:::.getBMFeatureMap()
    #  fm["symbol"] <- "external_gene_id"
    #  biomTrack <- BiomartGeneRegionTrack(genome = gen, featureMap=fm, biomart=martENSEMBL,
    #                                      chromosome = chr, start = start, 
    #                                    end = end,  name = "ENSEMBL",
    #                                      fontcolor="black",groupAnnotation = "group",
    #                                     just.group = "above",showId=showId,size=2)
    biomTrack <- BiomartGeneRegionTrack(genome = genTrunk, biomart=martENSEMBL,
                                        chromosome = chr, start = start, 
                                        end = end,  name = "transcript ENSEMBL of interest",
                                        fontcolor="black",groupAnnotation = "group",
                                        just.group = "above",showId=showId,size=2,
                                        col.line = NULL, col = NULL, collapse= FALSE)
    
  } else {
    martENSEMBL=useMart("ensembl",dataset='hsapiens_gene_ensembl')
    biomTrack <- BiomartGeneRegionTrack(genome = genTrunk, biomart=martENSEMBL,
                                        chromosome = chr, start = start, 
                                        end = end,  name = "transcript ENSEMBL of interest",
                                        fontcolor="black", groupAnnotation = "group",
                                        just.group = "above",showId=showId,size=2,
                                        col.line = NULL, col = NULL, collapse= FALSE)
  }
  
  #cat("data",gen,"\t",chr,"\t",start,"\t",end,"\n")
  
  
  #stacking="dense"
  #Color of interest features
  myfeatures <- feature(biomTrack)
  for( f in 1:length(myfeatures)){
   # cat("change features ",f,"\n")
    for (i in 1:nrow(as.data.frame(interestfeatures))){
  #    cat("change interest features ",i,"\n")
      if (start(biomTrack)[f] == as.numeric(interestfeatures[i,1]) & end(biomTrack)[f] == as.numeric(interestfeatures[i,2])){
    #    cat("change color",f,":",i,"(",start(biomTrack)[f],",",end(biomTrack)[f],") \n")
        feature(biomTrack)[f] <- as.character(interestfeatures[i,3])
      }
    }
  }
  
  displayPars(biomTrack) <- interestcolor
  
  biomTrack
}



#-------------------- CREATION track all type of chromatineHMM from UCSC ------------------
chromatinHMMAll<-function(gen,chr,start,end,mySession,track.name="Broad ChromHMM",pattern=NULL,table.name=NULL){
  if(is.null(chr)){
    stop("Invalid in function chromatinHMMAll :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid in function chromatinHMMAll :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid in function chromatinHMMAll :end null:\n")
  }
  if(is.null(gen)){
    stop("Invalid in function chromatinHMMAll :gen null:\n")
  }
  
  genTrunk <- gsub("\\..*","",gen)
  
  if(is.null(track.name) & (!is.na((match(genTrunk, tolower(c("hg19","grch37"))))))){
    track.name="Broad ChromHMM"
    
  }else if(is.null(track.name) & genTrunk != "hg19"){
    stop("Invalid in function chromatinHMMAll :track.namenull:\n")
  }
  tablestrack<-tableNames(ucscTableQuery(mySession, track=track.name))
  if(is.null(pattern)) {
    patterntable<-1:length(tablestrack)
  } else{
    patterntable<-grep(pattern, tablestrack,ignore.case=TRUE)
  }
  
  lltrack=list()
  for(i in patterntable){
    table.name<-tablestrack[i]
    tmp<-chromatinHMMOne(genTrunk,chr,start,end,mySession,track.name,table.name)
    if(!is.null(tmp) | length(feature(tmp)) > 0){
      lltrack=c(lltrack,tmp)
    }
  }
  if(length(lltrack) == 0){
    data_trackfunc <- AnnotationTrack()
    chromosome(data_trackfunc) <- chr
    start(data_trackfunc) <- start
    end(data_trackfunc) <- end
    lltrack=list(data_trackfunc)
  }
  lltrack
}

#-------------------- CREATION track one type of chromaHMM from UCSC ------------------
chromatinHMMOne<-function(gen,chr,start,end,mySession,track.name="Broad ChromHMM", table.name=NULL){
  if(is.null(chr)){
    stop("Invalid in function chromatinHMMOne :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid in function chromatinHMMOne :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid in function chromatinHMMOne :end null:\n")
  }
  if(is.null(gen)){
    stop("Invalid in function chromatinHMMOne :gen null:\n")
  }
  if(is.null(track.name) & gen == "hg19"){
    track.name="Broad ChromHMM"
  }else if(is.null(track.name)){
    stop("Invalid in function chromatinHMMOne :track.namenull:\n")
  }
  if(is.null(table.name)){
    table.name="wgEncodeBroadHmmHsmmHMM"
  }else if(is.null(table.name)){
    stop("Invalid in function chromatinHMMOne :gen null:\n")
  }
  mygrange <- GRanges(chr, IRanges(start, end))
  dataUCSC <- getTable(ucscTableQuery (mySession, range=mygrange, 
                                       track=track.name, table=table.name))
  
  data_trackfunc <- AnnotationTrack()
  if(nrow(dataUCSC) > 0) {
    data_trackfunc <- AnnotationTrack(chromosome=dataUCSC[,"chrom"],strand ="*",
                                      start=dataUCSC[,"chromStart"],end=dataUCSC[,"chromEnd"],
                                      feature=dataUCSC[,"name"],group=dataUCSC[,"name"],
                                      id=dataUCSC[,"name"], name = "Broad chromatinHMM",
                                      stacking="dense", col.line = NULL, col = NULL)
    chromosome(data_trackfunc) <- chr
    
    displayPars(data_trackfunc) <- list(
      "1_Active_Promoter"= "firebrick1",
      "2_Weak_Promoter"="darksalmon"  ,
      "3_Poised_Promoter"="blueviolet",
      "4_Strong_Enhancer"= "Orange",
      "5_Strong_Enhancer"= "coral",
      "6_Weak_Enhancer"="yellow",
      "7_Weak_Enhancer"="gold",
      "8_Insulator"="cornflowerblue",
      "9_Txn_Transition"="darkolivegreen",
      "10_Txn_Elongation"="forestgreen",
      "11_Weak_Txn"="darkseagreen1",
      "12_Repressed"="gainsboro",
      "13_Heterochrom/lo"="gray74",
      "14_Repetitive/CNV"="gray77",
      "15_Repetitive/CNV"="gray86")
  } else {
    data_trackfunc <- AnnotationTrack()
    chromosome(data_trackfunc) <- chr
    start(data_trackfunc) <- start
    end(data_trackfunc) <- end
  }
  
  data_trackfunc
}

#-------------------- CREATION track all types of Histone density from UCSC ------------------
HistoneAll<-function(gen,chr,start,end,mySession,pattern=NULL,track.name="Broad Histone",table.name=NULL){
  if(is.null(chr)){
    stop("Invalid in function HistoneAll :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid in function HistoneAll :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid in function HistoneAll :end null:\n")
  }
  if(is.null(gen)){
    stop("Invalid in function HistoneAll :gen null:\n")
  }
  if(is.null(track.name) & (!is.na(match(gen, tolower(c("hg19","grch37")))))){
    track.name="Broad Histone"
  }else if(is.null(track.name) & gen != "hg19"){
    stop("Invalid in function HistoneAll :track.namenull:\n")
  }
  tablestrack<-tableNames(ucscTableQuery(mySession, track=track.name))
  
  #patternPk="PkV*[0-9]*$"
  patternPk="Pk$"
  patterntablePk<-grep(patternPk, tablestrack,ignore.case=TRUE)
  
  patternPkV2="PkV2$"
  patterntablePkV2 <- grep(patternPkV2, tablestrack,ignore.case=TRUE)
  
  patterntablePktotal <- c(patterntablePk,patterntablePkV2)
  tablesub<-tablestrack[patterntablePktotal]
  if(!is.null(pattern)) {
    patterntable<-grep(pattern,tablesub ,ignore.case=TRUE)
  }
  lltrack=list()
  for(i in patterntable){
    table.name<-tablesub[i]
    print(table.name)
    tmp<-HistoneOne(gen,chr,start,end,mySession,track.name,table.name)
    if(!is.null(tmp)| length(feature(tmp)) > 0 ){
      print(tmp)
      lltrack=c(lltrack,tmp)
      print(length(lltrack))
    }
    
  }
  
  if(length(lltrack) == 0){
    data_trackfunc <- AnnotationTrack()
    chromosome(data_trackfunc) <- chr
    start(data_trackfunc) <- start
    end(data_trackfunc) <- end
    lltrack=list(data_trackfunc)
  }
  lltrack
}

#-------------------- CREATION track one type of Histone density from UCSC ------------------
HistoneOne<-function(gen,chr,start,end,mySession,track.name="Broad Histone",table.name=NULL){
  if(is.null(chr)){
    stop("Invalid in function HistoneOne :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid in function HistoneOne :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid in function HistoneOne :end null:\n")
  }
  if(is.null(gen)){
    stop("Invalid in function HistoneOne :gen null:\n")
  }
  if(is.null(track.name) & ( !is.na(match(gen, tolower(c("hg19","grch37")))))){
    track.name="Broad Histone"
  }else if(is.null(track.name)){
    stop("Invalid in function HistoneOne :track.namenull:\n")
  }
  if(is.null(table.name)){
    table.name="wgEncodeBroadHistoneGm12878H3k36me3StdPk"
  }else if(is.null(table.name)){
    stop("Invalid in function HistoneOne :gen null:\n")
  }
  mygrange <- GRanges(chr, IRanges(start, end))
  dataUCSC <- getTable(ucscTableQuery (mySession, range=mygrange, 
                                       track=track.name, table=table.name))
  
  
  data_trackfunc <- AnnotationTrack()
  if(nrow(dataUCSC) > 0) {
    data_trackfunc <- AnnotationTrack(chromosome=dataUCSC[,"chrom"],strand ="*",
                                      start=dataUCSC[,"chromStart"],end=dataUCSC[,"chromEnd"],
                                      feature=dataUCSC[,"score"],group=dataUCSC[,"name"],
                                      id=dataUCSC[,"name"], name = "Broad Histone"
                                      ,stacking="dense", col.line = NULL, col = NULL)
    chromosome(data_trackfunc) <- chr
    a<-0:1000
    b<-gray(0:1000 /1000)
    v=(a=b)
    displayPars(data_trackfunc) <- list(v)
  } else {
    data_trackfunc <- AnnotationTrack()
    chromosome(data_trackfunc) <- chr
    start(data_trackfunc) <- start
    end(data_trackfunc) <- end
  }
  
  data_trackfunc
}

#-------------------- CREATION track DNA cluster from UCSC ------------------
DNAseUCSC<-function(gen,chr,start,end,mySession,track.name="DNase Clusters",table.name=NULL){
  if(is.null(chr)){
    stop("Invalid in function DNAseUCS :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid in function DNAseUCS :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid in function DNAseUCS :end null:\n")
  }
  if(is.null(gen)){
    stop("Invalid in function DNAseUCS :gen null:\n")
  }
  if(is.null(track.name) & (!is.na(match(gen, tolower(c("hg19","grch37")))))){
    track.name="DNase Clusters"
  }else if(is.null(track.name)){
    stop("Invalid in function DNAseUCS :track.namenull:\n")
  }
  if(is.null(table.name)){
    table.name="wgEncodeRegDnaseClusteredV3"
  }else if(is.null(table.name)){
    stop("Invalid in function DNAseUCS :gen null:\n")
  }
  mygrange <- GRanges(chr, IRanges(start, end))
  dataUCSC <- getTable(ucscTableQuery (mySession, range=mygrange, 
                                       track=track.name, table=table.name))
  
  data_trackfunc <- AnnotationTrack()
  if(nrow(dataUCSC) > 0) {
    data_trackfunc <- AnnotationTrack(chromosome=dataUCSC[,"chrom"],strand ="*",
                                      start=dataUCSC[,"chromStart"],end=dataUCSC[,"chromEnd"],
                                      feature=dataUCSC[,"score"],group=dataUCSC[,"name"],
                                      id=dataUCSC[,"name"], name = "DNA cluster",
                                      stacking = "dense", col.line = NULL, col = NULL)
    chromosome(data_trackfunc) <- chr
    if(nrow(dataUCSC) > 0) {
      a<-0:1000
      b<-gray(0:1000 /1000)
      v=(a=b)
      displayPars(data_trackfunc) <- list(v)
    } 
    
  } else {
    data_trackfunc <- AnnotationTrack()
    chromosome(data_trackfunc) <- chr
    start(data_trackfunc) <- start
    end(data_trackfunc) <- end
  }
  data_trackfunc
}

#-------------------- CREATION track GC content from UCSC ------------------
gcContent <- function(gen,chr,start,end){
  if(is.null(chr)){
    stop("Invalid in function gcContent :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid in function gcContent :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid in function gcContent :end null:\n")
  }
  if(is.null(gen)){
    stop("Invalid in function gcContent :gen null:\n")
  }
  
  genTrunk <- gsub("\\..*","",gen)
  
  UcscTrack(genome = genTrunk, chromosome = chr, track = "GC Percent", table = "gc5Base", 
            from = start,	to = end, trackType = "DataTrack", start = "start", 
            end = "end", data = "score", type = "hist", window = -1,	windowSize = 1500, 
            fill.histogram = "black",	col.histogram = "red", ylim = c(30, 70), 
            name = "GC Percent", col.line = NULL, col = NULL)
}

#-------------------- CREATION track Known genes from UCSC ------------------
knownGenesUCSC<-function(gen,chr,start,end,showId=TRUE){
  if(is.null(chr)){
    stop("Invalid in function knownGenesUCSC :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid in function knownGenesUCSC :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid in function knownGenesUCSC :end null:\n")
  }
  if(is.null(gen)){
    stop("Invalid in function knownGenesUCSC :gen null:\n")
  }
  
  genTrunk <- gsub("\\..*","",gen)
  
  if(showId){
    UcscTrack(genome = genTrunk, chromosome = chr,track = "knownGene", from = start, to = end, 
              trackType = "GeneRegionTrack", rstarts = "exonStarts", rends = "exonEnds", 
              gene = "name", symbol = "name", transcript = "name", strand = "strand", 
              fill = "#8282d2", name = "UCSC Genes",stacking="squish", group="name",
              fontcolor="black", groupAnnotation = "group", just.group = "above",
              size=2, showId=TRUE,col.line = NULL, col = NULL)
  } else {
    UcscTrack(genome = genTrunk, chromosome = chr,track = "knownGene", from = start, to = end, 
              trackType = "GeneRegionTrack", rstarts = "exonStarts", rends = "exonEnds", 
              gene = "name", symbol = "name", transcript = "name", strand = "strand", 
              fill = "#8282d2", name = "UCSC Genes",stacking="dense",size=2,
              col.line = NULL, col = NULL)
  }
  
}

#-------------------- CREATION track ref Genes from UCSC ------------------
xenorefGenesUCSC<-function(gen,chr,start,end,showId=FALSE){
  if(is.null(chr)){
    stop("Invalid in function refGenesUCSC :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid in function refGenesUCSC :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid in function refGenesUCSC :end null:\n")
  }
  if(is.null(gen)){
    stop("Invalid in function refGenesUCSC :gen null:\n")
  }
  
  genTrunk <- gsub("\\..*","",gen)
  
  if(showId){
    UcscTrack(genome = genTrunk, chromosome = chr, track = "xenoRefGene", 
              from = start, to = end, trackType = "GeneRegionTrack", 
              rstarts = "exonStarts", rends = "exonEnds", gene = "name", 
              symbol = "name2", transcript = "name", strand = "strand", 
              fill = "#8282d2", stacking="squish", name = "Other RefSeq", group="name",
              groupAnnotation = "group", just.group = "above", showId=TRUE,
              col.line = NULL, col = NULL, fontcolor="black")
  }else {
    UcscTrack(genome = genTrunk, chromosome = chr, track = "xenoRefGene", 
              from = start, to = end, trackType = "GeneRegionTrack", 
              rstarts = "exonStarts", rends = "exonEnds", gene = "name", 
              symbol = "name2", transcript = "name", strand = "strand", 
              fill = "#8282d2", stacking="dense", name = "Other RefSeq", 
              col.line = NULL, col = NULL)
  }
  
}

#-------------------- CREATION track CpG Island from UCSC ------------------
cpgIslandsUCSC <-function(gen,chr,start,end){
  if(is.null(chr)){
    stop("Invalid in function cpgIslandsUCSC:chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid in function cpgIslandsUCSC :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid in function cpgIslandsUCSC :end null:\n")
  }
  if(is.null(gen)){
    stop("Invalid in function cpgIslandsUCSC :gen null:\n")
  }
  
  genTrunk <- gsub("\\..*","",gen)
  
  UcscTrack(genome = genTrunk, chromosome = chr, track = "cpgIslandExt", 
            from = start, to = end, trackType = "AnnotationTrack", 
            start = "chromStart", end = "chromEnd", id = "name", shape = "box",
            fill = "#006400", name = "CpG Islands UCSC",stacking="dense",
            col.line = NULL, col = NULL)
}

#-------------------- CREATION track SNPs from UCSC ------------------
snpLocationsUCSC <-function(gen,chr,start,end,track){ 
  if(is.null(chr)){
    stop("Invalid in function snpLocationsUCSC:chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid in function snpLocationsUCSC :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid in function snpLocationsUCSC :end null:\n")
  }
  if(is.null(gen)){
    stop("Invalid in function snpLocationsUCSC :gen null:\n")
  }
  if(is.null(track) & (gen== "hg19" | gen == "grch37")){
    track="snp142"
  }else if(is.null(track)){
    stop("Invalid in function snpLocationsUCSC :track null:\n")
  }
  
  genTrunk <- gsub("\\..*","",gen)
  
  UcscTrack(genome = genTrunk, chromosome = chr, track = track, from = start, to = end, 
            trackType = "AnnotationTrack", start = "chromStart", end = "chromEnd", 
            id = "name", feature = "func", strand = "strand", shape = "box", 
            stacking="dense", fill = "black", name = "SNPs UCSC",  
            col.line = NULL, col = NULL)
}

#-------------------- CREATION track Regulation from ENSEMBL ------------------
regulationBiomart <- function(gen, chr, start, end) {
  if(is.null(gen)){
    stop("Invalid in function regulationBiomart :gen null:\n")
  }
  if(is.null(chr)){
    stop("Invalid in function regulationBiomart :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid in function regulationBiomart :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid in function regulationBiomart :end null:\n")
  }
  
  genTrunk <- gsub("\\..*","",gen)
  chrEnsembl=chrUCSC2ENSEMBL(chr)
  martfunc=NULL
  dataset="hsapiens_motif_feature"
  
  if(!is.na(match(tolower(genTrunk), c("grch37","hg19")))){
    martfunc=useMart(host='grch37.ensembl.org', biomart='ENSEMBL_MART_FUNCGEN',
                     dataset='hsapiens_regulatory_feature')
  } else {
    
    martfunc <- useMart('functional_genomics',dataset='hsapiens_regulatory_feature')
  }
  ensfunc <- getBM(c("chromosome_name","chromosome_start","chromosome_end",
                     "feature_type_name"),
                   filters = c("chromosome_name","start","end"),
                   values = list(chrEnsembl, start, end), mart=martfunc) 
  
  data_trackfunc <- AnnotationTrack()
  if(nrow(ensfunc) > 0) {
    data_trackfunc <- AnnotationTrack(chromosome=chrEnsembl,strand ="*",start=ensfunc[,2],
                                      end=ensfunc[,3],
                                      feature=ensfunc[,4],group=ensfunc[,1],id=ensfunc[,1], 
                                      name = "Regulation ENSEMBL",stacking="dense", 
                                      col.line = NULL, col = NULL)
    displayPars(data_trackfunc) <- list(
      "Promoter Associated"="darkolivegreen",
      "CTCF Binding Site" = "cadetblue1",
      "Gene Associated" = "coral",
      "Non-Gene Associated" = "darkgoldenrod1",
      "Predicted Transcribed Region" = "greenyellow",
      "PolIII Transcription Associated" = "purple",
      "Enhancer" = "gold",
      "Transcription Factor Binding Site" = "darkorchid1",
      "Predicted Weak enhancer/Cis-reg element" = "yellow",
      "Heterochromatin" = "wheat4",
      "Open Chromatin" = "snow3",
      "Promoter Flank" = "tomato",
      "Repressed/Low Activity" ="snow4",
      "Unclassified" = "aquamarine")
    
  } else {
    data_trackfunc <- AnnotationTrack()
    chromosome(data_trackfunc) <- chr
    start(data_trackfunc) <- start
    end(data_trackfunc) <- end
  }
  
  data_trackfunc
}

#-------------------- CREATION track Short Variation from ENSEMBL ------------------
snpBiomart <- function(chr, start, end, dataset, showId=FALSE, title=NULL) {
  if(is.null(chr)){
    stop("Invalid in function snpBiomart :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid in function snpBiomart :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid in function snpBiomart :end null:\n")
  }
  if(is.null(dataset)){
    stop("Invalid in function snpBiomart :dataset null:\n")
  }
  if(is.null(title)){
    title="Short Variation"
  }
  
  chrEnsembl=chrUCSC2ENSEMBL(chr)
  martsnp=useMart("snp",dataset=dataset)
  ens_snp <- getBM(c("refsnp_id","chrom_start","chrom_strand","chr_name"),
                   filters = c("chr_name","start","end"),
                   values = list(chrEnsembl, start, end), mart=martsnp) 
  
  data_tracksnp <- AnnotationTrack()
  if(nrow(ens_snp) > 0) {
    data_tracksnp <- AnnotationTrack(chromosome=ens_snp[,4],strand =ens_snp[,3],start=ens_snp[,2],
                                     end=ens_snp[,2],feature="snp",group=ens_snp[,1],
                                     id=ens_snp[,1], name = title,stacking="dense",
                                     showId=showId,  col.line = NULL, col = NULL)
    displayPars(data_tracksnp) <- list(snp="red",
                                       insertion ="blueviolet",
                                       deletion = "orange",
                                       indel="darkgoldenrod1",
                                       substitution="dodgerblue2")
    
  } else {
    data_tracksnp <- AnnotationTrack()
    chromosome(data_tracksnp) <- chr
    start(data_tracksnp) <- start
    end(data_tracksnp) <- end
  }
  
  data_tracksnp
  
}

#-------------------- CREATION track Structural Variation from ENSEMBL ------------------
structureBiomart <- function(chr, start, end, strand, dataset,showId=FALSE,title=NULL) {
  if(is.null(chr)){
    stop("Invalid in function snpBiomart :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid in function snpBiomart :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid in function snpBiomart :end null:\n")
  }
  if(is.null(dataset)){
    stop("Invalid in function snpBiomart :dataset null:\n")
  }
  if(is.null(title)){
    title="Structural Variation"
  }
  chrEnsembl=chrUCSC2ENSEMBL(chr)
  martstruct=useMart("snp",dataset=dataset)
  ens <- getBM(c("sv_accession","chrom_start","chrom_end","seq_region_strand","chr_name",
                 "sv_variant_type","dgva_study_accession"),
               filters = c("chr_name","start","end"),
               values = list(chrEnsembl, start, end), mart=martstruct) 
  
  data_track <- AnnotationTrack()
  if(nrow(ens) > 0) {
    data_track <- AnnotationTrack(chromosome=chr,strand ="*",start=ens[,2],end=ens[,3],
                                  feature=ens[,6],group=ens[,1],id=ens[,1], 
                                  name = "Structural variation",stacking="squish",
                                  showId=showId,  col.line = NULL, col = NULL)
    displayPars(data_track) <- list(copy_number_variation="cornsilk",
                                    inversion="darkolivegreen",
                                    translocation="cyan",
                                    sequence_alteration="coral",
                                    snp="red",
                                    insertion ="blueviolet",
                                    deletion = "orange",
                                    indel="darkgoldenrod1",
                                    substitution="dodgerblue2")
  } else {
    data_track <- AnnotationTrack()
    chromosome(data_track) <- chr
    start(data_track) <- start
    end(data_track) <- end
  }
  
  data_track
}

#-------------------- CREATION track ClinVar Variants  Main from UCSC ------------------
ClinVarMainTrack <-function(gen,chr,start,end,showId=FALSE){ 
  if(is.null(chr)){
    stop("Invalid in function ClinVarUCSC:chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid in function ClinVarUCSC :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid in function ClinVarUCSC :end null:\n")
  }
  if(is.null(gen)){
    stop("Invalid in function ClinVarUCSC :gen null:\n")
  }
  
  genTrunk <- gsub("\\..*","",gen)
  
  if(showId){
    UcscTrack(genome = genTrunk, chromosome = chr, from = start, to = end,
              track="ClinVar Variants", table="clinvarMain", 
              trackType = "AnnotationTrack", start = "chromStart", end = "chromEnd", 
              id = "name", feature = "func", strand = "*", shape = "box", 
              stacking="squish", fill = "black", name = "ClinVar Variants", group="name",
              groupAnnotation = "group", just.group = "above",  showId=TRUE,
              col.line = NULL, col = NULL)
  } else {
    UcscTrack(genome = genTrunk, chromosome = chr, from = start, to = end,
              track="ClinVar Variants", table="clinvarMain", 
              trackType = "AnnotationTrack", start = "chromStart", end = "chromEnd", 
              id = "name", feature = "func", strand = "*", shape = "box", 
              stacking="dense", fill = "black", name = "ClinVar Variants", 
              col.line = NULL, col = NULL)
  }
  
}


#-------------------- CREATION track ClinVar Variants CNV from UCSC ------------------
ClinVarCnvTrack <-function(gen,chr,start,end,showId=FALSE){ 
  if(is.null(chr)){
    stop("Invalid in function ClinVarCnv:chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid in function ClinVarCnv :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid in function ClinVarCnv :end null:\n")
  }
  if(is.null(gen)){
    stop("Invalid in function ClinVarCnv :gen null:\n")
  }
  
  genTrunk <- gsub("\\..*","",gen)
  
  if(showId){
    UcscTrack(genome = genTrunk, chromosome = chr, from = start, to = end,
              track="ClinVar Variants", table="clinvarCnv", 
              trackType = "AnnotationTrack", start = "chromStart", end = "chromEnd", 
              id = "name", feature = "func", strand = "*", shape = "box", 
              stacking="squish", fill = "black", name = "ClinVar Variants", group="name",
              groupAnnotation = "group", just.group = "above",  showId=TRUE,
              col.line = NULL, col = NULL)
  } else {
    UcscTrack(genome = genTrunk, chromosome = chr, from = start, to = end,
              track="ClinVar Variants", table="clinvarCnv", 
              trackType = "AnnotationTrack", start = "chromStart", end = "chromEnd", 
              id = "name", feature = "func", strand = "*", shape = "box", 
              stacking="dense", fill = "black", name = "ClinVar Variants", 
              col.line = NULL, col = NULL)
  }
  
}

#-------------------- CREATION track Coriell CNV from UCSC ------------------
CoreillCNVTrack <-function(gen,chr,start,end,showId=FALSE){ 
  if(is.null(chr)){
    stop("Invalid in function CoreillCNV:chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid in function CoreillCNV :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid in function CoreillCNV :end null:\n")
  }
  if(is.null(gen)){
    stop("Invalid in function CoreillCNV :gen null:\n")
  }
  
  genTrunk <- gsub("\\..*","",gen)
  
  if(showId){
    UcscTrack(genome = genTrunk, chromosome = chr, from = start, to = end,
              track="Coriell CNVs", table="coriellDelDup", 
              trackType = "AnnotationTrack", start = "chromStart", end = "chromEnd", 
              id = "name", feature = "func", strand = "*", shape = "box", 
              stacking="squish", fill = "blue", name = "Coriell CNVs", group="name",
              groupAnnotation = "group", just.group = "above",  showId=TRUE,
              col.line = NULL, col = NULL)
    
  }else {
    UcscTrack(genome = genTrunk, chromosome = chr, from = start, to = end,
              track="Coriell CNVs", table="coriellDelDup", 
              trackType = "AnnotationTrack", start = "chromStart", end = "chromEnd", 
              id = "name", feature = "func", strand = "*", shape = "box", 
              stacking="dense", fill = "blue", name = "Coriell CNVs", 
              col.line = NULL, col = NULL)
  }
  
}

#-------------------- CREATION track COSMIC from UCSC ------------------
COSMICTrack <-function(gen,chr,start,end,showId=FALSE){ 
  if(is.null(chr)){
    stop("Invalid in function COSMICUCSC:chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid in function COSMICUCSC :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid in function COSMICUCSC :end null:\n")
  }
  if(is.null(gen)){
    stop("Invalid in function COSMICUCSC :gen null:\n")
  }
  
  genTrunk <- gsub("\\..*","",gen)
  
  if(showId) {
    UcscTrack(genome = genTrunk, chromosome = chr, from = start, to = end,
              track="COSMIC", table="cosmic", 
              trackType = "AnnotationTrack", start = "chromStart", end = "chromEnd", 
              id = "name", feature = "func", strand = "*", shape = "box", 
              stacking="squish", fill = "firebrick1", name = "COSMIC", group="name",
              groupAnnotation = "group", just.group = "above",  showId=TRUE,
              col.line = NULL, col = NULL)
  } else {
    UcscTrack(genome = genTrunk, chromosome = chr, from = start, to = end,
              track="COSMIC", table="cosmic", 
              trackType = "AnnotationTrack", start = "chromStart", end = "chromEnd", 
              id = "name", feature = "func", strand = "*", shape = "box", 
              stacking="dense", fill = "firebrick1", name = "COSMIC", 
              col.line = NULL, col = NULL)
  }
  
}

#-------------------- CREATION track GAD from UCSC ------------------
GADTrack <-function(gen,chr,start,end,showId=FALSE){ 
  if(is.null(chr)){
    stop("Invalid in function GAD:chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid in function GAD:start null:\n")
  }
  if(is.null(end)){
    stop("Invalid in function GAD:end null:\n")
  }
  if(is.null(gen)){
    stop("Invalid in function GAD:gen null:\n")
  }
  
  genTrunk <- gsub("\\..*","",gen)
  
  if(showId){
    UcscTrack(genome = genTrunk, chromosome = chr, from = start, to = end,
              track="GAD View", table="gad", 
              trackType = "AnnotationTrack", start = "chromStart", end = "chromEnd", 
              id = "name", group ="name", feature = "func", strand = "*", shape = "box", 
              stacking="squish", fill = "darkslategray1", name = "GAD",groupAnnotation = "group",
              just.group = "above", showId=TRUE, col.line = NULL, col = NULL)
  } else {
    UcscTrack(genome = genTrunk, chromosome = chr, from = start, to = end,
              track="GAD View", table="gad", 
              trackType = "AnnotationTrack", start = "chromStart", end = "chromEnd", 
              id = "name", feature = "func", strand = "*", shape = "box", 
              stacking="dense", fill = "darkslategray1", name = "GAD", 
              col.line = NULL, col = NULL)
  }
  
}

#-------------------- CREATION track raw GWAS Catalog from UCSC ------------------
GWASTrack <-function(gen,chr,start,end,showId=FALSE){ 
  if(is.null(chr)){
    stop("Invalid in function GWAS:chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid in function GWAS:start null:\n")
  }
  if(is.null(end)){
    stop("Invalid in function GWAS:end null:\n")
  }
  if(is.null(gen)){
    stop("Invalid in function GWAS:gen null:\n")
  }
  
  genTrunk <- gsub("\\..*","",gen)
  
  if(showId){
    UcscTrack(genome = genTrunk, chromosome = chr, from = start, to = end,
              track="GWAS Catalog", table="gwasCatalog", 
              trackType = "AnnotationTrack", start = "chromStart", end = "chromEnd", 
              id = "name", group="name", feature = "func", strand = "*", shape = "box", 
              stacking="squish", fill = "black", name = "GWAS Catalog",groupAnnotation = "group",
              just.group = "above",  col.line = NULL, col = NULL, showId=TRUE)
  }else {
    UcscTrack(genome = genTrunk, chromosome = chr, from = start, to = end,
              track="GWAS Catalog", table="gwasCatalog", 
              trackType = "AnnotationTrack", start = "chromStart", end = "chromEnd", 
              id = "name", feature = "func", strand = "*", shape = "box", 
              stacking="dense", fill = "black", name = "GWAS Catalog",  col.line = NULL, col = NULL)
  }
  
}

#-------------------- CREATION track GeneReviews from UCSC ------------------
GeneReviewsTrack <-function(gen,chr,start,end,showId=FALSE){ 
  if(is.null(chr)){
    stop("Invalid in function GeneReviews:chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid in function GeneReviews:start null:\n")
  }
  if(is.null(end)){
    stop("Invalid in function GeneReviews:end null:\n")
  }
  if(is.null(gen)){
    stop("Invalid in function GeneReviews:gen null:\n")
  }
  
  genTrunk <- gsub("\\..*","",gen)
  
  if(showId){
    UcscTrack(genome = genTrunk, chromosome = chr, from = start, to = end,
              track="GeneReviews", table="geneReviews", 
              trackType = "AnnotationTrack", start = "chromStart", end = "chromEnd", 
              id = "name", group="name", feature = "func", strand = "*", shape = "box", 
              stacking="squish", fill = "red", name = "GeneReviews",
              groupAnnotation = "group",just.group = "above", showId=TRUE)
  }else {
    UcscTrack(genome = genTrunk, chromosome = chr, from = start, to = end,
              track="GeneReviews", table="geneReviews", 
              trackType = "AnnotationTrack", start = "chromStart", end = "chromEnd", 
              id = "name", feature = "func", strand = "*", shape = "box", 
              stacking="dense", fill = "red", name = "GeneReviews")
  }
  
}

#-------------------- CREATION track ISCA from UCSC ------------------
ISCATrack <-function(gen,chr,start,end,mySession,table.name,showId=FALSE){ 
  if(is.null(chr)){
    stop("Invalid in function ISCA:chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid in function ISCA:start null:\n")
  }
  if(is.null(end)){
    stop("Invalid in function ISCA:end null:\n")
  }
  if(is.null(gen)){
    stop("Invalid in function ISCA:gen null:\n")
  }
  
  genTrunk <- gsub("\\..*","",gen)
  
  if((is.null(table.name) | ! exists(table.name))& (genTrunk == "hg19" | genTrunk == "grch37")){
    table.name="iscaPathogenic"
  }else if(is.null(table.name)){
    stop("Invalid in function ISCA:table null (possible table : iscaBenign
         , iscaCuratedBenign, iscaCuratedPathogenic, iscaLikelyBenign, 
         iscaLikelyPathogenic, iscaPathGainCum, iscaPathLossCum, 
         iscaPathogenic, iscaUncertain )\n")
  }
  
  #  UcscTrack(genome = gen, chromosome = chr, from = start, to = end,
  #           track="ISCA", table= table.name, 
  #          trackType = "AnnotationTrack", start = "chromStart", end = "chromEnd", 
  #         id = "name", feature = "func", strand = "*", shape = "box", 
  #        stacking="squish", fill = "deeppink3", name = "ISCA",showId=showId)
  
  mygrange <- GRanges(chr, IRanges(start, end))
  dataUCSC <- getTable(ucscTableQuery (mySession, range=mygrange, 
                                       track="ISCA", table=table.name))
  
  data_trackfunc <- AnnotationTrack()
  if(nrow(dataUCSC) > 0) {
    data_trackfunc <- AnnotationTrack(chromosome=dataUCSC[,"chrom"],strand ="*",
                                      start=dataUCSC[,"chromStart"],end=dataUCSC[,"chromEnd"],
                                      feature=dataUCSC[,"score"],group=dataUCSC[,"name"],
                                      id=dataUCSC[,"name"], name = "ISCA",
                                      stacking = "squish", showId=TRUE,
                                      groupAnnotation = "group",just.group = "above",
                                      col.line = NULL, col = NULL)
    chromosome(data_trackfunc) <- chr
    
    
    if(table.name == "iscaPathogenic") {
      displayPars(data_trackfunc) <- list(fill="purple")
    }
    if(table.name == "iscaPathGainCum") {
      displayPars(data_trackfunc) <- list(fill="red")
    }
    if(table.name == "iscaPathLossCum") {
      displayPars(data_trackfunc) <- list(fill="blue")
    }
    if(table.name == "iscaCuratedPathogenic") {
      displayPars(data_trackfunc) <- list(fill="purple")
    }
    if(table.name == "iscaLikelyPathogenic") {
      displayPars(data_trackfunc) <- list(fill="lightpurple")
    }
    if(table.name == "iscaUncertain") {
      displayPars(data_trackfunc) <- list(fill="lightgrey")
    }
    if(table.name == "iscaBenign") {
      displayPars(data_trackfunc) <- list(fill="black")
    }
    if(table.name == "iscaCuratedBenign") {
      displayPars(data_trackfunc) <- list(fill="black")
    }
    if(table.name == "iscaLikelyBenign") {
      displayPars(data_trackfunc) <- list(fill="black")
    }
    
  } else{
    data_trackfunc <- AnnotationTrack()
    chromosome(data_trackfunc) <- chr
    start(data_trackfunc) <- start
    end(data_trackfunc) <- end
  }
  data_trackfunc
}

#-------------------- CREATION track repeatMasker from UCSC ------------------
RepeatMaskerTrack <-function(gen,chr,start,end,showId=FALSE){ 
  if(is.null(chr)){
    stop("Invalid in function repeatMasker:chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid in function repeatMasker:start null:\n")
  }
  if(is.null(end)){
    stop("Invalid in function repeatMasker:end null:\n")
  }
  if(is.null(gen)){
    stop("Invalid in function repeatMasker:gen null:\n")
  }
  
  genTrunk <- gsub("\\..*","",gen)
  
  UcscTrack(genome = genTrunk, chromosome = chr, from = start, to = end,
            track="RepeatMasker", table="rmsk", 
            trackType = "AnnotationTrack", start = "genoStart", end = "genoEnd", 
            id = "repName", group="repName", feature = "repClass", strand = "*", 
            shape = "box", stacking="full", fill = "grey", name = "RepeatMasker",
            groupAnnotation = "group",just.group = "above", 
            showId=showId, col.line = NULL, col = NULL)
}

## New database in Regulation ENSEMBL
#-------------------- CREATION track for all other regulatory regions from ENSEMBL or a list of them ------------
OtherRegulatoryRegions <- function (gen, chr, start, end, featureDisplay = "all",datasetEnsembl = NULL) {
  if(is.null(gen)){
    stop("Invalid function OtherRegulatoryRegions :gen null:\n")
  }
  if(is.null(chr)){
    stop("Invalid function OtherRegulatoryRegions :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid function OtherRegulatoryRegions :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid function OtherRegulatoryRegions :end null:\n")
  }
  
  options(ucscChromosomeNames=FALSE)
  chrEnsembl <- chrUCSC2ENSEMBL(chr)
  
  biomTrack <- NULL
  martfunc <- NULL
  biomTrackDislay <- NULL
  
  genTrunk <- gsub("\\..*","",gen)
  
  if(!is.na(match(tolower(genTrunk), c("grch37","hg19")))){
    martfunc <- useMart(host='grch37.ensembl.org', biomart='ENSEMBL_MART_FUNCGEN',
                        dataset='hsapiens_external_feature')
    
  } else if(!is.na(match(tolower(genTrunk), c("grch38","hg38")))) {
    martfunc <- useMart("regulation", dataset = "hsapiens_external_feature")
    
  } else if(!is.null(datasetEnsembl) && !is.na(datasetEnsembl)) {
    martfunc <- useMart("regulation", dataset = datasetEnsembl)
    
  } else {
    stop("Invalid function OtherRegulatoryRegions :genome not recognised:\n")
  }
  
  biomTrack <- getBM(c("chromosome_name", "chromosome_start","chromosome_end",
                       "feature_type_class"), 
                     filters = c("chromosome_name","start","end"),
                     values = list(chrEnsembl, start, end), mart = martfunc)
  
  biomTrackDislay <- biomTrack
  
  if( !("all" %in% featureDisplay) ) {
    biomTrackDislay <- biomTrack[which(biomTrack$feature_type_class %in% featureDisplay),]
  }
  
  data_trackfunc <- AnnotationTrack(genome = genTrunk,chromosome=chrEnsembl,strand ="*",start=biomTrackDislay[,2],
                                    end=biomTrackDislay[,3],
                                    feature=biomTrackDislay[,4],group=biomTrackDislay[,1],
                                    id=biomTrackDislay[,1], name = "Other Regulatory Regions ENSEMBL",stacking="dense", 
                                    col.line = "black", col = NULL, collapse= FALSE)
  displayPars(data_trackfunc) <- list("Enhancer" = "#e41a1c", "Transcription Start Site" = "#4daf4a"
  )
  
  data_trackfunc
  
} 

#-------------------- CREATION track for all Regulatory features from ENSEMBL or a list of them ----------------
RegulatoryFeaturesBiomart  <- function (gen, chr, start, end, featureDisplay = "all",datasetEnsembl = NULL) {
  if(is.null(gen)){
    stop("Invalid function RegulatoryFeaturesBiomart :gen null:\n")
  }
  if(is.null(chr)){
    stop("Invalid function RegulatoryFeaturesBiomart :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid function RegulatoryFeaturesBiomart :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid function RegulatoryFeaturesBiomart :end null:\n")
  }
  
  options(ucscChromosomeNames=FALSE)
  chrEnsembl <- chrUCSC2ENSEMBL(chr)
  
  biomTrack <- NULL
  martfunc <- NULL
  biomTrackDislay <- NULL
  
  genTrunk <- gsub("\\..*","",gen)
  
  if(!is.na(match(tolower(genTrunk), c("grch37","hg19")))){
    martfunc <- useMart(host='grch37.ensembl.org', biomart='ENSEMBL_MART_FUNCGEN',
                        dataset='hsapiens_regulatory_feature')
    
  } else if(!is.na(match(tolower(genTrunk), c("grch38","hg38")))) {
    martfunc <- useMart("regulation", dataset = "hsapiens_regulatory_feature")
    
  } else if(!is.null(datasetEnsembl) && !is.na(datasetEnsembl)) {
    martfunc <- useMart("regulation", dataset = datasetEnsembl)
    
  } else {
    stop("Invalid function RegulatoryFeaturesBiomart :genome not recognised:\n")
  }
  
  
  
  biomTrack <- getBM(c("chromosome_name", "chromosome_start","chromosome_end",
                       "feature_type_name"), 
                     filters = c("chromosome_name","start","end"),
                     values = list(chrEnsembl, start, end), mart = martfunc)
  
  biomTrackDislay <- biomTrack
  
  if( !("all" %in% featureDisplay) ) {
    biomTrackDislay <- biomTrack[which(biomTrack$feature_type_class %in% featureDisplay),]
  }
  
  data_trackfunc <- AnnotationTrack(genome = genTrunk,chromosome=chrEnsembl,strand ="*",start=biomTrackDislay[,2],
                                    end=biomTrackDislay[,3],
                                    feature=biomTrackDislay[,4],group=biomTrackDislay[,1],
                                    id=biomTrackDislay[,1], name = "Other Regulatory Regions ENSEMBL",stacking="dense", 
                                    col.line = "black", col = NULL, collapse= FALSE)
  displayPars(data_trackfunc) <- list("Promoter" = "#1b9e77", "TF binding site" = "#d95f02", 
                                      "Open chromatin" = "#7570b3", "Promoter Flanking Region" =  "#e7298a", 
                                      "CTCF Binding Site" = "#66a61e", "Enhancer" = "#e6ab02")
  data_trackfunc
  
} 

#-------------------- CREATION track for all Regulatory segments from ENSEMBL or a list of them ----------------
RegulatorySegmentsBiomart  <- function (gen, chr, start, end, featureDisplay = "all",datasetEnsembl = NULL) {
  if(is.null(gen)){
    stop("Invalid function RegulatorySegmentsBiomart :gen null:\n")
  }
  if(is.null(chr)){
    stop("Invalid function RegulatorySegmentsBiomart :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid function RegulatorySegmentsBiomart :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid function RegulatorySegmentsBiomart :end null:\n")
  }
  
  options(ucscChromosomeNames=FALSE)
  chrEnsembl <- chrUCSC2ENSEMBL(chr)
  
  biomTrack <- NULL
  martfunc <- NULL
  biomTrackDislay <- NULL
  
  genTrunk <- gsub("\\..*","",gen)
  
  if(!is.na(match(tolower(genTrunk), c("grch37","hg19")))){
    martfunc <- useMart(host='grch37.ensembl.org', biomart='ENSEMBL_MART_FUNCGEN',
                        dataset='hsapiens_segmentation_feature')
    
  } else if(!is.na(match(tolower(genTrunk), c("grch38","hg38")))) {
    martfunc <- useMart("regulation", dataset = "hsapiens_segmentation_feature")
    
  } else if(!is.null(datasetEnsembl) && !is.na(datasetEnsembl)) {
    martfunc <- useMart("regulation", dataset = datasetEnsembl)
    
  } else {
    stop("Invalid function RegulatorySegmentsBiomart: genome not recognised:\n")
  }
  
  
  
  biomTrack <- getBM(c("chromosome_name", "chromosome_start","chromosome_end",
                       "feature_type_name"), 
                     filters = c("chromosome_name","start","end"),
                     values = list(chrEnsembl, start, end), mart = martfunc)
  
  biomTrackDislay <- biomTrack
  
  if( !("all" %in% featureDisplay) ) {
    biomTrackDislay <- biomTrack[which(biomTrack$feature_type_class %in% featureDisplay),]
  }
  
  data_trackfunc <- AnnotationTrack(genome = genTrunk,chromosome=chrEnsembl,strand ="*",start=biomTrackDislay[,2],
                                    end=biomTrackDislay[,3],
                                    feature=biomTrackDislay[,4],group=biomTrackDislay[,1],
                                    id=biomTrackDislay[,1], name = "Other Regulatory Regions ENSEMBL",stacking="dense", 
                                    col.line = "black", col = NULL, collapse= FALSE)
  displayPars(data_trackfunc) <- list("Predicted Promoter with TSS" = "#7fc97f","CTCF enriched" = "#beaed4",
                                      "Predicted Promoter Flank" = "#fdc086", "Predicted Enhancer" = "#ffff99",
                                      "Predicted Repressed" = "#386cb0", "Predicted low activity" = "#f0027f",
                                      "Predicted Transcribed Region" = "#bf5b17", "Predicted heterochomatin" = "#666666")
  
  data_trackfunc
} 

#-------------------- CREATION track for all binding motifs from ENSEMBL or a list of them ----------------

BindingMotifsBiomart <- function (gen, chr, start, end, featureDisplay = "all",datasetEnsembl = NULL) {
  if(is.null(gen)){
    stop("Invalid function BindingMotifsBiomart :gen null:\n")
  }
  if(is.null(chr)){
    stop("Invalid function BindingMotifsBiomart :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid function BindingMotifsBiomart :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid function BindingMotifsBiomart :end null:\n")
  }
  
  options(ucscChromosomeNames=FALSE)
  chrEnsembl <- chrUCSC2ENSEMBL(chr)
  
  biomTrack <- NULL
  martfunc <- NULL
  biomTrackDislay <- NULL
  
  genTrunk <- gsub("\\..*","",gen)
  
  if(!is.na(match(tolower(genTrunk), c("grch37","hg19")))){
    martfunc <- useMart(host='grch37.ensembl.org', biomart='ENSEMBL_MART_FUNCGEN',
                        dataset='hsapiens_motif_feature')
    
  } else if(!is.na(match(tolower(genTrunk), c("grch38","hg38")))) {
    martfunc <- useMart("regulation", dataset = "hsapiens_motif_feature")
    
  } else if(!is.null(datasetEnsembl) && !is.na(datasetEnsembl)) {
    martfunc <- useMart("regulation", dataset = datasetEnsembl)
    
  } else {
    stop("Invalid function BindingMotifsBiomart :genome not recognised:\n")
  }
  
  
  
  biomTrack <- getBM(c("chromosome_name", "chromosome_start","chromosome_end",
                       "feature_type_name"), 
                     filters = c("chromosome_name","start","end"),
                     values = list(chrEnsembl, start, end), mart = martfunc)
  
  biomTrackDislay <- biomTrack
  
  if( !("all" %in% featureDisplay) ) {
    biomTrackDislay <- biomTrack[which(biomTrack$feature_type_name %in% featureDisplay),]
  }
  
  data_trackfunc <- AnnotationTrack(genome = genTrunk,chromosome=chrEnsembl,strand ="*",start=biomTrackDislay[,2],
                                    end=biomTrackDislay[,3],
                                    feature=biomTrackDislay[,4],group=biomTrackDislay[,1],
                                    id=biomTrackDislay[,1], name = "Regulation ENSEMBL",stacking="dense", 
                                    col.line = "black", col = NULL, collapse= FALSE)
  displayPars(data_trackfunc) <- list(
    "Egr1" = "#a6cee3", "CTCF" = "#1f78b4", "Cjun" = "#b2df8a", "USF1" = "#33a02c", "PU1" = "#fb9a99",
    "Gabp" = "#e31a1c", "JUN::FOS" = "#fdbf6f", "Jund" = "#ff7f00", "Znf263" = "#cab2d6", 
    "FOXA1" = "#6a3d9a", "E2F4" = "#ffff99", "SP1" = "#b15928", "Yy1" = "#62725b","Srf" = "#383838",
    "HNF4A" = "#AFE6DC", "Nrsf" = "#3B3E19", "FOSL2" = "#F05868", "MYC::MAX" = "#8B2323", "Max" = "#5BBAAE",
    "E2F6" = "#FFD3D7", "EBF1" = "#FFCA08", "Nrf1" = "#A0A0A0", "MEF2A" = "#3088F0", "ELF1" = "#0BFCFF",
    "Cfos" = "#8E6363", "ZBTB33" = "#2F6568", "Cmyc" = "#FFF803", "FOSL1" = "#8FB247", "Tcf12" = "#FC03FF", 
    "FOXA2" = "#03FFAC", "CTCFL" = "#FF5703", "SP2" = "#BB9C36", "ZEB1" = "#A036BB", "Pax5" = "#36BB98",
    "NFKB" = "#C4FF00", "RXRA" = "#18044F", "HNF4G" = "#6F41F0", "IRF4" = "#1F9433", "Gata2" = "#86941F",
    "Tr4" = "#1A5E45", "PPARG::RXRA" = "#5E1A52", "Junb" = "#8A0A66", "RXR::RAR_DR5" = "#0A328A", 
    "Tal1::Gata1" = "#4A8A0A", "POU2F2" = "#8A520A", "RXRA::VDR" = "#94DB48", "Gata1" = "#74C2D6",
    "ETS1" = "#57C716", "Nr1h3::Rxra" = "#2E16B5", "EcR::usp" = "#B51660", "THAP1" = "#9BB516",
    "BHLHE40" = "#157BAB", "MEF2C" = "#5C728C", "Pbx3" = "#E31D42","E2F1" = "#CCCC00", "SRebp1" = "#A33333",
    "SRebp2" = "#B2B2F0")
  
  data_trackfunc
} 

#-------------------- CREATION track for all Segmental duplications from UCSC  ----------------
SegmentalDupsUCSC <- function(gen, chr, start, end, showId=FALSE){
  if(is.null(gen)){
    stop("Invalid function SegmentalDupsUCSC :gen null:\n")
  }
  if(is.null(chr)){
    stop("Invalid function SegmentalDupsUCSC :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid function SegmentalDupsUCSC :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid function SegmentalDupsUCSC :end null:\n")
  }
  
  Dupregions <- UcscTrack(genome = gen, chromosome = chr, from = start, to = end, 
                          track = "Segmental Dups", table = "genomicSuperDups", 
                          trackType = "AnnotationTrack", start = "chromStart", 
                          end="chromEnd", id = "name", name = "Segmental Dups UCSC",
                          showID=showId)
  
  Dupregions
  
}

#-------------------- CREATION track for CG content  ----------------
# CGcontent <-  function(bsgenome, chr, start, end, tilewidth)
# {
# #  dna <- bsgenome[[chr]]
#   dna <- getSeq(bsgenome, chr, start, end)
#   ## CpG on the plus and minus strand (?)
#   islands <- matchPDict(DNAStringSet(c("GC", "CG")), dna)
#   cvg <- coverage(islands)    # CpG island coverage
#   
#   if(typewidth == "tiles") {
#     tiles <- tileGenome(seqlengths(bsgenome)[chr], tilewidth=tilewidth,
#                          cut.last.tile.in.chrom=TRUE)
#     ## Average coverage in each tile
#     ## Divide by 2 so each CpG counts only once
#     v <- Views(cvg, ranges(tiles))
#     tiles$CpG <- viewSums(v) / width(v) / 2
#     tiles
#   } else {
#     region <- diff(cumsum(cvg), lag=slidewidth) / slidewidth / 2
#   }
#   
# }
