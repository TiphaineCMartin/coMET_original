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
genesName_ENSEMBL<-function(gen,chr,start,end,dataset){
  if(is.null(gen)){
    stop("Invalid in function genesNameENSEMBL :gen null:\n")
  }
  if(is.null(chr)){
    stop("Invalid in function genesNameENSEMBL :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid in function genesNameENSEMBL :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid in function genesNameENSEMBL :end null:\n")
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
genes_ENSEMBL<-function(gen,chr,start,end,showId=FALSE){
  if(is.null(gen)){
    stop("Invalid in function genesENSEMBL :gen null:\n")
  }
  if(is.null(chr)){
    stop("Invalid in function genesENSEMBL :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid in function genesENSEMBL :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid in function genesENSEMBL :end null:\n")
  }
  
  #cat("data",gen,"\t",chr,"\t",start,"\t",end,"\n")
  
  genTrunk <- gsub("\\..*","",gen)
  
  chrEnsembl=chrUCSC2ENSEMBL(chr)
  
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
                                        chromosome = chrEnsembl, start = start, 
                                        end = end,  name = "genes ENSEMBL",
                                        groupAnnotation = "group",
                                        just.group = "above",
                                        fontcolor="black",showId=showId,size=2,
                                        col.line = NULL, col = NULL)
    
    
  } else {
    martENSEMBL=useMart("ensembl",dataset='hsapiens_gene_ensembl')
    biomTrack <- BiomartGeneRegionTrack(genome = genTrunk, biomart=martENSEMBL,
                                        chromosome = chrEnsembl, start = start, 
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
interestGenes_ENSEMBL<-function(gen,chr,start,end,interestfeatures,interestcolor,showId=FALSE){
  if(is.null(gen)){
    stop("Invalid in function interestGenesENSEMBL :gen null:\n")
  }
  if(is.null(chr)){
    stop("Invalid in function interestGenesENSEMBL :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid in function interestGenesENSEMBL :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid in function interestGenesENSEMBL :end null:\n")
  }
  if(is.null(interestfeatures)){
    stop("Invalid in function interestGenesENSEMBL: interestfeatures null:\n")
  }
  if(is.null(interestcolor)){
    stop("Invalid in function interestGenesENSEMBL: interestcolor null:\n")
  }
  #cat("data",gen,"\t",chr,"\t",start,"\t",end,"\n")
  
  genTrunk <- gsub("\\..*","",gen)
  
  chrEnsembl=chrUCSC2ENSEMBL(chr)
  
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
                                        chromosome = chrEnsembl, start = start, 
                                        end = end,  name = "genes ENSEMBL of interest",
                                        groupAnnotation = "group",
                                        just.group = "above",
                                        fontcolor="black",showId=showId,size=2,
                                        col.line = NULL, col = NULL)
    
    
  } else {
    martENSEMBL=useMart("ensembl",dataset='hsapiens_gene_ensembl')
    biomTrack <- BiomartGeneRegionTrack(genome = genTrunk, biomart=martENSEMBL,
                                        chromosome = chrEnsembl, start = start, 
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
transcript_ENSEMBL<-function(gen,chr,start,end,showId=FALSE){
  if(is.null(gen)){
    stop("Invalid in function transcriptENSEMBL :gen null:\n")
  }
  if(is.null(chr)){
    stop("Invalid in function transcriptENSEMBL :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid in function transcriptENSEMBL :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid in function transcriptENSEMBL :end null:\n")
  }
  
  genTrunk <- gsub("\\..*","",gen)
  
  chrEnsembl=chrUCSC2ENSEMBL(chr)
  
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
                                        chromosome = chrEnsembl, start = start, 
                                        end = end,  name = "transcripts ENSEMBL",
                                        fontcolor="black",groupAnnotation = "group",
                                        just.group = "above",showId=showId,size=2,
                                        col.line = NULL, col = NULL)
    
  } else {
    martENSEMBL=useMart("ensembl",dataset='hsapiens_gene_ensembl')
    biomTrack <- BiomartGeneRegionTrack(genome = genTrunk, biomart=martENSEMBL,
                                        chromosome = chrEnsembl, start = start, 
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
interestTranscript_ENSEMBL<-function(gen,chr,start,end,interestfeatures,interestcolor,showId=FALSE){
  if(is.null(gen)){
    stop("Invalid in function interestTranscriptENSEMBL: gen null:\n")
  }
  if(is.null(chr)){
    stop("Invalid in function interestTranscriptENSEMBL :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid in function interestTranscriptENSEMBL: start null:\n")
  }
  if(is.null(end)){
    stop("Invalid in function interestTranscriptENSEMBL: end null:\n")
  }
  if(is.null(interestfeatures)){
    stop("Invalid in function interestTranscriptENSEMBL: interestfeatures null:\n")
  }
  if(is.null(interestcolor)){
    stop("Invalid in function interestTranscriptENSEMBL: interestcolor null:\n")
  }
  
  genTrunk <- gsub("\\..*","",gen)
  
  chrEnsembl=chrUCSC2ENSEMBL(chr)
  
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
                                        chromosome = chrEnsembl, start = start, 
                                        end = end,  name = "transcript ENSEMBL of interest",
                                        fontcolor="black",groupAnnotation = "group",
                                        just.group = "above",showId=showId,size=2,
                                        col.line = NULL, col = NULL, collapse= FALSE)
    
  } else {
    martENSEMBL=useMart("ensembl",dataset='hsapiens_gene_ensembl')
    biomTrack <- BiomartGeneRegionTrack(genome = genTrunk, biomart=martENSEMBL,
                                        chromosome = chrEnsembl, start = start, 
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
chromatinHMMAll_UCSC<-function(gen,chr,start,end,mySession,color= 'coMET', pattern=NULL,table.name=NULL){
  if(is.null(gen)){
    stop("Invalid in function chromatinHMMAll :gen null:\n")
  }
  if(is.null(chr)){
    stop("Invalid in function chromatinHMMAll :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid in function chromatinHMMAll :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid in function chromatinHMMAll :end null:\n")
  }
  
  
  genTrunk <- gsub("\\..*","",gen)
  
  if(!is.na(match(tolower(gen), c("grch38","hg38")))){
    stop("Invalid in function chromatinHMMOne :genome not supported, only available for hg19\n")
  } else if( !is.na((match(tolower(genTrunk), c("hg19","grch37"))))){
    track.name="Broad ChromHMM"
  }else if( genTrunk != "hg19"){
    stop("Invalid in function chromatinHMMOne :genome not supported, only available for hg19\n")
  }
  track.name="Broad ChromHMM"
  tablestrack<-tableNames(ucscTableQuery(mySession, track=track.name))
  if(is.null(pattern)) {
    patterntable<-1:length(tablestrack)
  } else{
    patterntable<-grep(pattern, tablestrack,ignore.case=TRUE)
  }
  
  lltrack=list()
  for(i in patterntable){
    table.name<-tablestrack[i]
    tmp<-chromatinHMMOne_UCSC(genTrunk, chr, start, end, mySession, color, table.name)
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
chromatinHMMOne_UCSC<-function(gen,chr,start,end,mySession, color= 'coMET', table.name=NULL){
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
  if(!is.na(match(tolower(gen), c("grch38","hg38")))){
    stop("Invalid in function chromatinHMMOne :genome not supported, only available for hg19:\n")
  } else if( gen == "hg19"){
    track.name="Broad ChromHMM"
  }else if(gen != "hg19"){
    stop("Invalid in function chromatinHMMOne :genome not supported, only available for hg19:\n")
  }
  if(is.null(table.name)){
    table.name="wgEncodeBroadHmmHsmmHMM"
  }else if(is.null(table.name)){
    stop("Invalid in function chromatinHMMOne :gen null:\n")
  }
  
  track.name="Broad ChromHMM"
  colorcase <- tolower(color)
  mygrange <- GRanges(chr, IRanges(start, end))
  dataUCSC <- getTable(ucscTableQuery (mySession, range=mygrange, 
                                       track=track.name, table=table.name))
  
  data_trackfunc <- AnnotationTrack()
  if(nrow(dataUCSC) > 0) {
    data_trackfunc <- AnnotationTrack(genome=gen,chromosome=dataUCSC[,"chrom"],strand ="*",
                                      start=dataUCSC[,"chromStart"],end=dataUCSC[,"chromEnd"],
                                      feature=dataUCSC[,"name"],group=dataUCSC[,"name"],
                                      id=dataUCSC[,"name"], name = "Broad chromatinHMM",
                                      stacking="dense", col.line = NULL, col = NULL)
    chromosome(data_trackfunc) <- chr
    
    if(colorcase == "ucsc"){
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
      
    }else if (colorcase == "comet"){
      displayPars(data_trackfunc) <- list(
        "1_Active_Promoter"= "#E31A1C",
        "2_Weak_Promoter"="#FB9A99"  ,
        "3_Poised_Promoter"="#6A3D9A",
        "4_Strong_Enhancer"= "#FF7F00",
        "5_Strong_Enhancer"= "#CAB2D6",
        "6_Weak_Enhancer"="#FFFF99",
        "7_Weak_Enhancer"="#FDBF6F",
        "8_Insulator"="#1F78B4",
        "9_Txn_Transition"="#B2DF8A",
        "10_Txn_Elongation"="#33A02C",
        "11_Weak_Txn"="#00E1EF",
        "12_Repressed"="#FF00FF",
        "13_Heterochrom/lo"="#806000",
        "14_Repetitive/CNV"="#808080",
        "15_Repetitive/CNV"="#BFBFBF")
    }else {
      stop("Invalid in function chromatinHMMOne :color choice invalid :\n")
    }
    
  } else {
    data_trackfunc <- AnnotationTrack()
    chromosome(data_trackfunc) <- chr
    start(data_trackfunc) <- start
    end(data_trackfunc) <- end
  }
  
  data_trackfunc
}

#-------------------- CREATION track all types of Histone density from UCSC ------------------
HistoneAll_UCSC<-function(gen,chr,start,end,mySession,pattern=NULL,track.name="Broad Histone",table.name=NULL){
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
  if(is.null(track.name) & (!is.na(match(tolower(gen), c("hg19","grch37"))))){
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
    tmp<-HistoneOne_UCSC(gen,chr,start,end,mySession,track.name,table.name)
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
HistoneOne_UCSC<-function(gen,chr,start,end,mySession,track.name="Broad Histone",table.name=NULL){
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
  if(is.null(track.name) & ( !is.na(match(tolower(gen), c("hg19","grch37"))))){
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
    data_trackfunc <- AnnotationTrack(genome=gen,chromosome=dataUCSC[,"chrom"],strand ="*",
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
DNAse_UCSC<-function(gen,chr,start,end,mySession,track.name="DNase Clusters",table.name=NULL){
  if(is.null(gen)){
    stop("Invalid in function DNAseUCS :gen null:\n")
  }
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
  if(is.null(track.name) & (!is.na(match(tolower(gen), c("hg19","grch37"))))){
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
    data_trackfunc <- AnnotationTrack(genome=gen, chromosome=dataUCSC[,"chrom"],strand ="*",
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
gcContent_UCSC <- function(gen,chr,start,end){
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
            from = start,    to = end, trackType = "DataTrack", start = "start", 
            end = "end", data = "score", type = "hist", window = -1,    windowSize = 1500, 
            fill.histogram = "black",    col.histogram = "red", ylim = c(30, 70), 
            name = "GC Percent", col.line = NULL, col = NULL)
}

#-------------------- CREATION track Known genes from UCSC ------------------
knownGenes_UCSC<-function(gen,chr,start,end,showId=TRUE){
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
    knowngenestrack <- UcscTrack(genome = genTrunk, chromosome = chr,track = "knownGene", from = start, to = end, 
                                 trackType = "GeneRegionTrack", rstarts = "exonStarts", rends = "exonEnds", 
                                 gene = "name", symbol = "name", transcript = "name", strand = "strand", 
                                 fill = "#8282d2", name = "UCSC known Genes",stacking="squish", group="name",
                                 fontcolor="black", groupAnnotation = "group", just.group = "above",
                                 size=2, showId=TRUE,col.line = NULL, col = NULL)
  } else {
    knowngenestrack <- UcscTrack(genome = genTrunk, chromosome = chr,track = "knownGene", from = start, to = end, 
                                 trackType = "GeneRegionTrack", rstarts = "exonStarts", rends = "exonEnds", 
                                 gene = "name", symbol = "name", transcript = "name", strand = "strand", 
                                 fill = "#8282d2", name = "UCSC Genes",stacking="dense",size=2,
                                 col.line = NULL, col = NULL)
  }
  
  knowngenestrack
}

#-------------------- CREATION track reference genes from UCSC ------------------
refGenes_UCSC <-function(gen,chr,start,end, IdType = "Ref", showId=TRUE){
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
  
  IdTypecase <- tolower(IdType)
  
  if (IdTypecase == "name"){
    IDShow <- "name2"
  }
  else if (IdTypecase == "ref"){
    IDShow <- "name"
  }
  else{
    stop("Invalid in function refGenesUCSC :Invalid IdType:\n")
  }
  genTrunk <- gsub("\\..*","",gen)
  
  if(showId){
    refgenestrack <- UcscTrack(genome = genTrunk, chromosome = chr,track = "refGene", from = start, to = end,
                               trackType = "GeneRegionTrack", rstarts = "exonStarts", rends = "exonEnds",
                               gene = IDShow, symbol = IDShow, transcript = IDShow, strand = "strand",
                               fill = "#8282d2", name = "Ref Genes UCSC",stacking="squish", group=IDShow,
                               fontcolor="black", groupAnnotation = "group", just.group = "above",
                               size=2, showId=TRUE,col.line = NULL, col = NULL)
  } else {
    refgenestrack <- UcscTrack(genome = genTrunk, chromosome = chr,track = "knownGene", from = start, to = end,
                               trackType = "GeneRegionTrack", rstarts = "exonStarts", rends = "exonEnds",
                               gene = IDShow, symbol = IDShow, transcript = IDShow, strand = "strand",
                               fill = "#8282d2", name = "ref Genes",stacking="dense",size=2,
                               col.line = NULL, col = NULL)
  }
  
  refgenestrack 
}

#-------------------- CREATION track ref Genes from UCSC ------------------
xenorefGenes_UCSC<-function(gen,chr,start,end,showId=FALSE){
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
cpgIslands_UCSC <-function(gen,chr,start,end){
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
snpLocations_UCSC <-function(gen,chr,start,end,track="All SNPs(142)"){ 
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
regulationBiomart_ENSEMBL <- function(gen, chr, start, end) {
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
    data_trackfunc <- AnnotationTrack(genome=genTrunk, chromosome=chrEnsembl,strand="*",start=ensfunc[,2],
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
snpBiomart_ENSEMBL <- function(gen,chr, start, end, dataset, showId=FALSE, title_track=NULL) {
  if(is.null(gen)){
    stop("Invalid in function snpBiomart :gen null:\n")
  }
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
  if(is.null(title_track)){
    title_track="Short Variation"
  }
  
  genTrunk <- gsub("\\..*","",gen)
  chrEnsembl=chrUCSC2ENSEMBL(chr)
  martsnp=NULL
  
  if(!is.na(match(tolower(genTrunk), c("grch37","hg19")))){
    martsnp=useMart(host='grch37.ensembl.org', biomart='ENSEMBL_MART_SNP',
                    dataset=dataset)
  } else {
    
    martsnp <- useMart('snp',dataset=dataset)
  }
  
  ens_snp <- getBM(c("refsnp_id","chrom_start","chrom_strand","chr_name"),
                   filters = c("chr_name","start","end"),
                   values = list(chrEnsembl, start, end), mart=martsnp) 
  
  data_tracksnp <- AnnotationTrack()
  if(nrow(ens_snp) > 0) {
    data_tracksnp <- AnnotationTrack(genome=genTrunk,chromosome=ens_snp[,4],strand ="*",start=ens_snp[,2],
                                     end=ens_snp[,2],feature="snp",group=ens_snp[,1],
                                     id=ens_snp[,1], name = title_track,stacking="dense",
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
structureBiomart_ENSEMBL <- function(gen,chr, start, end, strand, dataset,showId=FALSE,title_track=NULL) {
  if(is.null(gen)){
    stop("Invalid in function snpBiomart :gen null:\n")
  }
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
  if(is.null(title_track)){
    title_track="Structural Variation"
  }
  genTrunk <- gsub("\\..*","",gen)
  chrEnsembl=chrUCSC2ENSEMBL(chr)
  martstruct=NULL
  
  if(!is.na(match(tolower(genTrunk), c("grch37","hg19")))){
    martstruct=useMart(host='grch37.ensembl.org', biomart='ENSEMBL_MART_SNP',
                       dataset=dataset)
  } else {
    
    martstruct <- useMart('snp',dataset=dataset)
  }
  
  ens <- getBM(c("sv_accession","chrom_start","chrom_end","seq_region_strand","chr_name",
                 "sv_variant_type","dgva_study_accession"),
               filters = c("chr_name","start","end"),
               values = list(chrEnsembl, start, end), mart=martstruct) 
  
  data_track <- AnnotationTrack()
  if(nrow(ens) > 0) {
    data_track <- AnnotationTrack(genome=genTrunk, chromosome=chr,strand =ens[,4],start=ens[,2],end=ens[,3],
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
ClinVarMain_UCSC <-function(gen,chr,start,end,showId=FALSE){ 
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
ClinVarCnv_UCSC <-function(gen,chr,start,end,showId=FALSE){ 
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
CoreillCNV_UCSC <-function(gen,chr,start,end,showId=FALSE){ 
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
COSMIC_UCSC <-function(gen,chr,start,end,showId=FALSE){ 
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
GAD_UCSC <-function(gen,chr,start,end,showId=FALSE){ 
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
GWAScatalog_UCSC <-function(gen,chr,start,end,showId=FALSE){ 
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
GeneReviews_UCSC <-function(gen,chr,start,end,showId=FALSE){ 
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
ISCA_UCSC <-function(gen,chr,start,end,mySession,table.name,showId=FALSE){ 
  if(is.null(gen)){
    stop("Invalid in function ISCA:gen null:\n")
  }
  if(is.null(chr)){
    stop("Invalid in function ISCA:chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid in function ISCA:start null:\n")
  }
  if(is.null(end)){
    stop("Invalid in function ISCA:end null:\n")
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
    data_trackfunc <- AnnotationTrack(genome=genTrunk, chromosome=dataUCSC[,"chrom"],strand ="*",
                                      start=dataUCSC[,"chromStart"],end=dataUCSC[,"chromEnd"],
                                      feature=dataUCSC[,"score"],group=dataUCSC[,"name"],
                                      id=dataUCSC[,"name"], name = "ISCA",
                                      stacking = "squish", showId=TRUE,
                                      groupAnnotation = "group",just.group = "above",
                                      col.line = NULL, col = NULL)
    chromosome(data_trackfunc) <- chr
    
    if(table.name == "iscaPathogenic") {
      displayPars(data_trackfunc) <- list(fill="#E41A1C")
    }
    if(table.name == "iscaPathGainCum") {
      displayPars(data_trackfunc) <- list(fill="#377EB8")
    }
    if(table.name == "iscaPathLossCum") {
      displayPars(data_trackfunc) <- list(fill="#4DAF4A")
    }
    if(table.name == "iscaCuratedPathogenic") {
      displayPars(data_trackfunc) <- list(fill="#984EA3")
    }
    if(table.name == "iscaLikelyPathogenic") {
      displayPars(data_trackfunc) <- list(fill="#FF7F00")
    }
    if(table.name == "iscaUncertain") {
      displayPars(data_trackfunc) <- list(fill="#FFFF33")
    }
    if(table.name == "iscaBenign") {
      displayPars(data_trackfunc) <- list(fill="#A65628")
    }
    if(table.name == "iscaCuratedBenign") {
      displayPars(data_trackfunc) <- list(fill="#F781BF")
    }
    if(table.name == "iscaLikelyBenign") {
      displayPars(data_trackfunc) <- list(fill="#999999")
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
repeatMasker_UCSC <-function(gen,chr,start,end,showId=FALSE,type_stacking="full"){ 
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
            shape = "box", stacking=type_stacking, fill = "grey", name = "RepeatMasker",
            groupAnnotation = "group",just.group = "above", 
            showId=showId, col.line = NULL, col = NULL)
}

## New database in Regulation ENSEMBL
#-------------------- CREATION track for all other regulatory regions from ENSEMBL or a list of them ------------
otherRegulatoryRegions_ENSEMBL <- function (gen, chr, start, end, featureDisplay = "all",datasetEnsembl = "hsapiens_external_feature") {
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
  if(is.null(featureDisplay)){
    stop("Invalid function OtherRegulatoryRegions :featureDisplay null:\n")
  }
  if(is.null(datasetEnsembl)){
    stop("Invalid function OtherRegulatoryRegions :datasetEnsembl null:\n")
  }
  
  options(ucscChromosomeNames=FALSE)
  chrEnsembl <- chrUCSC2ENSEMBL(chr)
  
  biomTrack <- NULL
  martfunc <- NULL
  biomTrackDisplay <- NULL
  
  genTrunk <- gsub("\\..*","",gen)
  
  if(!is.na(match(tolower(genTrunk), c("grch37","hg19")))){
    martfunc <- useMart(host='grch37.ensembl.org', biomart='ENSEMBL_MART_FUNCGEN',
                        dataset='hsapiens_external_feature')
    
  } else if(!is.na(match(tolower(genTrunk), c("grch38","hg38")))) {
    martfunc <- useMart(biomart="regulation", dataset="hsapiens_external_feature")
    
  } else if(!is.null(datasetEnsembl) && !is.na(datasetEnsembl)) {
    martfunc <- useMart(biomart="regulation", dataset=datasetEnsembl)
    
  } else {
    stop("Invalid function OtherRegulatoryRegions :genome not recognised:\n")
  }
  
  biomTrack <- getBM(c("chromosome_name", "chromosome_start","chromosome_end",
                       "feature_type_class"), 
                     filters = c("chromosome_name","start","end"),
                     values = list(chrEnsembl, start, end), mart = martfunc)
  
  biomTrackDisplay <- biomTrack
  
  if( !("all" %in% featureDisplay) ) {
    biomTrackDisplay <- biomTrack[which(biomTrack$feature_type_class %in% featureDisplay),]
  }
  if (nrow(biomTrackDisplay) == 0){
    biomTrackDisplay <- data.frame(nrow=1,ncol=4)
    biomTrackDisplay[1,1] <- chr
    biomTrackDisplay[1,2] <- start
    biomTrackDisplay[1,3] <- end
    biomTrackDisplay[1,4] <- 'Empty'
  }
  
  data_trackfunc <- AnnotationTrack(genome = genTrunk,chromosome=chrEnsembl,strand ="*",start=biomTrackDisplay[,2],
                                    end=biomTrackDisplay[,3],
                                    feature=biomTrackDisplay[,4],group=biomTrackDisplay[,1],
                                    id=biomTrackDisplay[,1], name = "Other Regulatory Regions ENSEMBL",stacking="dense", 
                                    col.line = "black", col = NULL, collapse= FALSE)
  
  displayPars(data_trackfunc) <- list("Enhancer" = "#e41a1c", "Transcription Start Site" = "#4daf4a","Empty" = "#ffffff")
  
  #print(biomTrackDisplay)
  
  data_trackfunc
  
} 

#-------------------- CREATION track for all Regulatory Evidence from ENSEMBL or a list of them ----------------
regulatoryEvidenceBiomart_ENSEMBL <- function (gen, chr, start, end, featureDisplay = "all",datasetEnsembl = "hsapiens_annotated_feature") {
  if(is.null(gen)){
    stop("Invalid function RegulatoryEvidenceBiomart :gen null:\n")
  }
  if(is.null(chr)){
    stop("Invalid function RegulatoryEvidenceBiomart :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid function RegulatoryEvidenceBiomart :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid function RegulatoryEvidenceBiomart :end null:\n")
  }
  if(is.null(featureDisplay)){
    stop("Invalid function RegulatoryEvidenceBiomart :featureDisplay null:\n")
  }
  if(is.null(datasetEnsembl)){
    stop("Invalid function RegulatoryEvidenceBiomart :datasetEnsembl null:\n")
  }
  
  options(ucscChromosomeNames=FALSE)
  chrEnsembl <- chrUCSC2ENSEMBL(chr)
  
  biomTrack <- NULL
  martfunc <- NULL
  biomTrackDisplay <- NULL
  
  genTrunk <- gsub("\\..*","",gen)
  
  if(!is.na(match(tolower(genTrunk), c("grch37","hg19")))){
    martfunc <- useMart(host='grch37.ensembl.org', biomart='ENSEMBL_MART_FUNCGEN',
                        dataset='hsapiens_annotated_feature')
    
  } else if(!is.na(match(tolower(genTrunk), c("grch38","hg38")))) {
    martfunc <- useMart(biomart="regulation", dataset="hsapiens_annotated_feature")
    
  } else if(!is.null(datasetEnsembl) && !is.na(datasetEnsembl)) {
    martfunc <- useMart(biomart="regulation", dataset=datasetEnsembl)
    
  } else {
    stop("Invalid function RegulatoryEvidenceBiomart :genome not recognised:\n")
  }
  
  biomTrack <- getBM(c("chromosome_name", "chromosome_start","chromosome_end",
                       "feature_type_name"),
                     filters = c("chromosome_name","start","end"),
                     values = list(chrEnsembl, start, end), mart = martfunc)
  
  biomTrackDisplay <- biomTrack
  if( !("all" %in% featureDisplay) ) {
    biomTrackDisplay <- biomTrack[which(biomTrack$feature_type_name %in% featureDisplay),]
  }
  if (nrow(biomTrackDisplay) == 0){
    biomTrackDisplay <- data.frame(nrow=1,ncol=4)
    biomTrackDisplay[1,1] <- chr
    biomTrackDisplay[1,2] <- start
    biomTrackDisplay[1,3] <- end
    biomTrackDisplay[1,4] <- 'Empty'
  }
  
  data_trackfunc <- AnnotationTrack(genome = genTrunk,chromosome=chrEnsembl,strand ="*",start=biomTrackDisplay[,2],
                                    end=biomTrackDisplay[,3],
                                    feature=biomTrackDisplay[,4],group=biomTrackDisplay[,1],
                                    id=biomTrackDisplay[,1], name = "Other Regulatory Regions ENSEMBL",stacking="dense",
                                    col.line = "black", col = NULL, collapse= FALSE)
  
  displayPars(data_trackfunc) <- list('H3K23me2' = '#41B4EE', 'H2BK5ac' = '#CDB47B','H3K9me1' = '#D5DEF6', 'MEF2C' = '#184173',
                                      'PolIII' = '#737B7B', 'XRCC4' = '#5A4120', 'BHLHE40' = '#60A8BC', 'H3K23ac' = '#E3E186',
                                      'NR4A1' = '#EE6241', 'HDAC8' = '#BD3162', 'ZNF274' = '#00FF00', 'Junb' = '#17B103',
                                      'BAF170' = '#E5DE00', 'H3K79me1' = '#525262', 'H2BK15ac' = '#8B8BA4', 'H3K14ac' = '#98FFA2',
                                      'ZEB1' = '#624A31', 'BAF155' = '#39394A', 'GTF2B' = '#FFE699', 'Gata1' = '#FFCDBD', 'THAP1'= '#B35F41',
                                      'SP2' = '#C8C8FA', 'Nrf1' = '#E66294', 'HNF4G' = '#94ADFF', 'Nfe2' = '#4C31AF', 'POU2F2' = '#FF40FF',
                                      'H2AK5ac' = '#767582', 'Pbx3' = '#000000', 'ETS1' = '#415A20', 'NFKB' = '#BDE673', 'H2BK12ac' = '#5DCF8B',
                                      'Nanog' = '#9C8B31', 'BCL3' = '#319C73', 'ZBTB7A' = '#804000', 'H3K56ac' = '#BD4A73', 'RXRA' = '#628BBD',
                                      'POU5F1' = '#A83D4C', 'CTCFL' = '#6FE9FF', 'FOXA2' = '#FFC22C', 'BCLAF1' = '#BF7EFF', 'SETDB1' = '#83C944',
                                      'H2BK20ac' = '#EEB4CD', 'FOSL1' = '#8B1608', 'Brg1' = '#C80096', 'Znf263' = '#DDFF00', 'Pax5' = '#C54129',
                                      'ZBTB33' = '#006A62', 'IRF4' = '#41ACA4', 'ATF3' = '#AC5273', 'H4K91ac' = '#83ACEE', 'Ini1' = '#A4CDE6',
                                      'FOSL2' = '#DEC552', 'SIX5' = '#D0DAD3', 'H4K8ac' = '#7B6220', 'Tr4' = '#4A6A83', 'BCL11A' = '#817E7A',
                                      'Gata2' = '#A87E7F', 'HDAC2' = '#6A6A52', 'MEF2A' = '#FF8158', 'EBF1' = '#00FDFF', 'p300' = '#393997',
                                      'E2F4' = '#AC3929', 'Ap2alpha' = '#D5A4DE', 'Sin3Ak20' = '#A55593', 'Srf' = '#8B52AC', 'Tcf12' = '#62317B',
                                      'H3K4ac' = '#C59700', 'HEY1' = '#C55A6A', 'H2BK120ac' = '#00FCC4', 'H4K5ac' = '#56CA4B', 'TAF7' = '#4C4C4C',
                                      'H3K18ac' = '#BDB46A', 'Cfos' = '#623918', 'E2F1' = '#EEFF94', 'HNF4A' = '#526229', 'BATF' = '#EE9C39',
                                      'SP1' = '#5900FF', 'Ap2gamma' = '#AC8B41', 'H3K9me3' = '#552431', 'E2F6' = '#186A88', 'ELF1' = '#9C5A4A',
                                      'PU1' = '#67B339', 'FOXA1' = '#F6BD5A', 'Jund' = '#9F847F', 'Gabp' = '#104131', 'Egr1' = '#83DEA4', 'Nrsf' = '#8BD5EE',
                                      'USF1' = '#986664', 'Cjun' = '#EEAC9C', 'H4K20me1' = '#085A73', 'Cmyc' = '#FF93F0', 'Yy1' = '#7B2941', 'TAF1' = '#C5C5C5',
                                      'Max' = '#20206A', 'H3K79me2' = '#D54152', 'H2AZ' = '#CCEBC5', 'H3K4me1' = '#83B420', 'Rad21' = '#7B9CB4',
                                      'PolII' = '#0070C0', 'H3K9ac' = '#52834A', 'H3K27ac' = '#08ACD5', 'H3K4me2' = '#FF603D', 'H3K4me3' = '#9D49C7',
                                      'DNase1' = '#FFFF00', 'CTCF' = '#0432FF', 'H3K27me3' = '#8EFF0E', 'H3K36me3' = '#FF0000',"Empty" = "#ffffff")
  data_trackfunc
  
}


#-------------------- CREATION track for all Regulatory features from ENSEMBL or a list of them ----------------                     
regulatoryFeaturesBiomart_ENSEMBL  <- function (gen, chr, start, end, featureDisplay = "all",datasetEnsembl = "hsapiens_regulatory_feature") {
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
  if(is.null(featureDisplay)){
    stop("Invalid function RegulatoryFeaturesBiomart :featureDisplay null:\n")
  }
  if(is.null(datasetEnsembl)){
    stop("Invalid function RegulatoryFeaturesBiomart :datasetEnsembl null:\n")
  }
  
  options(ucscChromosomeNames=FALSE)
  chrEnsembl <- chrUCSC2ENSEMBL(chr)
  
  biomTrack <- NULL
  martfunc <- NULL
  biomTrackDisplay <- NULL
  
  genTrunk <- gsub("\\..*","",gen)
  
  if(!is.na(match(tolower(genTrunk), c("grch37","hg19")))){
    martfunc <- useMart(host='grch37.ensembl.org', biomart='ENSEMBL_MART_FUNCGEN',
                        dataset='hsapiens_regulatory_feature')
    
  } else if(!is.na(match(tolower(genTrunk), c("grch38","hg38")))) {
    martfunc <- useMart(biomart="regulation", dataset="hsapiens_regulatory_feature")
    
  } else if(!is.null(datasetEnsembl) && !is.na(datasetEnsembl)) {
    martfunc <- useMart(biomart="regulation", dataset=datasetEnsembl)
    
  } else {
    stop("Invalid function RegulatoryFeaturesBiomart :genome not recognised:\n")
  }
  
  biomTrack <- getBM(c("chromosome_name", "chromosome_start","chromosome_end",
                       "feature_type_name"), 
                     filters = c("chromosome_name","start","end"),
                     values = list(chrEnsembl, start, end), mart = martfunc)
  
  biomTrackDisplay <- biomTrack
  
  if( !("all" %in% featureDisplay) ) {
    biomTrackDisplay <- biomTrack[which(biomTrack$feature_type_name %in% featureDisplay),]
  }
  if (nrow(biomTrackDisplay) == 0){
    biomTrackDisplay <- data.frame(nrow=1,ncol=4)
    biomTrackDisplay[1,1] <- chr
    biomTrackDisplay[1,2] <- start
    biomTrackDisplay[1,3] <- end
    biomTrackDisplay[1,4] <- 'Empty'
  }
  
  data_trackfunc <- AnnotationTrack(genome = genTrunk,chromosome=chrEnsembl,strand="*",start=biomTrackDisplay[,2],
                                    end=biomTrackDisplay[,3],
                                    feature=biomTrackDisplay[,4],group=biomTrackDisplay[,1],
                                    id=biomTrackDisplay[,1], name = "Regulatory Features ENSEMBL",stacking="dense", 
                                    col.line = "black", col = NULL, collapse= FALSE)
  
  displayPars(data_trackfunc) <- list("Promoter" = "#1b9e77", "TF binding site" = "#d95f02",
                                      "Open chromatin" = "#7570b3", "Promoter Flanking Region" =  "#e7298a",
                                      "CTCF Binding Site" = "#66a61e", "Enhancer" = "#e6ab02","Empty" = "#ffffff")
  
  data_trackfunc
  
} 

#-------------------- CREATION track for all Regulatory segments from ENSEMBL or a list of them ----------------
regulatorySegmentsBiomart_ENSEMBL  <- function (gen, chr, start, end,featureDisplay = 'all', datasetEnsembl = "hsapiens_segmentation_feature") {
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
  if(is.null(featureDisplay)){
    stop("Invalid function RegulatorySegmentsBiomart :featureDisplay null:\n")
  }
  if(is.null(datasetEnsembl)){
    stop("Invalid function RegulatorySegmentsBiomart :datasetEnsembl null:\n")
  }
  
  options(ucscChromosomeNames=FALSE)
  chrEnsembl <- chrUCSC2ENSEMBL(chr)
  
  biomTrack <- NULL
  martfunc <- NULL
  biomTrackDisplay <- NULL
  
  genTrunk <- gsub("\\..*","",gen)
  
  if(!is.na(match(tolower(genTrunk), c("grch37","hg19")))){
    martfunc <- useMart(host='grch37.ensembl.org', biomart='ENSEMBL_MART_FUNCGEN',
                        dataset='hsapiens_segmentation_feature')
    
  } else if(!is.na(match(tolower(genTrunk), c("grch38","hg38")))) {
    martfunc <- useMart(biomart="regulation", dataset="hsapiens_segmentation_feature")
    
  } else if(!is.null(datasetEnsembl) && !is.na(datasetEnsembl)) {
    martfunc <- useMart(biomart="regulation", dataset=datasetEnsembl)
    
  } else {
    stop("Invalid function RegulatorySegementsBiomart :genome not recognised:\n")
  }
  
  
  
  biomTrack <- getBM(c("chromosome_name", "chromosome_start","chromosome_end",
                       "feature_type_name"), 
                     filters = c("chromosome_name","start","end"),
                     values = list(chrEnsembl, start, end), mart = martfunc)
  
  biomTrackDisplay <- biomTrack
  
  if( !("all" %in% featureDisplay) ) {
    biomTrackDisplay <- biomTrack[which(biomTrack$feature_type_name %in% featureDisplay),]
  }
  
  if (nrow(biomTrackDisplay) == 0){
    biomTrackDisplay <- data.frame(nrow=1,ncol=4)
    biomTrackDisplay[1,1] <- chr
    biomTrackDisplay[1,2] <- start
    biomTrackDisplay[1,3] <- end
    biomTrackDisplay[1,4] <- 'Empty'
  }
  
  data_trackfunc <- AnnotationTrack(genome = genTrunk,chromosome=chrEnsembl,strand="*",start=biomTrackDisplay[,2],
                                    end=biomTrackDisplay[,3],
                                    feature=biomTrackDisplay[,4],group=biomTrackDisplay[,1],
                                    id=biomTrackDisplay[,1], name = "Regulatory Segments ENSEMBL",stacking="dense", 
                                    col.line = "black", col = NULL, collapse= FALSE)
  
  displayPars(data_trackfunc) <- list('Predicted Promoter with TSS' = '#a6cee3','CTCF enriched' = '#1f78b4',
                                      'Predicted Poised' = '#b2df8a', 'Predicted Promoter Flank' = '#33a02c',
                                      'Predicted Enhancer' = '#fb9a99', 'Predicted Transcribed Region' = '#e31a1c',
                                      'Predicted low activity' = '#fdbf6f', 'Predicted Repressed' = '#ff7f00',
                                      'Predicted heterochomatin' = '#cab2d6',"Empty" = "#ffffff")
  
  data_trackfunc
  
} 

#-------------------- CREATION track for all binding motifs from ENSEMBL or a list of them ----------------
bindingMotifsBiomart_ENSEMBL <- function (gen, chr, start, end, featureDisplay = "all",datasetEnsembl = NULL) {
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
  biomTrackDisplay <- NULL
  
  genTrunk <- gsub("\\..*","",gen)
  
  if(!is.na(match(tolower(genTrunk), c("grch37","hg19")))){
    martfunc <- useMart(host='grch37.ensembl.org', biomart='ENSEMBL_MART_FUNCGEN',
                        dataset='hsapiens_motif_feature')
    
  } else if(!is.na(match(tolower(genTrunk), c("grch38","hg38")))) {
    martfunc <- useEnsembl(biomart="regulation", dataset="hsapiens_motif_feature")
    
  } else if(!is.null(datasetEnsembl) && !is.na(datasetEnsembl)) {
    martfunc <- useMart(biomart="regulation", dataset=datasetEnsembl)
    
  } else {
    stop("Invalid function BindingMotifsBiomart :genome not recognised:\n")
  }
  
  biomTrack <- getBM(c("chromosome_name", "chromosome_start","chromosome_end",
                       "feature_type_name"), 
                     filters = c("chromosome_name","start","end"),
                     values = list(chrEnsembl, start, end), mart = martfunc)
  
  biomTrackDisplay <- biomTrack
  
  if( !("all" %in% featureDisplay) ) {
    biomTrackDisplay <- biomTrack[which(biomTrack$feature_type_name %in% featureDisplay),]
  }
  
  if (nrow(biomTrackDisplay) == 0){
    biomTrackDisplay <- data.frame(nrow=1,ncol=4)
    biomTrackDisplay[1,1] <- chr
    biomTrackDisplay[1,2] <- start
    biomTrackDisplay[1,3] <- end
    biomTrackDisplay[1,4] <- 'Empty'
  }
  
  data_trackfunc <- AnnotationTrack(genome = genTrunk,chromosome=chrEnsembl,strand ="*",start=biomTrackDisplay[,2],
                                    end=biomTrackDisplay[,3],
                                    feature=biomTrackDisplay[,4],group=biomTrackDisplay[,1],
                                    id=biomTrackDisplay[,1], name = "Binding Motifs ENSEMBL",stacking="dense",
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
    "SRebp2" = "#B2B2F0","Empty" = "#ffffff")
  
  data_trackfunc
} 

#-------------------- CREATION track for all miRNA Target from ENSEMBL ----------------
miRNATargetRegionsBiomart_ENSEMBL  <- function (gen, chr, start, end, showId=FALSE, datasetEnsembl = "hsapiens_mirna_target_feature") {
  if(is.null(gen)){
    stop("Invalid function miRNATargetRegionsBiomart :gen null:\n")
  }
  if(is.null(chr)){
    stop("Invalid function miRNATargetRegionsBiomart :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid function miRNATargetRegionsBiomart :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid function miRNATargetRegionsBiomart :end null:\n")
  }
  if(is.null(datasetEnsembl)){
    stop("Invalid function miRNATargetRegionsBiomart :datasetEnsembl null:\n")
  }
  
  options(ucscChromosomeNames=FALSE)
  chrEnsembl <- chrUCSC2ENSEMBL(chr)
  
  biomTrack <- NULL
  martfunc <- NULL
  biomTrackDisplay <- NULL
  
  genTrunk <- gsub("\\..*","",gen)
  
  if(!is.na(match(tolower(genTrunk), c("grch37","hg19")))){
    martfunc <- useMart(host='grch37.ensembl.org', biomart='ENSEMBL_MART_FUNCGEN',
                        dataset='hsapiens_mirna_target_feature')
    
  } else if(!is.na(match(tolower(genTrunk), c("grch38","hg38")))) {
    martfunc <- useMart(biomart="regulation", dataset="hsapiens_mirna_target_feature")
    
  } else if(!is.null(datasetEnsembl) && !is.na(datasetEnsembl)) {
    martfunc <- useMart(biomart="regulation", dataset=datasetEnsembl)
    
  } else {
    stop("Invalid function miRNATargetRegionsBiomart :genome not recognised:\n")
  }
  
  
  
  biomTrack <- getBM(c("chromosome_name", "chromosome_start","chromosome_end",
                       "chromosome_strand","feature_type_class","xref_display_label"), 
                     filters = c("chromosome_name","start","end"),
                     values = list(chrEnsembl, start, end), mart = martfunc)
  
  data_trackfunc <- AnnotationTrack(genome = genTrunk,chromosome=chrEnsembl,strand =biomTrack[,4],start=biomTrack[,2],
                                    end=biomTrack[,3],
                                    feature=biomTrack[,5],group=biomTrack[,1],
                                    id=biomTrack[,6], name = "miRNA Target Regions ENSEMBL",stacking="dense", 
                                    col.line = "plum4", col = NULL, collapse= FALSE,showId=showId)
  displayPars(data_trackfunc) <- list(RNA="plum4")
  
  data_trackfunc
  
} 

#-------------------- CREATION track for all Segmental duplications from UCSC  ----------------
segmentalDups_UCSC <- function(gen, chr, start, end){
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
                          end="chromEnd", id = "name", name = "Segmental Dups UCSC")
  
  Dupregions
  
}

#------------------- DNaseI element of ROADMap visualisation data --------------
ChIPTF_ENCODE <- function(gen="hg19",chr,start, end, bedFilePath, featureDisplay='all', motifColorFile, type_stacking="dense",showId=FALSE,just_group="above") {
  if(is.null(gen)){
    stop("Invalid function ChIPTF_ENCODE :gen null:\n")
  }
  if(is.null(chr)){
    stop("Invalid function ChIPTF_ENCODE :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid function ChIPTF_ENCODE :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid function ChIPTF_ENCODE :end null:\n")
  }
  if(is.null(featureDisplay)){
    stop("Invalid function ChIPTF_ENCODE :end null:\n")
  }
  if(is.null(bedFilePath)){
    stop("Invalid function ChIPTF_ENCODE :bedFilePath null:\n")
  }
  if(is.null(motifColorFile)){
    stop("Invalid function ChIPTF_ENCODE :motifColorFile null:\n")
  }
  
  chrEnsembl <- chrUCSC2ENSEMBL(tolower(chr))
  
  bedtab <- read.table(bedFilePath,header = FALSE,sep=" ")
  RegionDisplay <- bedtab[which((bedtab[,3] >= start | bedtab[,4]  <= end )
                                & bedtab[,2]  == chr),]
  
  desiredRegionDisplay <- RegionDisplay
  
  if( !("all" %in% featureDisplay) ) {
    desiredRegionDisplay <- RegionDisplay[which(RegionDisplay[,1] %in% featureDisplay),]
  }
  
  if (nrow(desiredRegionDisplay) == 0){
    desiredRegionDisplay <- data.frame(nrow=1,ncol=5)
    desiredRegionDisplay[1,1] <- chr
    desiredRegionDisplay[1,2] <- start
    desiredRegionDisplay[1,3] <- end
    desiredRegionDisplay[1,4] <- 'Empty'
    desiredRegionDisplay[1,5] <- 'Empty'
  } 
  
  track <- AnnotationTrack(genome=gen,chromosome=chrEnsembl,strand =desiredRegionDisplay[,5],start=desiredRegionDisplay[,3],
                           end=desiredRegionDisplay[,4],
                           feature=desiredRegionDisplay[,1], name = "TF motifs ENCODE",
                           id=desiredRegionDisplay[,1],just.group = just_group,
                           stacking=type_stacking, group=desiredRegionDisplay[,1],
                           col.line = "black", col = NULL, collapse= FALSE,showId=showId)
  
  colortab <- read.table(motifColorFile,header = TRUE,sep="\t") 
  colortab[,2] <- paste("#",colortab[,2],sep="")
  colortab <- data.frame(lapply(colortab, as.character), stringsAsFactors=FALSE)
  emptytab <- c("Empty", "#ffffff")
  colortab <- rbind(colortab,emptytab)
  colorTrack <-lapply(seq_len(nrow(colortab)), function(i) colortab[i,2])
  names(colorTrack) <-  colortab[,1]
  displayPars(track) <-colorTrack
  
  track
}

#------------------- Hi-C data ------------------------------
#------------------ Create matrix from Rao data -------------
HiCdata2matrix <- function(chr,start, end, bedFilePath) {
  if(is.null(chr)){
    stop("Invalid function HiCdata2matrix :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid function HiCdata2matrix :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid function HiCdata2matrix :end null:\n")
  }
  if(is.null(bedFilePath)){
    stop("Invalid function HiCdata2matrix :bedFilePath null:\n")
  }
  
  chrEnsembl <- chrUCSC2ENSEMBL(tolower(chr))
  
  data_intrachr_HiC <- read.csv(bedFilePath, header=FALSE, sep = "\t", quote = "")
  data_intrachr_HiC <- as.matrix(data_intrachr_HiC)
  
  start_interest <- start
  end_interest <- end
  list_bins <- which(data_intrachr_HiC[,1] >= start_interest &
                       data_intrachr_HiC[,1] <= end_interest  &
                       data_intrachr_HiC[,2] >= start_interest &
                       data_intrachr_HiC[,2] <= end_interest)
  
  subdata_intrachr_HiC <- data_intrachr_HiC[list_bins, ]
  matrix_HiC <- NULL
  if(nrow(subdata_intrachr_HiC) >0) {
    regions <- NULL
    regions <- sort(unique(subdata_intrachr_HiC[,1]))
    matrix_HiC<-matrix(nrow=length(regions),ncol=length(regions))
    
    colnames(matrix_HiC) <- regions 
    rownames(matrix_HiC) <- regions 
    
    diag(matrix_HiC)<-1
    for(i in 1:nrow(subdata_intrachr_HiC)){
      num_feat1 <- which(rownames(matrix_HiC) == subdata_intrachr_HiC[i,1])
      num_feat2 <- which(rownames(matrix_HiC) == subdata_intrachr_HiC[i,2])
      matrix_HiC[num_feat1,num_feat2] <- subdata_intrachr_HiC[i,3]
      matrix_HiC[num_feat2,num_feat1] <- subdata_intrachr_HiC[i,3]
    }
  }
  
  matrix_HiC
}

#------------------- ROADMap visualisation data --------------
#------------------- chromHMM of ROADMap visualisation data --------------
chromHMM_RoadMap <- function(gen="hg19",chr,start, end, bedFilePath, featureDisplay = 'all', colorcase='roadmap15') {
  
  if(is.null(gen)){
    stop("Invalid function chromHMM_RoadMap :gen null:\n")
  }
  if(is.null(chr)){
    stop("Invalid function chromHMM_RoadMap :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid function chromHMM_RoadMap :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid function chromHMM_RoadMap :end null:\n")
  }
  if(is.null(bedFilePath)){
    stop("Invalid function chromHMM_RoadMap :bedFilePath null:\n")
  }
  
  chrEnsembl <- chrUCSC2ENSEMBL(tolower(chr))
  
  bedFile <- read.table(bedFilePath,header = FALSE,sep="\t")
  colnames(bedFile)<-c("chromosome_name","chromosome_start","chromosome_stop","feature_type_name")
  desiredRegion <- subset(bedFile, chromosome_stop > start & chromosome_start < end &  chromosome_name == chr)
  
  desiredRegionDisplay <- desiredRegion
  
  if( !("all" %in% featureDisplay) ) {
    desiredRegionDisplay <- desiredRegion[which(desiredRegion$feature_type_name %in% featureDisplay),]
  }
  
  if (nrow(desiredRegionDisplay) == 0){
    desiredRegionDisplay <- data.frame(nrow=1,ncol=4)
    desiredRegionDisplay[1,1] <- chr
    desiredRegionDisplay[1,2] <- start
    desiredRegionDisplay[1,3] <- end
    desiredRegionDisplay[1,4] <- 'Empty'
  }
  
  track <- AnnotationTrack(genome=gen,chromosome=chrEnsembl,strand ="*",start=desiredRegionDisplay[,2],
                           end=desiredRegionDisplay[,3],
                           feature=desiredRegionDisplay[,4], name = paste(colorcase," chromHMM RoadMap"),stacking="dense",
                           col.line = "black", col = NULL, collapse= FALSE)
  
  if(colorcase == "roadmap15"){
    displayPars(track) <- list("1_TssA" = "#FF0000", "2_TssAFlnk" = "#FF6E00","3_TxFlnk" = "#32CD32",
                               "4_Tx" = "#008000", "5_TxWk" = "#006400", "6_EnhG" = "#C2E105", "7_Enh" = "#FFFF00",
                               "8_ZNF/Rpts" = "#66CDAA", "9_Het" = "#8A91D0", "10_TssBiv" = "#CD5C5C",
                               "11_BivFlnk" = "#E9967A", "12_EnhBiv" = "#BDB76B", "13_ReprPC" = "#3A3838",
                               "14_ReprPCWk" = "#808080", "15_Quies" = "#DCDCDC","Empty" = "#ffffff")
    
  } else if(colorcase == "roadmap18"){
    displayPars(data_trackfunc) <- list("1_TssA" = "#FF0000", "2_TssFlnk" = "#FF4500", "3_TssFlnkU" = "#FF4500", 
                                        "4_TssFlnkD" = "#FF4500", "5_Tx" = "#008000", "6_TxWk" = "#006400",
                                        "7_EnhG1" = "#C2FF05", "8_EnhG2" = "#C2FF05", "9_EnhA1" = "#FFC34D",
                                        "10_EnhA2" = "#FFC34D", "11_EnhWk" = "#FFFF00", "12_ZNF/Rpts" = "#66CDAA",
                                        "13_Het" = "#8A91D0", "14_TssBiv" = "#CD5C5C", "15_EnhBiv" = "#BDB76B",
                                        "16_ReprPC" = "#808080", "17_ReprPC" = "#C0C0C0", "18_Quies" = "#FFFFFF")
    
  }else if (colorcase == "comet18"){
    displayPars(data_trackfunc) <- list("1_TssA" = "#FF0000", "2_TssFlnk" = "#FF6E00", "3_TssFlnkU" = "#FF9300", 
                                        "4_TssFlnkD" = "#DA7B08", "5_Tx" = "#008000", "6_TxWk" = "#006400",
                                        "7_EnhG1" = "#C2FF05", "8_EnhG2" = "#C2FFBD", "9_EnhA1" = "#FE00DB",
                                        "10_EnhA2" = "#FFA7D6", "11_EnhWk" = "#FFFF00", "12_ZNF/Rpts" = "#66CDAA",
                                        "13_Het" = "#8A91D0", "14_TssBiv" = "#CD5C5C", "15_EnhBiv" = "#BDB76B",
                                        "16_ReprPC" = "#323232", "17_ReprPC" = "#AFAFAF", "18_Quies" = "#DCDCDC")
  }else if(colorcase == "roadmap25"){
    displayPars(data_trackfunc) <- list("1_TssA" = "#FF0000", "2_PromU" = "#FF4500", "3_PromD1" = "#FF4500",
                                        "4_PromD2" = "#FF4500", "5_Tx5???" = "#008000", "6_Tx" = "#008000", 
                                        "7_Tx3???" = "#008000", "8_TxWk" = "#009600", "9_TxReg" = "#C2FF05",
                                        "10_TxEnh5???" = "#C2FF05", "11_TxEnh3???" = "#C2FF05", "12_TxEnhW" = "#C2FF05",
                                        "13_EnhA1" = "#FFC34D", "14_EnhA2" = "#FFC34D", "15_EnhAF" = "#FFC34D",
                                        "16_EnhW1" = "#FFFF00", "17_EnhW2" = "#FFFF00", "18_EnhAc" = "#FFFF00",
                                        "19_DNase" = "#FFFF66", "20_ZNF/Rpts" = "#66CDAA", "21_Het" = "#8A91D0",
                                        "22_PromP" = "#E6B8B7", "23_PromBiv" = "#7030A0", "24_ReprPC" = "#808080",
                                        "25_Quies" = "#FFFFFF")
    
    
  }else if (colorcase == "comet25"){
    displayPars(data_trackfunc) <- list("1_TssA" = "#FF0000", "2_PromU" = "#FC6D00", "3_PromD1" = "#DD8100",
                                        "4_PromD2" = "#AD7622", "5_Tx5???" = "#008000", "6_Tx" = "#004D00", 
                                        "7_Tx3???" = "#009462", "8_TxWk" = "#00FE00", "9_TxReg" = "#00FFFF",
                                        "10_TxEnh5???" = "#009FFF", "11_TxEnh3???" = "#0028FF", "12_TxEnhW" = "#0000AE",
                                        "13_EnhA1" = "#FF00FF", "14_EnhA2" = "#FFB2FF", "15_EnhAF" = "#FFD8FF",
                                        "16_EnhW1" = "#FFFF00", "17_EnhW2" = "#E3FF8C", "18_EnhAc" = "#FFD500",
                                        "19_DNase" = "#FFFFC2", "20_ZNF/Rpts" = "#66CDAA", "21_Het" = "#8A91D0",
                                        "22_PromP" = "#E6B8B7", "23_PromBiv" = "#7030A0", "24_ReprPC" = "#646464",
                                        "25_Quies" = "#DCDCDC")
  }else {
    #change function name (if neeeded)
    stop("Invalid in function RoadMap :color choice invalid :\n")
  }
  
  track
}

#------------------- DNaseI element of ROADMap visualisation data --------------
DNaseI_RoadMap <- function(gen="hg19",chr,start, end, bedFilePath, featureDisplay='promotor',showId=TRUE, type_stacking="dense") {
  if(is.null(gen)){
    stop("Invalid function DNaseI_RoadMap :gen null:\n")
  }
  if(is.null(chr)){
    stop("Invalid function DNaseI_RoadMap :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid function DNaseI_RoadMap :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid function DNaseI_RoadMap :end null:\n")
  }
  if(is.null(bedFilePath)){
    stop("Invalid function DNaseI_RoadMap :bedFilePath null:\n")
  }
  
  chrEnsembl <- chrUCSC2ENSEMBL(tolower(chr))
  
  bedtab <- read.table(bedFilePath,header = FALSE,sep="\t")
  desiredRegionDisplay <- bedtab[which(bedtab[,2] > start & bedtab[,3]  < end & bedtab[,1]  == chr), 1:4]
  
  if (nrow(desiredRegionDisplay) == 0){
    desiredRegionDisplay <- data.frame(nrow=1,ncol=5)
    desiredRegionDisplay[1,1] <- chr
    desiredRegionDisplay[1,2] <- start
    desiredRegionDisplay[1,3] <- end
    desiredRegionDisplay[1,4] <- ' '
    desiredRegionDisplay[1,5] <- 'Empty'
  } else {
    desiredRegionDisplay$feature_type <- featureDisplay
  }
  
  track <- AnnotationTrack(genome=gen,chromosome=chrEnsembl,strand ="*",start=desiredRegionDisplay[,2],
                           end=desiredRegionDisplay[,3],
                           feature=desiredRegionDisplay[,5], name = paste(featureDisplay,"RoadMap"),
                           id=desiredRegionDisplay[,4],
                           stacking=type_stacking,showId=showId,
                           col.line = "black", col = NULL, collapse= FALSE)
  
  displayPars(track) <- list("promoter" = "#FF0000", "enhancer" = "#006400",  
                             "dyadic" = "#8A91D0","Empty" = "#ffffff")
  
  track
}

#-------------------  DG footptints of ROADMap visualisation data --------------
dgfootprints_RoadMap <- function(gen="hg19",chr,start, end, bedFilePath, tissueGroupDisplay='Blood & T-cell',showId=FALSE, type_stacking="dense") {
  if(is.null(gen)){
    stop("Invalid function dgfootprints_RoadMap :gen null:\n")
  }
  if(is.null(chr)){
    stop("Invalid function dgfootprints_RoadMap :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid function dgfootprints_RoadMap :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid function dgfootprints_RoadMap :end null:\n")
  }
  if(is.null(bedFilePath)){
    stop("Invalid function dgfootprints_RoadMap :bedFilePath null:\n")
  }
  
  chrEnsembl <- chrUCSC2ENSEMBL(tolower(chr))
  
  bedtab <- read.table(bedFilePath,header = FALSE,sep="\t")
  desiredRegion <- bedtab[which(bedtab[,2] > start & bedtab[,3]  < end & bedtab[,1]  == chr), 1:4]
  
  desiredRegionDisplay <- desiredRegion
  
  if( !("all" %in% tissueGroupDisplay) ) {
    desiredRegionDisplay <- desiredRegion[which(desiredRegion$feature_type_name %in% tissueGroupDisplay),]
  }
  
  if (nrow(desiredRegionDisplay) == 0){
    desiredRegionDisplay <- data.frame(nrow=1,ncol=5)
    desiredRegionDisplay[1,1] <- chr
    desiredRegionDisplay[1,2] <- start
    desiredRegionDisplay[1,3] <- end
    desiredRegionDisplay[1,4] <- 'Empty'
    desiredRegionDisplay[1,5] <- 'Empty'
  } else {
    desiredRegionDisplay$feature_type <- tissueGroupDisplay
  }
  
  track <- AnnotationTrack(genome=gen,chromosome=chrEnsembl,strand ="*",start=desiredRegionDisplay[,2],
                           end=desiredRegionDisplay[,3],
                           feature=desiredRegionDisplay[,5], name = paste(tissueGroupDisplay,"DGFP RoadMap"),
                           id=desiredRegionDisplay[,4], just.group="above",
                           stacking=type_stacking, group=desiredRegionDisplay[,4],
                           col.line = "black", col = NULL, collapse= FALSE,showId=showId)
  
  displayPars(track) <- list("Neurosph"="#FFD924","Epithelial"="#FF9D0C",
                             "IMR90"="#E41A1C","Thymus"="#DAB92E",
                             "Heart"="#D56F80","Brain"="#C5912B","Digestive"="#C58DAA",
                             "Muscle"="#C2655D","Other"="#999999","iPSC"="#69608A",
                             "HSC & B-cell"="#678C69","Blood & T-cell"="#55A354",
                             "ES-deriv"="#4178AE","Empty" = "#ffffff")
  
  track
}

#-------------------  FANTOM5 --------------
#------------------- DNaseI element of ROADMap visualisation data --------------
DNaseI_FANTOM <- function(gen="hg19",chr,start, end, bedFilePath, featureDisplay='enhancer', stacking_type="dense") {
  if(is.null(gen)){
    stop("Invalid function DNaseI_FANTOM :gen null:\n")
  }
  if(is.null(chr)){
    stop("Invalid function DNaseI_FANTOM :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid function DNaseI_FANTOM :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid function DNaseI_FANTOM :end null:\n")
  }
  if(is.null(bedFilePath)){
    stop("Invalid function DNaseI_FANTOM :bedFilePath null:\n")
  }
  if(is.null(featureDisplay)){
    stop("Invalid function DNaseI_FANTOM :featureDisplay null:\n")
  }
  
  chrEnsembl <- chrUCSC2ENSEMBL(tolower(chr))
  
  bedtab <- read.table(bedFilePath,header = FALSE,sep="\t")
  desiredRegionDisplay <- bedtab[which(bedtab[,2] > start & bedtab[,3]  < end 
                                       & bedtab[,1]  == chr), 1:4]
  
  if (nrow(desiredRegionDisplay) == 0){
    desiredRegionDisplay <- data.frame(nrow=1,ncol=4)
    desiredRegionDisplay[1,1] <- chr
    desiredRegionDisplay[1,2] <- start
    desiredRegionDisplay[1,3] <- end
    desiredRegionDisplay[1,4] <- 'Empty'
  } else {
    desiredRegionDisplay$feature_type <- featureDisplay
  }
  
  track <- AnnotationTrack(genome=gen,chromosome=chrEnsembl,strand ="*",start=desiredRegionDisplay[,2],
                           end=desiredRegionDisplay[,3],
                           feature=desiredRegionDisplay[,5], name = paste(featureDisplay,"RoadMap"),
                           id=desiredRegionDisplay[,4], group= desiredRegionDisplay[,5],
                           stacking=stacking_type, just.group="above",
                           col.line = "black", col = NULL, collapse= FALSE)
  
  displayPars(track) <- list("promoter" = "#FF0000", "enhancer" = "#006400",  
                             "dyadic" = "#8A91D0","Empty" = "#ffffff")
  
  track
}

#------------------- TFBS motif visualisation data --------------
TFBS_FANTOM <- function(gen,chr,start, end, bedFilePath) {
  if(is.null(gen)){
    stop("Invalid function TFBS_FANTOM :gen null:\n")
  }
  if(is.null(chr)){
    stop("Invalid function TFBS_FANTOM :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid function TFBS_FANTOM :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid function TFBS_FANTOM :end null:\n")
  }
  if(is.null(bedFilePath)){
    stop("Invalid function TFBS_FANTOM :bedFilePath null:\n")
  }
  
  
  chrEnsembl <- chrUCSC2ENSEMBL(tolower(chr))
  
  bedtab <- read.table(bedFilePath,header = TRUE,sep="\t")
  desiredRegionDisplay <- bedtab[which(bedtab$start.0base > start & bedtab$end  < end &
                                         bedtab$chrom  == chr),]
  
  if (nrow(desiredRegionDisplay) == 0){
    desiredRegionDisplay <- data.frame(nrow=1,ncol=6)
    desiredRegionDisplay[1,1] <- chr
    desiredRegionDisplay[1,2] <- start
    desiredRegionDisplay[1,3] <- end
    desiredRegionDisplay[1,4] <- 'Empty'
    desiredRegionDisplay[1,5] <- 'Empty'
    desiredRegionDisplay[1,6] <- "*"
  } 
  
  track <- AnnotationTrack(genome=gen,chromosome=chrEnsembl,strand=desiredRegionDisplay[,6],
                           start=desiredRegionDisplay[,2],end=desiredRegionDisplay[,3],
                           feature=desiredRegionDisplay[,4], 
                           name = paste(desiredRegionDisplay[1,4],"TF motif FANTOM5"),
                           id=desiredRegionDisplay[,4],
                           stacking="dense",
                           col.line = "black", col = NULL, collapse= FALSE)
  
  track
}

#------------------- metQTL visualisation data --------------
metQTL <- function(gen,chr,start, end, bedFilePath, featureDisplay = 'all', showId=FALSE,type_stacking="squish",just_group="above" ) {
  if(is.null(gen)){
    stop("Invalid function metQTL :gen null:\n")
  }
  if(is.null(chr)){
    stop("Invalid function metQTL :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid function metQTL :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid function metQTL :end null:\n")
  }
  
  if(is.null(bedFilePath)){
    stop("Invalid function metQTL :bedFilePath null:\n")
  }
  
  # chrEnsembl <- chrUCSC2ENSEMBL(tolower(chr))
  
  bedFile <- read.table(bedFilePath,header = TRUE,sep="\t")
  desiredRegion <- subset(bedFile, chromosome_stop >= start & chromosome_start <= end &  chromosome_name == chr)
  
  desiredRegionDisplay <- desiredRegion
  
  if( !("all" %in% featureDisplay) ) {
    desiredRegionDisplay <- desiredRegion[which(desiredRegion$feature_type %in% featureDisplay),]
  }
  
  if (nrow(desiredRegionDisplay) == 0){
    desiredRegionDisplay <- data.frame(nrow=1,ncol=4)
    desiredRegionDisplay[1,1] <- chr
    desiredRegionDisplay[1,2] <- start
    desiredRegionDisplay[1,3] <- end
    desiredRegionDisplay[1,4] <- 'Empty'
  }
  
  desiredRegionDisplay <- data.frame(lapply(desiredRegionDisplay, as.character), stringsAsFactors=FALSE)
  track <- AnnotationTrack(genome=gen,chromosome=chr,strand =desiredRegionDisplay[,4],start=desiredRegionDisplay[,2],
                           end=desiredRegionDisplay[,3],
                           feature=desiredRegionDisplay[,5], group=desiredRegionDisplay[,7],
                           id=desiredRegionDisplay[,7], name = "metQTL", groupAnnotation = "group",
                           just.group = just_group, stacking=type_stacking,
                           col.line = "black", col = NULL, collapse= FALSE,showId=showId)
  
  displayPars(track) <- list("SNP_pheno" = "#A6CEE3", "SNP" = "#1F78B4", "CpG_pheno" = "#B2DF8A",
                             "CpG" ="#33A02C", "cis_local_metQTL" = "#FB9A99",
                             "trans_local_metQTL" = "#E31A1C", "distal_metQTL" = "#FDBF6F",
                             "cis_local_metQTL_pheno" = "#FF7F00",
                             "trans_local_metQTL_pheno" = "#CAB2D6",
                             "distal_metQTL_pheno" = "#6A3D9A","Empty" = "#ffffff")
  
  track
}

#------------------- etQTL visualisation data --------------
eQTL <- function(gen,chr,start, end, bedFilePath, featureDisplay = 'all', showId=FALSE, type_stacking="squish",just_group="above" ) {
  if(is.null(gen)){
    stop("Invalid function eQTL :gen null:\n")
  }
  if(is.null(chr)){
    stop("Invalid function eQTL :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid function eQTL :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid function eQTL :end null:\n")
  }
  if(is.null(bedFilePath)){
    stop("Invalid function eQTL :bedFilePath null:\n")
  }
  
  
  # chrEnsembl <- chrUCSC2ENSEMBL(tolower(chr))
  
  bedFile <- read.table(bedFilePath,header = TRUE,sep="\t")
  desiredRegion <- subset(bedFile, chromosome_stop >= start & chromosome_start <= end &  chromosome_name == chr)
  
  desiredRegionDisplay <- desiredRegion
  
  if( !("all" %in% featureDisplay) ) {
    desiredRegionDisplay <- desiredRegion[which(desiredRegion$feature_type %in% featureDisplay),]
  }
  if (nrow(desiredRegionDisplay) == 0){
    desiredRegionDisplay <- data.frame(nrow=1,ncol=4)
    desiredRegionDisplay[1,1] <- chr
    desiredRegionDisplay[1,2] <- start
    desiredRegionDisplay[1,3] <- end
    desiredRegionDisplay[1,4] <- 'Empty'
  }
  
  track <- AnnotationTrack(genome=gen,chromosome=chr,strand=desiredRegionDisplay[,4],start=desiredRegionDisplay[,2],
                           end=desiredRegionDisplay[,3],
                           feature=desiredRegionDisplay[,5], group=desiredRegionDisplay[,7],
                           id=desiredRegionDisplay[,7], name = "eQTL", groupAnnotation = "group",
                           just.group = just_group, stacking=type_stacking,
                           col.line = "black", col = NULL, collapse= FALSE,showId=showId)
  
  displayPars(track) <- list("SNP_pheno" = "#A6CEE3", "SNP" = "#1F78B4","exon" = "#33A02C",
                             "exon_pheno" = "#B2DF8A", "mRNA" = "#E31A1C",
                             "mRNA_pheno" = "#FB9A99",
                             "cis_local_eQTL" = "#FDBF6F",
                             "trans_local_eQTL" = "#FF7F00", "distal_eQTL" = "#CAB2D6",
                             "cis_local_eQTL_pheno" = "6A3D9A",
                             "trans_local_eQTL_pheno" = "#FFFF99",
                             "distal_eQTL_pheno" = "#B15928","Empty" = "#ffffff")
  
  track
}

#---------------------   GTEx data -------------------------
#------------------- eQTL visualisation data --------------
eQTL_GTEx <- function(gen= "hg19",chr,start, end, bedFilePath, featureDisplay = 'all', showId=FALSE, type_stacking="squish",just_group="above" ) {
  if(is.null(gen)){
    stop("Invalid function eQTL_GTEx :gen null:\n")
  }
  if(is.null(chr)){
    stop("Invalid function eQTL_GTEx :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid function eQTL_GTEx :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid function eQTL_GTEx :end null:\n")
  }
  if(is.null(bedFilePath)){
    stop("Invalid function eQTL_GTEx :bedFilePath null:\n")
  }
  if(is.null(featureDisplay)){
    stop("Invalid function eQTL_GTEx :featureDisplay null:\n")
  }
  
  chrEnsembl <- chrUCSC2ENSEMBL(tolower(chr))
  
  bedFile <- read.table(bedFilePath,header = TRUE,sep="\t")
  desiredRegion_snpgenes <- subset(bedFile, (snp_pos >= start & snp_pos <= end &  snp_chrom == chrEnsembl) 
                                   | ((gene_start >= start | gene_stop <= end) &  gene_chr == chrEnsembl) ,
                                   select= c("snp","snp_pos","snp_chrom","rs_id_dbSNP142_GRCh37p13",
                                             "gene_name", "gene_chr","gene_start","gene_stop","orientation"))
  
  desiredRegionDisplay <- NULL
  if(nrow(desiredRegion_snpgenes) >0) {
    desiredRegion_snpgenes$name <- paste(desiredRegion_snpgenes$snp,desiredRegion_snpgenes$gene_name,sep="_")
    
    snp_desiredRegion <- NULL
    snp_desiredRegion <- desiredRegion_snpgenes[,c("rs_id_dbSNP142_GRCh37p13","snp_chrom","snp_pos","snp_pos")]
    snp_desiredRegion$strand <- "*"
    snp_desiredRegion$type <- "SNP"
    snp_desiredRegion$group <- "cis_local_eQTL"
    snp_desiredRegion$name <- desiredRegion_snpgenes[,c("name")]
    colnames(snp_desiredRegion)<- c("Name_feature","chr","start","end","strand","feature_type","feature_group","name")
    
    gene_desiredRegion <- NULL
    gene_desiredRegion <- desiredRegion_snpgenes[,c("gene_name","gene_chr","gene_start","gene_stop")]
    gene_desiredRegion$strand <- desiredRegion_snpgenes[,c("orientation")]
    gene_desiredRegion$type <- "gene"
    gene_desiredRegion$group <- "cis_local_eQTL"
    gene_desiredRegion$name <- desiredRegion_snpgenes[,c("name")]
    colnames(gene_desiredRegion)<- c("Name_feature","chr","start","end","strand","feature_type","feature_group","name")
    
    desiredRegion <- rbind(snp_desiredRegion,gene_desiredRegion)
    
    desiredRegionDisplay <- desiredRegion
  } 
  
  if( !("all" %in% featureDisplay) ) {
    desiredRegionDisplay <- desiredRegion[which(desiredRegion$feature_type %in% featureDisplay),]
  }
  
  if (is.null(desiredRegionDisplay) ){
    desiredRegionDisplay <- data.frame(nrow=1,ncol=8)
    desiredRegionDisplay[1,1] <- 'Empty'
    desiredRegionDisplay[1,2] <- chr
    desiredRegionDisplay[1,3] <- start
    desiredRegionDisplay[1,4] <- end
    desiredRegionDisplay[1,5] <- '*'
    desiredRegionDisplay[1,6] <- 'Empty'
    desiredRegionDisplay[1,7] <- 'Empty'
    desiredRegionDisplay[1,8] <- 'Empty'
  } else if(nrow(desiredRegionDisplay) == 0){
    desiredRegionDisplay <- data.frame(nrow=1,ncol=8)
    desiredRegionDisplay[1,1] <- 'Empty'
    desiredRegionDisplay[1,2] <- chr
    desiredRegionDisplay[1,3] <- start
    desiredRegionDisplay[1,4] <- end
    desiredRegionDisplay[1,5] <- '*'
    desiredRegionDisplay[1,6] <- 'Empty'
    desiredRegionDisplay[1,7] <- 'Empty'
    desiredRegionDisplay[1,8] <- 'Empty'
  }
  
  track <- AnnotationTrack(genome=gen, chromosome=chrEnsembl,strand=desiredRegionDisplay[,5],start=desiredRegionDisplay[,3],
                           end=desiredRegionDisplay[,4],
                           feature=desiredRegionDisplay[,6], group=desiredRegionDisplay[,8],
                           id=desiredRegionDisplay[,8], name = "eQTL GTEX", groupAnnotation = "group",
                           just.group = just_group , stacking=type_stacking,
                           col.line = "black", col = NULL, collapse= FALSE,showId=showId)
  
  displayPars(track) <- list("SNP" = "#A6CEE3", "CNV" = "#1F78B4",
                             "other_genetic" = "#33A02C","exon" = "#B2DF8A", 
                             "gene" = "#E31A1C", "mRNA" = "#FB9A99",
                             "cis_local_eQTL" = "#FDBF6F",
                             "trans_local_eQTL" = "#FF7F00", "distal_eQTL" = "#CAB2D6",
                             "cis_local_eQTL_pheno" = "6A3D9A",
                             "trans_local_eQTL_pheno" = "#FFFF99",
                             "distal_eQTL_pheno" = "#B15928","Empty" = "#ffffff")
  
  track
}

#------------------- psiQTL visualisation data --------------
psiQTL_GTEx <- function(gen,chr,start, end, bedFilePath, featureDisplay = 'all', showId=FALSE, type_stacking="squish",just_group="above" ) {
  if(is.null(gen)){
    stop("Invalid function psiQTL_GTEx :gen null:\n")
  }
  if(is.null(chr)){
    stop("Invalid function psiQTL_GTEx :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid function psiQTL_GTEx :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid function psiQTL_GTEx :end null:\n")
  }
  
  if(is.null(bedFilePath)){
    stop("Invalid function psiQTL_GTEx :bedFilePath null:\n")
  }
  
  chrEnsembl <- chrUCSC2ENSEMBL(tolower(chr))
  
  bedFile <- read.table(bedFilePath,header=FALSE,sep="\t")
  colnames(bedFile) <- c("snpID","snpID2","gene_name","chr_start_end","chr_snp","chr_exon","pos_snp","pos_middle_exon",
                         "dist_snp_exon","spearman_corr","pvalue","log_pvalue")
  desiredRegion_snpgenes <- subset(bedFile, (pos_snp >= start & pos_snp <= end &  chr_snp == chrEnsembl) 
                                   | ((pos_middle_exon >= start | pos_middle_exon <= end) &  chr_exon == chrEnsembl))
  desiredRegion_snpgenes$name <- paste(desiredRegion_snpgenes$snpID,desiredRegion_snpgenes$gene_name,sep="_")
  
  desiredRegion_snpgenes$type_group <-  "cis_local_psiQTL"
  if(length(which( desiredRegion_snpgenes$chr_snp != desiredRegion_snpgenes$chr_exon)) > 0) {
    desiredRegion_snpgenes[which( desiredRegion_snpgenes$chr_snp != desiredRegion_snpgenes$chr_exon), "type_group" ]<- "distal_psiQTL"
  } else {
    desiredRegion_snpgenes[which( abs(desiredRegion_snpgenes$chr_snp - desiredRegion_snpgenes$chr_exon ) > 5000 ), "type_group" ]<- "distal_psiQTL"
  }
  
  snp_desiredRegion <- desiredRegion_snpgenes[,c("snpID","chr_snp","pos_snp")]
  snp_desiredRegion$type <- "SNP"
  snp_desiredRegion$group <- desiredRegion_snpgenes[,c("type_group")]
  snp_desiredRegion$name <- desiredRegion_snpgenes[,c("name")]
  colnames(snp_desiredRegion)<- c("Name_feature","chr","start","feature_type","feature_group","name")
  
  psi_desiredRegion <- desiredRegion_snpgenes[,c("gene_name","chr_exon","pos_middle_exon")]
  psi_desiredRegion$type <- "exon"
  psi_desiredRegion$group <- desiredRegion_snpgenes[,c("type_group")]
  psi_desiredRegion$name <- desiredRegion_snpgenes[,c("name")]
  colnames(psi_desiredRegion)<- c("Name_feature","chr","start","feature_type","feature_group","name")
  
  desiredRegion <- rbind(snp_desiredRegion,psi_desiredRegion)
  
  desiredRegionDisplay <- desiredRegion
  if( !("all" %in% featureDisplay) ) {
    desiredRegionDisplay <- desiredRegion[which(desiredRegion$feature_type %in% featureDisplay),]
  }
  if (nrow(desiredRegionDisplay) == 0){
    desiredRegionDisplay <- data.frame(nrow=1,ncol=6)
    desiredRegionDisplay[1,1] <- 'Empty'
    desiredRegionDisplay[1,2] <- chr
    desiredRegionDisplay[1,3] <- start
    desiredRegionDisplay[1,4] <- 'Empty'
    desiredRegionDisplay[1,5] <- 'Empty'
    desiredRegionDisplay[1,6] <- 'Empty'
  }
  
  track <- AnnotationTrack(genome=gen,chromosome=chrEnsembl,strand="*",start=desiredRegionDisplay[,3],
                           end=desiredRegionDisplay[,3],
                           feature=desiredRegionDisplay[,4], group=desiredRegionDisplay[,6],
                           id=desiredRegionDisplay[,6], name = "psiQTL GTEX", groupAnnotation = "group",
                           just.group = just_group , stacking=type_stacking,
                           col.line = "black", col = NULL, collapse= FALSE,showId=showId)
  
  displayPars(track) <- list("SNP" = "#A6CEE3", "exon" = "#FB9A99", 
                             "cis_local_psiQTL" = "#FDBF6F",
                             "trans_local_psiQTL" = "#FF7F00", 
                             "distal_psiQTL" = "#CAB2D6",
                             "Empty" = "#ffffff")
  
  track
}


#------------------- sQTL visualisation data --------------
# #need to improve
# sQTL_Altrans_GTEx <- function(gen,chr,start, end, bedFilePath, featureDisplay = 'all', showId=FALSE, type_stacking="squish",just_group="above" ) {
#   if(is.null(gen)){
#     stop("Invalid function sQTL_Altrans_GTEx :gen null:\n")
#   }
#   if(is.null(chr)){
#     stop("Invalid function sQTL_Altrans_GTEx :chr null:\n")
#   }
#   if(is.null(start)){
#     stop("Invalid function sQTL_Altrans_GTEx :start null:\n")
#   }
#   if(is.null(end)){
#     stop("Invalid function sQTL_Altrans_GTEx :end null:\n")
#   }
#   
#   if(is.null(bedFilePath)){
#     stop("Invalid function sQTL_Altrans_GTEx :bedFilePath null:\n")
#   }
#   
#   chrEnsembl <- chrUCSC2ENSEMBL(tolower(chr))
#   
#   bedFile <- read.table(bedFilePath,header = TRUE,sep="\t")
#   
#   gtfFile <- read.table(gtfFilePath,sep="\t")
#   
#   attribu_gene <- gtfFile[which( grepl("ENSG00000001036", gtfFile[,9]) & gtfFile[,3] == "transcript"),1]
#   
#   desiredRegion_snpgenes <- subset(bedFile, (snp_pos >= start & snp_pos <= end &  snp_chrom == chrEnsembl) 
#                                    | ((gene_start >= start | gene_stop <= end) &  gene_chr == chrEnsembl) ,
#                                    select= c("snp","snp_pos","snp_chrom","rs_id_dbSNP142_GRCh37p13",
#                                              "gene_name", "gene_chr","gene_start","gene_stop","orientation"))
#   desiredRegion_snpgenes$name <- paste(desiredRegion_snpgenes$snp,desiredRegion_snpgenes$gene_name,sep="_")
#   
#   snp_desiredRegion <- desiredRegion_snpgenes[,c("rs_id_dbSNP142_GRCh37p13","snp_chrom","snp_pos","snp_pos")]
#   snp_desiredRegion$strand <- "*"
#   snp_desiredRegion$type <- "SNP"
#   snp_desiredRegion$group <- "cis_local_eQTL"
#   snp_desiredRegion$name <- desiredRegion_snpgenes[,c("name")]
#   colnames(snp_desiredRegion)<- c("Name_feature","chr","start","end","strand","feature_type","feature_group","name")
#   
#   gene_desiredRegion <- desiredRegion_snpgenes[,c("gene_name","gene_chr","gene_start","gene_stop")]
#   gene_desiredRegion$strand <- desiredRegion_snpgenes[,c("orientation")]
#   gene_desiredRegion$type <- "gene"
#   gene_desiredRegion$group <- "cis_local_eQTL"
#   gene_desiredRegion$name <- desiredRegion_snpgenes[,c("name")]
#   colnames(gene_desiredRegion)<- c("Name_feature","chr","start","end","strand","feature_type","feature_group","name")
#   
#   desiredRegion <- rbind(snp_desiredRegion,gene_desiredRegion)
#   
#   desiredRegionDisplay <- desiredRegion
#   if( !("all" %in% featureDisplay) ) {
#     desiredRegionDisplay <- desiredRegion[which(desiredRegion$feature_type %in% featureDisplay),]
#   }
#   if (nrow(desiredRegionDisplay) == 0){
#     desiredRegionDisplay <- data.frame(nrow=1,ncol=8)
#     desiredRegionDisplay[1,1] <- 'Empty'
#     desiredRegionDisplay[1,2] <- chr
#     desiredRegionDisplay[1,3] <- start
#     desiredRegionDisplay[1,4] <- end
#     desiredRegionDisplay[1,5] <- '*'
#     desiredRegionDisplay[1,6] <- 'Empty'
#     desiredRegionDisplay[1,7] <- 'Empty'
#     desiredRegionDisplay[1,8] <- 'Empty'
#   }
#   
#   track <- AnnotationTrack(genome=gen,chromosome=chrEnsembl,strand=desiredRegionDisplay[,5],start=desiredRegionDisplay[,3],
#                            end=desiredRegionDisplay[,4],
#                            feature=desiredRegionDisplay[,6], group=desiredRegionDisplay[,8],
#                            id=desiredRegionDisplay[,8], name = "eQTL GTEX", groupAnnotation = "group",
#                            just.group = just_group , stacking=type_stacking,
#                            col.line = "black", col = NULL, collapse= FALSE,showId=showId)
#   
#   displayPars(track) <- list("SNP" = "#A6CEE3","exon" = "#B2DF8A","Empty" = "#ffffff")
#   
#   track
# }

#------------------- Visualisation of geneExpression data --------------
geneExpression_GTEx <- function(chr,start, end, gtfFilePath, genexpressionFilePath, tissue="" ) {
  if(is.null(chr)){
    stop("Invalid function geneExpression_GTEx :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid function geneExpression_GTEx :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid function geneExpression_GTEx :end null:\n")
  }
  
  if(is.null(gtfFilePath)){
    stop("Invalid function geneExpression_GTEx :gtfFilePath null:\n")
  }
  
  if(is.null(genexpressionFilePath)){
    stop("Invalid function geneExpression_GTEx :genexpressionFilePath null:\n")
  }
  
  chrEnsembl <- chrUCSC2ENSEMBL(tolower(chr))
  
  geneExp <- read.table(genexpressionFilePath,header = TRUE,sep="\t")
  
  gffRangedData<-import.gff(gtfFilePath)
  myGranges<-as(gffRangedData, "GRanges")
  
  if (nrow(desiredRegionDisplay) == 0){
    desiredRegionDisplay <- data.frame(nrow=1,ncol=8)
    desiredRegionDisplay[1,1] <- 'Empty'
    desiredRegionDisplay[1,2] <- chr
    desiredRegionDisplay[1,3] <- start
    desiredRegionDisplay[1,4] <- end
    desiredRegionDisplay[1,5] <- '*'
    desiredRegionDisplay[1,6] <- 'Empty'
    desiredRegionDisplay[1,7] <- 'Empty'
    desiredRegionDisplay[1,8] <- 'Empty'
  }
  
  track <- DataTrack(chromosome=chr,strand=desiredRegionDisplay[,5],start=desiredRegionDisplay[,3],
                     end=desiredRegionDisplay[,4],
                     feature=desiredRegionDisplay[,6])
  
  track
}

#------------------- Visualisation of geneExpression data --------------
imprintedGenes_GTEx <- function(gen="hg19",chr,start, end, tissues="all", classification="all",showId=FALSE) {
  if(is.null(gen)){
    stop("Invalid function ImprintedGenes_GTEx :gen null:\n")
  }
  if(is.null(chr)){
    stop("Invalid function ImprintedGenes_GTEx :chr null:\n")
  }
  if(is.null(start)){
    stop("Invalid function ImprintedGenes_GTEx :start null:\n")
  }
  if(is.null(end)){
    stop("Invalid function ImprintedGenes_GTEx :end null:\n")
  }
  
  if(is.null(tissues)){
    stop("Invalid function ImprintedGenes_GTEx :tissue null:\n")
  }
  if(is.null(classification)){
    stop("Invalid function ImprintedGenes_GTEx :classification null:\n")
  }
  
  chrEnsembl <- chrUCSC2ENSEMBL(tolower(chr))
  geneTrunk="hg19"
  
  imprintedGenesGTEx <- NULL
  data(imprintedGenesGTEx)
  if ( !(exists("imprintedGenesGTEx") && is.data.frame(get("imprintedGenesGTEx"))) ) {
    extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
    ImprintedGenesGTExFile <- file.path(extdata, "/GTEX/GTEx_imprinted_genes_v3_2015.csv")
    imprintedGenesGTEx <- read.table(ImprintedGenesGTExFile,header = TRUE,sep=",")
  }
  
  genesTrack <- genes_ENSEMBL(gen,chr,start,end)
  
  desiredImprintedGenesGTEx <- imprintedGenesGTEx
  if( !("all" %in% tissues) ) {
    desiredImprintedGenesGTEx <- desiredImprintedGenesGTEx[which(desiredImprintedGenesGTEx$Tissue.Name %in% tissues),]
  } 
  if( !("all" %in% classification) ) {
    desiredImprintedGenesGTEx <- desiredImprintedGenesGTEx[which(desiredImprintedGenesGTEx$Classification %in% classification),]
  }
  
  biomTrackDisplay <- NULL
  
  if(length(desiredImprintedGenesGTEx) > 0 ) {
    desiredImprintedGenes <- unique(desiredImprintedGenesGTEx$Ensembl.ID)
    geneTrackShow <- genesTrack[gene(genesTrack) %in% desiredImprintedGenes,]
    
    if (length(feature(geneTrackShow)) == 0){
      biomTrackDisplay <- matrix(data = NA, nrow = 1, ncol = 7)
      biomTrackDisplay[1,1] <- 'Empty'
      biomTrackDisplay[1,2] <- chrEnsembl
      biomTrackDisplay[1,3] <- start
      biomTrackDisplay[1,4] <- end
      biomTrackDisplay[1,5] <- '*'
      biomTrackDisplay[1,6] <- ' '
      biomTrackDisplay[1,7] <- 'Empty'
    } else {
      genesList <- gene(geneTrackShow) 
      chromList <- rep(chromosome(geneTrackShow),length(genesList))
      startList <- start(geneTrackShow) 
      endList <- end(geneTrackShow) 
      strandList <- strand(geneTrackShow) 
      symbolList <- symbol(geneTrackShow)
      featuresList <- cbind(genesList,chromList,startList,endList,strandList,symbolList)
      targets <- desiredImprintedGenesGTEx[,c(1,4)]
      colnames(targets)[1]<-"genesList"
      biomTrackDisplay <- merge(featuresList,targets,by="genesList",all.x=TRUE)
    }
    
  } else {
    biomTrackDisplay <- matrix(data = NA, nrow = 1, ncol = 7)
    biomTrackDisplay[1,1] <- 'Empty'
    biomTrackDisplay[1,2] <- chrEnsembl
    biomTrackDisplay[1,3] <- start
    biomTrackDisplay[1,4] <- end
    biomTrackDisplay[1,5] <- '*'
    biomTrackDisplay[1,6] <- ' '
    biomTrackDisplay[1,7] <- 'Empty'
  }
  
  geneTrackShow <- AnnotationTrack(genome=gen,chromosome=chrEnsembl,strand=biomTrackDisplay[,5],
                                   start=biomTrackDisplay[,3],end=biomTrackDisplay[,4],
                                   feature=biomTrackDisplay[,7],group=biomTrackDisplay[,6],
                                   id=biomTrackDisplay[,6], name = "Imprinted genes GTEx",stacking="squish", 
                                   just.group="above",
                                   col.line = "black", col = NULL, collapse= FALSE,showId=showId)
  
  displayPars(geneTrackShow) <- list("consistent with biallelic" = "#1b9e77", "NC" = "#7570b3", 
                                     "imprinted" =  "#e7298a","consistent with imprinting" = "#66a61e", 
                                     "biallelic" = "#e6ab02","Empty" = "#ffffff")
  
  geneTrackShow
  
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
