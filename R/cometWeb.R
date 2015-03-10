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
# File: cometWeb.R
# Author: Tiphaine Martin
# Email: tiphaine.martin@kcl.ac.uk
# Purpose: coMET allows the display of p-values from association
#           with a correlation heatmap.
# Version : 0.99.9
###########################################################################

#coMET allows the display of p-values from association with a correlation heatmap. 
#This version coMET.web is used for web site.

# The values given such as attributes of function erase ones of configuration file.
comet.web <- function(mydata.file = NULL,
                      mydata.format = c("site","region","site_asso","region_asso"),
                      mydata.large.file = NULL,
                      mydata.large.format = c("site","region","site_asso","region_asso"),
                      cormatrix.file = NULL,
                      cormatrix.method = c("spearman","pearson","kendall"),
                      cormatrix.format= c("cormatrix","raw","raw_rev"),
                      cormatrix.color.scheme = "heat",
                      cormatrix.conf.level=0.05,
                      cormatrix.sig.level= 1,
                      cormatrix.adjust="none",
                      mydata.ref = NULL,
                      genome="hg19",
                      start = NULL,
                      end = NULL,
                      zoom = FALSE,
                      lab.Y = "log",
                      pval.threshold = 10e-8,
                      disp.pval.threshold = 1,
                      disp.association = FALSE,
                      disp.association.large = FALSE,
                      disp.region = FALSE,
                      disp.region.large = FALSE, 
                      symbols = "circle-fill",
                      symbols.large = NA,
                      sample.labels = NULL,
                      sample.labels.large = NULL,
                      use.colors = TRUE,
                      disp.color.ref = TRUE,
                      color.list = NULL,
                      color.list.large = NULL,
                      biofeat.user.file= NULL,
                      biofeat.user.type= c("GeneRegion","Annotation","Data"),
                      biofeat.user.type.plot = NULL,
                      list.tracks = "geneENSEMBL,CGI,ChromHMM,DNAse,RegENSEMBL,SNP",
                      pattern.regulation="GM12878",
                      image.title = NULL,
                      image.name = "coMET",
                      image.type = c("pdf", "eps"),
                      image.size = 3.5,
                      print.image = FALSE,
                      config.file = NULL,
                      verbose = FALSE) {
  
  #-------------------MAIN FUNCTION----------------------------
  #LIST.TRACK DEFAULT VALUE geneENSEMBL,CGI,ChromHMM,DNAse,RegENSEMBL,SNP
  #DEBUG STATEMENT
  if (verbose) cat("START COMET  VERSION WEB\n")
  
  #-------------------GLOBAL VARIABLES BEGINS------------------
  
  #create three variables (global variables), which are lists of all other variables
  #used throughout the program to be passed to functions with values changed as
  #necessary.
  gbl.var <- NULL
  
  #DATA VARIABLES
  presence.mydata<- 0 # 0 presence of omic-WAS results (object mydata), 1 no presence of data
  general.data <- NULL
  split.biofeature.data.user.file <- NULL
  split.biofeature.data.user.type <- NULL
  split.biofeature.data.user.type.plot <- NULL
  
  mydata.data <- NULL
  mydata.num <- 0
  mydata.names <- NULL
  mydata.chr <- NULL
  mydata.gen <- NULL
  mydata.ref.pos <- NULL
  mydata.best.position <- NULL
  
  sorted.mydata.names <- NULL
  sorted.mydata.pos <- NULL
  sorted.mydata.pos.zoom <- NULL
  sorted.mydata.large.pos.zoom <- NULL
  
  cur.sample <- NULL
  mydata.samples <- NULL
  
  cur.sample.large <- NULL
  mydata.samples.large <- NULL
  
  mydata.hash.names.pos <- new.env(hash=TRUE)
  mydata.hash.names.start <- new.env(hash=TRUE)
  mydata.hash.names.end <- new.env(hash=TRUE)
  mydata.hash.pos.names <- new.env(hash=TRUE)
  mydata.large.hash.pos.names <- new.env(hash=TRUE)
  mydata.large.hash.names.start <- new.env(hash=TRUE)
  mydata.large.hash.names.end <- new.env(hash=TRUE)
  mydata.large.hash.names.pos <- new.env(hash=TRUE)
  
  split.mydata.file <- NULL
  split.cormatrix.file <- NULL
  cormatrix.pvalue.data <- NULL
  
  #GRAPH BOUNDARY VARIABLES
  
  pval.flag <- TRUE
  
  cur.exp <- 0
  
  total.dist <- NULL
  min.dist <- NULL
  max.dist <- NULL
  
  min.x <- NULL
  min.user.x <- NULL
  min.y <- NULL
  
  max.x <- NULL
  max.user.x <- NULL
  max.y <- NULL
  
  r.x <- NULL
  r.y <- NULL
  
  axis.x <- NULL
  axis.y <- NULL
  
  equidis.pos <- NULL
  
  #FORMATTING VARIABLES
  
  split.sample.labels <- NULL
  split.color.list <- NULL
  color.list <- NULL
  symbol.list <- NULL
  fill.list <- NULL
  split.format <- NULL
  split.association <- NULL
  split.region <- NULL
  split.type <- NULL
  
  large.split.sample.labels <- NULL
  large.split.color.list <- NULL
  large.color.list <- NULL
  large.symbol.list <- NULL
  large.fill.list <- NULL
  large.split.format <- NULL
  large.split.association <- NULL
  large.split.region <- NULL
  large.split.type <- NULL
  ref <- NULL
  samples <- NULL
  large.samples <- NULL
  
  
  palette.size <- NULL
  cex.factor <- 0.5
  cex.factor.symbol <- 0.25
  font.size <- NULL
  line.width <- 0.5
  
  cormatrix.data <- NULL
  split.cormatrix.type <- NULL
  matrix.data <- NULL
  cormatrix.data.full <- NULL
  cormatrix.pvalue.data <- NULL
  cormatrix.pvalue.data.full <- NULL
  cormatrix.CI.data <- NULL
   
  mySession <- NULL
  listtracks_gviz <- NULL
  listtracks_user <- NULL
  listtracks_ggbio <- NULL
  listtracks_trackviewer <- NULL
  split.list.tracks <- NULL
  
  gbl.var <- list(mydata.data = mydata.data,
                  general.data = general.data,
                  presence.mydata=presence.mydata,
                  split.mydata.file =split.mydata.file,
                  split.type = split.type,
                  split.biofeature.data.user.file = split.biofeature.data.user.file,
                  split.biofeature.data.user.type = split.biofeature.data.user.type,
                  split.biofeature.data.user.type.plot = split.biofeature.data.user.type.plot,
                  split.format = split.format,
                  split.association =  split.association,
                  split.region = split.region,
                  large.split.format = large.split.format,
                  large.split.association =  large.split.association,
                  large.split.region = large.split.region,
                  large.split.type = large.split.type,
                  mySession = mySession,
                  listtracks_gviz = listtracks_gviz,
                  listtracks_ggbio = listtracks_ggbio,
                  listtracks_trackviewer = listtracks_trackviewer,
                  split.list.tracks = split.list.tracks,
                  mydata.num = mydata.num,
                  mydata.names = mydata.names,
                  mydata.chr = mydata.chr,
                  mydata.gen = mydata.gen,
                  mydata.best.position = mydata.best.position,
                  sorted.mydata.names = sorted.mydata.names,
                  sorted.mydata.pos = sorted.mydata.pos,
                  sorted.mydata.pos.zoom = sorted.mydata.pos.zoom,
                  mydata.ref.pos = mydata.ref.pos,
                  cur.sample = cur.sample,
                  mydata.samples = mydata.samples,
                  cur.sample.large = cur.sample.large,
                  mydata.samples.large = mydata.samples.large,
                  mydata.hash.names.pos = mydata.hash.names.pos,
                  mydata.hash.names.start = mydata.hash.names.start,
                  mydata.hash.names.end = mydata.hash.names.end,
                  mydata.large.hash.names.pos = mydata.large.hash.names.pos,
                  mydata.large.hash.names.start = mydata.large.hash.names.start,
                  mydata.large.hash.names.end = mydata.large.hash.names.end,
                  mydata.hash.pos.names = mydata.hash.pos.names,
                  cormatrix.data = cormatrix.data,
                  split.cormatrix.type = split.cormatrix.type,
                  cormatrix.pvalue.data = cormatrix.pvalue.data,
                  cormatrix.pvalue.data.full = cormatrix.pvalue.data.full,
                  split.cormatrix.file = split.cormatrix.file,
                  matrix.data = matrix.data,
                  cormatrix.data.full = cormatrix.data.full,
                  cormatrix.data.full = cormatrix.data.full,
                  cormatrix.CI.data = cormatrix.CI.data,
                  samples = samples,
                  large.samples = large.samples,
                  pval.flag = pval.flag,
                  cur.exp = cur.exp,
                  total.dist = total.dist,
                  min.dist = min.dist,
                  max.dist = max.dist,
                  min.x = min.x,
                  min.user.x = min.user.x,
                  min.y = min.y,
                  max.x = max.x,
                  max.user.x = max.user.x,
                  max.y = max.y,
                  r.x = r.x,
                  r.y = r.y,
                  axis.x = axis.x,
                  axis.y = axis.y,
                  equidis.pos = equidis.pos,
                  split.sample.labels = split.sample.labels,
                  split.color.list = split.color.list,
                  color.list = color.list,
                  symbol.list = symbol.list,
                  fill.list = fill.list,
                  large.split.sample.labels = large.split.sample.labels,
                  large.split.color.list = large.split.color.list,
                  large.color.list = large.color.list,
                  large.symbol.list = large.symbol.list,
                  large.fill.list = large.fill.list,
                  palette.size = palette.size,
                  cex.factor = cex.factor,
                  cex.factor.symbol = cex.factor.symbol,
                  font.size = font.size,
                  line.width = line.width,
                  verbose = verbose
  )
  
  #-------------------GLOBAL VARIABLES ENDS------------------
  
  #-------------------CONFIGURATION VARIABLES BEGINS---------
  mydata.type = "file"
  mydata.large.type = "listfile"
  cormatrix.type = "listfile"
  disp.cormatrixmap = TRUE
  disp.pvalueplot = TRUE
  disp.mydata.names = TRUE
  disp.connecting.lines = TRUE
  disp.mydata = TRUE
  disp.type = "symbol"
  biofeat.user.type.plot="histogram"
  tracks.gviz = NULL
  tracks.ggbio = NULL
  tracks.trackviewer = NULL
  biofeat.user.file= NULL
  palette.file = NULL
  disp.color.bar = TRUE
  disp.phys.dist = TRUE
  disp.legend = TRUE
  disp.marker.lines = TRUE
  disp.mult.lab.X = FALSE
  connecting.lines.factor = 1.5
  connecting.lines.adj = 0.01
  connecting.lines.vert.adj = -1
  connecting.lines.flex = 0
  font.factor = NULL
  color.list = "red"
  symbol.factor = NULL
  dataset.gene = "hsapiens_gene_ensembl"
  DATASET.SNP="hsapiens_snp"
  VERSION.DBSNP="snp142Common"
  DATASET.SNP.STOMA="hsapiens_snp_som"
  DATASET.REGULATION="hsapiens_feature_set"
  DATASET.STRU="hsapiens_structvar"
  DATASET.STRU.STOMA="hsapiens_structvar_som"
  BROWSER.SESSION="UCSC"
  
  
  #-------------------UPDATE CONFIGURATION VARIABLES---------
  config.var <- list(mydata.file = mydata.file,
                     mydata.format = mydata.format,
                     mydata.ref = mydata.ref,
                     mydata.type = mydata.type,
                     mydata.large.file = mydata.large.file,
                     mydata.large.format = mydata.large.format,
                     mydata.large.type = mydata.large.type,
                     biofeat.user.file = biofeat.user.file,
                     biofeat.user.type = biofeat.user.type,
                     biofeat.user.type.plot = biofeat.user.type.plot,
                     genome = genome,
                     palette.file = palette.file,
                     lab.Y = lab.Y,
                     disp.pval.threshold = disp.pval.threshold,
                     start = start,
                     end = end,
                     zoom = zoom,
                     disp.association = disp.association,
                     disp.association.large = disp.association.large,
                     disp.region = disp.region,
                     disp.region.large = disp.region.large,
                     dataset.gene = dataset.gene,
                     DATASET.SNP = DATASET.SNP,
                     DATASET.SNP.STOMA = DATASET.SNP.STOMA,
                     DATASET.REGULATION = DATASET.REGULATION,
                     DATASET.STRU = DATASET.STRU,
                     DATASET.STRU.STOMA = DATASET.STRU.STOMA,
                     VERSION.DBSNP = VERSION.DBSNP,
                     pattern.regulation = pattern.regulation, 
                     BROWSER.SESSION = BROWSER.SESSION,
                     tracks.gviz = tracks.gviz,
                     tracks.ggbio = tracks.ggbio,
                     list.tracks = list.tracks,
                     tracks.trackviewer = tracks.trackviewer,
                     symbols = symbols,
                     symbols.large = symbols.large,
                     pval.threshold = pval.threshold,
                     disp.type = disp.type,
                     disp.cormatrixmap = disp.cormatrixmap,
                     disp.mydata = disp.mydata,
                     disp.marker.lines = disp.marker.lines,
                     disp.connecting.lines = disp.connecting.lines,
                     disp.mydata.names = disp.mydata.names,
                     disp.pvalueplot = disp.pvalueplot,
                     cormatrix.color.scheme = cormatrix.color.scheme,
                     color.list = color.list,
                     color.list.large = color.list.large,
                     use.colors = use.colors,
                     disp.color.ref = disp.color.ref,
                     cormatrix.method = cormatrix.method,
                     cormatrix.file = cormatrix.file,
                     cormatrix.type = cormatrix.type,
                     cormatrix.conf.level=cormatrix.conf.level,
                     cormatrix.sig.level= cormatrix.sig.level,
                     cormatrix.adjust=cormatrix.adjust,
                     disp.color.bar = disp.color.bar,
                     disp.phys.dist = disp.phys.dist,
                     image.title = image.title,
                     disp.legend = disp.legend,
                     sample.labels = sample.labels,
                     sample.labels.large = sample.labels.large,
                     image.type = image.type,
                     image.size = image.size,
                     disp.mult.lab.X = disp.mult.lab.X,
                     image.name = image.name,
                     print.image = print.image,
                     connecting.lines.factor = connecting.lines.factor,
                     connecting.lines.adj = connecting.lines.adj,
                     connecting.lines.vert.adj = connecting.lines.vert.adj,
                     connecting.lines.flex = connecting.lines.flex,
                     font.factor = font.factor,
                     symbol.factor = symbol.factor,
                     verbose = verbose)
  
  #-------------------CONFIGURATION VARIABLES BEGINS from config file---------
  
  gbl.var$verbose <- config.var$verbose
  
  if(!is.null(config.file)) {
    config.var <- read.config(config.file, config.var)
  }
  
  if(!is.null(mydata.file) ){
    config.var$mydata.file <- mydata.file
  }
  if(!is.null(mydata.format) & length(mydata.format) == 1){
    config.var$mydata.format <- mydata.format
  }
  if(!is.null(mydata.large.format) & length(mydata.large.format) == 1){
    config.var$mydata.large.format <- mydata.large.format
  }
  if(!is.null(cormatrix.file)){
    config.var$cormatrix.file <- cormatrix.file
  }
  if(!is.null(cormatrix.format) & length(cormatrix.format) == 1){
    config.var$cormatrix.format <- cormatrix.format
  }
  if(!is.null(cormatrix.method) & length(cormatrix.method) == 1){
    config.var$cormatrix.method <- cormatrix.method
  }
  
  if(!is.null(biofeat.user.file)){
    config.var$biofeat.user.file <- biofeat.user.file
  }
  if(!is.null(palette.file)){
    config.var$palette.file <- palette.file
  }
  if(!is.null(start)){
    config.var$start <- start
  }
  if(!is.null(end)){
    config.var$end <- end
  }
  if(!is.null(genome)){
    config.var$genome <- genome
  }
  if(is.null(config.var$mydata.file)) {
    stop("Invalid MYDATA data file: ", config.var$mydata.file, "\n")
  }
  
  #-------------- CHECK if all parameters have values -------------------
  check.configVar(config.var)
  
  #------------- connection to database
  if(!is.null(config.var$genome) & !is.null(config.var$BROWSER.SESSION)){
    mySession <- browserSession(config.var$BROWSER.SESSION)
    genome(mySession) <- config.var$genome
    gbl.var$mySession <- mySession
  }
  
  
  #------------- READ DATA for ANNOTATION TRACKS
  if(!is.null(config.var$biofeat.user.file)){
    split.biofeature.data.user.file <- strsplit(config.var$biofeat.user.file, ",")
    split.biofeature.data.user.type <- strsplit(config.var$biofeat.user.type, ",")
    
    gbl.var$split.biofeature.data.user.file <- split.biofeature.data.user.file
    gbl.var$split.biofeature.data.user.type <- split.biofeature.data.user.type
    
    if(!is.null(config.var$biofeat.user.type.plot)){
      split.biofeature.data.user.type.plot <- strsplit(config.var$biofeat.user.type.plot, ",")
      gbl.var$split.biofeature.data.user.type.plot <- split.biofeature.data.user.type.plot
    }
  }
  
  #------------ SELECTION ANNOTATION TRACKS in comet.web
  if(!is.null(config.var$list.tracks)){
    split.list.tracks <- strsplit(config.var$list.tracks, ",")
    
    someenv.list.tracks<-hash()
    for(i in 1:length(split.list.tracks[[1]]))
    {
      info.track <- split.list.tracks[[1]][i]
      # if (config.var$verbose) cat(" long",info.track,"\n")
      someenv.list.tracks[[ info.track ]]<- info.track
    } 
    gbl.var$split.list.tracks <- someenv.list.tracks
  }
  
  
  #------------- READ DATA and UPDATE VARIABLES
  if(!is.null(config.var$cormatrix.file)){
    if(!is.null(config.var$cormatrix.type)){
      split.cormatrix.type <- strsplit(config.var$cormatrix.type, ",")
      comatrix.type.length <- length(split.cormatrix.type[[1]])
      gbl.var$split.cormatrix.type <- split.cormatrix.type
    } else {
      stop("Need to define the format of your correlation matrix from FILE or from MATRIX")
    }
    
    split.cormatrix.file <- strsplit(config.var$cormatrix.file, ",")
    comatrix.file.length <- length(split.cormatrix.file[[1]])
    gbl.var$split.cormatrix.file <- split.cormatrix.file
  } else {
    warning("No visualisation of correlation matrice")
    disp.mydata <- FALSE
  }
  
  if(!is.null(config.var$mydata.file)){
    #--- TYPE of DATA
    if(!is.null(config.var$mydata.type)){
      split.type <- strsplit(config.var$mydata.type, ",")
      split.type.length <- length(split.type[[1]])
      gbl.var$split.type <- split.type
    } else{
      stop("Need to define the type of data: FILE, MATRIX\n")
    }
    
    split.mydata.file <- strsplit(config.var$mydata.file, ",")
    mydata.file.length <- length(split.mydata.file[[1]])
    gbl.var$split.mydata.file <- split.mydata.file
    
    gbl.var$mydata.samples <- mydata.file.length
    gbl.var$hap.samples <- 0
    samples <- mydata.file.length
    gbl.var$samples <- samples
    
    #REFERENCE DATA
    split.color.list <- create.color.list(config.var,gbl.var)
    gbl.var$color.list <- split.color.list[[1]]
    
    split.symbol.fill.lists <- create.symbol.list(config.var, split.color.list, gbl.var)
    
    gbl.var$symbol.list <- split.symbol.fill.lists$split.symbol.list
    gbl.var$fill.list <- split.symbol.fill.lists$split.fill.list[[1]]
    
    
    #create sample labels
    #split sample.labels on comma, then take the top 5 or less entries and create a list
    if(!is.null(config.var$sample.labels)) {
      split.tmp.sample.list <- strsplit(config.var$sample.labels, ",")
    } else {
      split.tmp.sample.list <- rep("sample", samples)
    }
    
    if(length(split.tmp.sample.list[[1]]) < samples) {
      warning("sample.labels contains less labels than the number of samples. Labels must be separated by commas without spaces.\n")
    }
    
    #--- FORMAT of DATA
    if(!is.null(config.var$mydata.format)){
      split.format <- strsplit(config.var$mydata.format, ",")
      split.format.length <- length(split.format[[1]])
      gbl.var$split.format <- split.format
    }
    
    #-- Visualisation of ASSOCIATION
    if(!is.null(config.var$disp.association)){
      split.association <- strsplit(config.var$disp.association, ",")
      split.association.length <- length(split.association[[1]])
      gbl.var$split.association <- split.association
    }
    
    #---VISUALISATION of region
    if(!is.null(config.var$disp.region)){
      split.region <- strsplit(config.var$disp.region, ",")
      split.region.length <- length(split.region[[1]])
      gbl.var$split.region <- split.region
    }
    
    split.tmp.sample.list.truncated <- substr(split.tmp.sample.list[[1]], 1, 25)
  } else {
    stop("No data visualise")
  }
  
  #DEBUG STATEMENT
  #if (config.var$verbose) cat("samples ", samples, "\n")
  
  
  gbl.var$split.sample.labels <- list(head(split.tmp.sample.list.truncated, samples))
  
  #---LARGE REFERENCE DATA
  if(!is.null(config.var$mydata.large.file)){
    #--- TYPE of LARGE DATA
    if(!is.null(config.var$mydata.large.type)){
      large.split.type <- strsplit(config.var$mydata.large.type, ",")
      large.split.type.length <- length(large.split.type[[1]])
      gbl.var$large.split.type <- large.split.type
    }else{
      stop("Need to define the type of extra data (LISTFILE,LISTMATRIX)\n")
    }
    
    large.split.mydata.file <- strsplit(config.var$mydata.large.file, ",")
    large.mydata.file.length <- length(large.split.mydata.file[[1]])
    gbl.var$large.split.mydata.file <- large.split.mydata.file
    
    gbl.var$large.mydata.samples <- large.mydata.file.length
    gbl.var$large.hap.samples <- 0
    large.samples <- large.mydata.file.length
    gbl.var$large.samples <- large.samples
    
    large.split.color.list <- create.color.list.large(config.var,gbl.var)
    gbl.var$large.color.list <- large.split.color.list[[1]]
    
    large.split.symbol.fill.lists <- create.symbol.list.large(config.var, large.split.color.list, gbl.var)
    
    gbl.var$large.symbol.list <- large.split.symbol.fill.lists$large.split.symbol.list
    gbl.var$large.fill.list <- large.split.symbol.fill.lists$large.split.fill.list[[1]]
    
    #create sample labels
    #split sample.labels on comma, then take the top 5 or less entries and create a list
    large.split.tmp.sample.list <- strsplit(config.var$sample.labels.large, ",")
    
    if(length(large.split.tmp.sample.list[[1]]) < large.samples) {
      warning("LARGE.sample.labels contains less labels than the number of samples. Labels must be separated by commas without spaces.\n")
    }
    
    large.split.tmp.sample.list.truncated <- substr(large.split.tmp.sample.list[[1]], 1, 25)
    
    #DEBUG STATEMENT
    #if (config.var$verbose) cat("slarge.amples ", large.samples, "\n") 
    
    gbl.var$large.split.sample.labels <- list(head(large.split.tmp.sample.list.truncated, large.samples))
    
    #--- FORMAT of LARGE DATA
    if(!is.null(config.var$mydata.large.format)){
      #if (config.var$verbose) cat("ASSO", config.var$mydata.large.format,"\n") 
      large.split.format <- strsplit(config.var$mydata.large.format, ",")
      large.split.format.length <- length(large.split.format[[1]])
      gbl.var$large.split.format <- large.split.format
      # if (config.var$verbose) cat("L",large.split.format.length,"  value ", gbl.var$large.split.format[[1]][1],"\n") 
    }
    
    #-- Visualisation of ASSOCIATION
    if(!is.null(config.var$disp.association.large)){
      large.split.association <- strsplit(config.var$disp.association.large, ",")
      large.split.association.length <- length(large.split.association[[1]])
      gbl.var$large.split.association <- large.split.association
    }
    
    #---VISUALISATION of region
    if(!is.null(config.var$disp.region.large)){
      large.split.region <- strsplit(config.var$disp.region.large, ",")
      large.split.region.length <- length(large.split.region[[1]])
      gbl.var$large.split.region <- large.split.region
    }
  } else {
    gbl.var$large.split.mydata.file <- NULL
    gbl.var$large.mydata.samples <- NULL
    gbl.var$large.hap.samples <- 0
    gbl.var$large.samples <- NULL
    gbl.var$large.color.list <- NULL
    gbl.var$large.symbol.list <- NULL
    gbl.var$large.fill.list <- NULL
    gbl.var$large.split.sample.labels <- NULL
    gbl.var$large.split.sample.format <- NULL
    gbl.var$large.split.association <- NULL
    gbl.var$large.split.region <- NULL
  }
  
  #---- FIX GENERAL VARIABLES
  fix.var <- retrieve.data(config.var, gbl.var)
  
  config.var <- fix.var$config.var
  gbl.var <- fix.var$gbl.var
  
  #############DRAW Picture
  
  #---------------- DRAW DIFFERENT ANNOTATION TRACK ---------
  gbl.var <- create.tracks.web(config.var,gbl.var)
  
  #------ DRAW the STRUCTURE COMET	
  if(print.image == FALSE || is.null(config.var$image.name)){
    gbl.var <- draw.plot.comet.web(config.var, gbl.var,newpage=TRUE)
  } else{
    printPlot.comet.web(config.var, gbl.var)
  }
  
  #END FUNCTION LIST
  
  invisible(list(config.var, gbl.var))
  #DEBUG STATEMENT
  if (config.var$verbose) cat("FINISH WEB COMET \n")
}


