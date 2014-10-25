
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
# Version : 0.99.0
###########################################################################

#coMET allows the display of p-values from association with a correlation heatmap. 
#This version coMET.web is used for web site.

# The values given such as attributes of function erase ones of configuration file.
comet.web <- function(MYDATA.FILE = NULL,
                      MYDATA.FORMAT = c("SITE","REGION","SITE_ASSO","REGION_ASSO"),
                      MYDATA.LARGE.FILE = NULL,
                      MYDATA.LARGE.FORMAT = c("SITE","REGION","SITE_ASSO","REGION_ASSO"),
                      CORMATRIX.FILE = NULL,
                      CORMATRIX.METHOD = c("spearman","pearson","kendall"),
                      CORMATRIX.FORMAT= c("CORMATRIX","RAW","RAW_REV"),
                      CORMATRIX.COLOR.SCHEME = "heat",
                      MYDATA.REF = NULL,
                      GENOME="hg19",
                      START = NULL,
                      END = NULL,
                      ZOOM = FALSE,
                      LAB.Y = "log",
                      PVAL.THRESHOLD = 10e-8,
                      DISP.PVAL.THRESHOLD = 1,
                      DISP.ASSOCIATION = FALSE,
                      DISP.ASSOCIATION.LARGE = FALSE,
                      DISP.REGION = FALSE,
                      DISP.REGION.LARGE = FALSE, 
                      SYMBOLS = "circle-fill",
                      SYMBOLS.LARGE = NA,
                      SAMPLE.LABELS = NULL,
                      SAMPLE.LABELS.LARGE = NULL,
                      USE.COLORS = TRUE,
                      DISP.COLOR.REF = TRUE,
                      COLOR.LIST = NULL,
                      COLOR.LIST.LARGE = NULL,
                      BIOFEAT.USER.FILE= NULL,
                      BIOFEAT.USER.TYPE= c("GeneRegion","Annotation","Data"),
                      BIOFEAT.USER.TYPE.PLOT = NULL,
                      LIST.TRACKS = "geneENSEMBL,CGI,ChromHMM,DNAse,RegENSEMBL,SNP",
                      PATTERN.REGULATION="GM12878",
                      IMAGE.TITLE = NULL,
                      IMAGE.NAME = "coMET",
                      IMAGE.TYPE = c("pdf", "eps"),
                      IMAGE.SIZE = 3.5,
                      PRINT.IMAGE = FALSE,
                      config.file = NULL,
                      VERBOSE = FALSE) {
  
  #-------------------MAIN FUNCTION----------------------------
  #LIST.TRACK DEFAULT VALUE geneENSEMBL,CGI,ChromHMM,DNAse,RegENSEMBL,SNP
  #DEBUG STATEMENT
  if (VERBOSE) cat("START COMET  VERSION WEB\n")
  
  #-------------------GLOBAL VARIABLES BEGINS------------------
  
  #create three variables (global variables), which are lists of all other variables
  #used throughout the program to be passed to functions with values changed as
  #necessary.
  gbl.var <- NULL
  
  #DATA VARIABLES
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
  
  large.split.sample.labels <- NULL
  large.split.color.list <- NULL
  large.color.list <- NULL
  large.symbol.list <- NULL
  large.fill.list <- NULL
  large.split.format <- NULL
  large.split.association <- NULL
  large.split.region <- NULL
  ref <- NULL
  samples <- NULL
  large.samples <- NULL
  
  
  palette.size <- NULL
  cex.factor <- 0.5
  cex.factor.symbol <- 0.25
  font.size <- NULL
  line.width <- 0.5
  
  cormatrix.data <- NULL
  matrix.data <- NULL
  cormatrix.data.full <- NULL
  
  mySession <- NULL
  listtracks_gviz <- NULL
  listtracks_user <- NULL
  listtracks_ggbio <- NULL
  listtracks_trackviewer <- NULL
  split.list.tracks <- NULL
  verbose <- FALSE
  
  gbl.var <- list(mydata.data = mydata.data,
                  general.data = general.data,
                  split.mydata.file =split.mydata.file,
                  split.biofeature.data.user.file = split.biofeature.data.user.file,
                  split.biofeature.data.user.type = split.biofeature.data.user.type,
                  split.biofeature.data.user.type.plot = split.biofeature.data.user.type.plot,
                  split.format = split.format,
                  split.association =  split.association,
                  split.region = split.region,
                  large.split.format = large.split.format,
                  large.split.association =  large.split.association,
                  large.split.region = large.split.region,
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
                  split.cormatrix.file = split.cormatrix.file,
                  matrix.data = matrix.data,
                  cormatrix.data.full = cormatrix.data.full,
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
  DISP.CORMATRIXMAP = TRUE
  DISP.PVALUEPLOT = TRUE
  DISP.MYDATA.NAMES = TRUE
  DISP.CONNECTING.LINES = TRUE
  DISP.MYDATA = TRUE
  DISP.TYPE = "symbol"
  BIOFEAT.USER.TYPE.PLOT="histogram"
  TRACKS.GVIZ = NULL
  TRACKS.GGBIO = NULL
  TRACKS.TRACKVIEWER = NULL
  BIOFEAT.USER.FILE= NULL
  PALETTE.FILE = NULL
  DISP.COLOR.BAR = TRUE
  DISP.PHYS.DIST = TRUE
  DISP.LEGEND = TRUE
  DISP.MARKER.LINES = TRUE
  DISP.MULT.LAB.X = FALSE
  CONNECTING.LINES.FACTOR = 1.5
  CONNECTING.LINES.ADJ = 0.01
  CONNECTING.LINES.VERT.ADJ = -1
  CONNECTING.LINES.FLEX = 0
  FONT.FACTOR = NULL
  COLOR.LIST = "red"
  SYMBOL.FACTOR = NULL
  DATASET.GENE = "hsapiens_gene_ensembl"
  DATASET.SNP="hsapiens_snp"
  VERSION.DBSNP="snp138"
  DATASET.SNP.STOMA="hsapiens_snp_som"
  DATASET.REGULATION="hsapiens_feature_set"
  DATASET.STRU="hsapiens_structvar"
  DATASET.STRU.STOMA="hsapiens_structvar_som"
  BROWSER.SESSION="UCSC"
  
  
  #-------------------UPDATE CONFIGURATION VARIABLES---------
  config.var <- list(MYDATA.FILE = MYDATA.FILE,
                     MYDATA.FORMAT = MYDATA.FORMAT,
                     MYDATA.REF = MYDATA.REF,
                     MYDATA.LARGE.FILE = MYDATA.LARGE.FILE,
                     MYDATA.LARGE.FORMAT = MYDATA.LARGE.FORMAT,
                     BIOFEAT.USER.FILE = BIOFEAT.USER.FILE,
                     BIOFEAT.USER.TYPE = BIOFEAT.USER.TYPE,
                     BIOFEAT.USER.TYPE.PLOT = BIOFEAT.USER.TYPE.PLOT,
                     GENOME = GENOME,
                     PALETTE.FILE = PALETTE.FILE,
                     LAB.Y = LAB.Y,
                     DISP.PVAL.THRESHOLD = DISP.PVAL.THRESHOLD,
                     START = START,
                     END = END,
                     ZOOM = ZOOM,
                     DISP.ASSOCIATION = DISP.ASSOCIATION,
                     DISP.ASSOCIATION.LARGE = DISP.ASSOCIATION.LARGE,
                     DISP.REGION = DISP.REGION,
                     DISP.REGION.LARGE = DISP.REGION.LARGE,
                     DATASET.GENE = DATASET.GENE,
                     DATASET.SNP = DATASET.SNP,
                     DATASET.SNP.STOMA = DATASET.SNP.STOMA,
                     DATASET.REGULATION = DATASET.REGULATION,
                     DATASET.STRU = DATASET.STRU,
                     DATASET.STRU.STOMA = DATASET.STRU.STOMA,
                     VERSION.DBSNP = VERSION.DBSNP,
                     PATTERN.REGULATION = PATTERN.REGULATION, 
                     BROWSER.SESSION = BROWSER.SESSION,
                     TRACKS.GVIZ = TRACKS.GVIZ,
                     TRACKS.GGBIO = TRACKS.GGBIO,
                     LIST.TRACKS = LIST.TRACKS,
                     TRACKS.TRACKVIEWER = TRACKS.TRACKVIEWER,
                     SYMBOLS = SYMBOLS,
                     SYMBOLS.LARGE = SYMBOLS.LARGE,
                     PVAL.THRESHOLD = PVAL.THRESHOLD,
                     DISP.TYPE = DISP.TYPE,
                     DISP.CORMATRIXMAP = DISP.CORMATRIXMAP,
                     DISP.MYDATA = DISP.MYDATA,
                     DISP.MARKER.LINES = DISP.MARKER.LINES,
                     DISP.CONNECTING.LINES = DISP.CONNECTING.LINES,
                     DISP.MYDATA.NAMES = DISP.MYDATA.NAMES,
                     DISP.PVALUEPLOT = DISP.PVALUEPLOT,
                     CORMATRIX.COLOR.SCHEME = CORMATRIX.COLOR.SCHEME,
                     COLOR.LIST = COLOR.LIST,
                     COLOR.LIST.LARGE = COLOR.LIST.LARGE,
                     USE.COLORS = USE.COLORS,
                     DISP.COLOR.REF = DISP.COLOR.REF,
                     CORMATRIX.METHOD = CORMATRIX.METHOD,
                     CORMATRIX.FILE = CORMATRIX.FILE,
                     DISP.COLOR.BAR = DISP.COLOR.BAR,
                     DISP.PHYS.DIST = DISP.PHYS.DIST,
                     IMAGE.TITLE = IMAGE.TITLE,
                     DISP.LEGEND = DISP.LEGEND,
                     SAMPLE.LABELS = SAMPLE.LABELS,
                     SAMPLE.LABELS.LARGE = SAMPLE.LABELS.LARGE,
                     IMAGE.TYPE = IMAGE.TYPE,
                     IMAGE.SIZE = IMAGE.SIZE,
                     DISP.MULT.LAB.X = DISP.MULT.LAB.X,
                     IMAGE.NAME = IMAGE.NAME,
                     PRINT.IMAGE = PRINT.IMAGE,
                     CONNECTING.LINES.FACTOR = CONNECTING.LINES.FACTOR,
                     CONNECTING.LINES.ADJ = CONNECTING.LINES.ADJ,
                     CONNECTING.LINES.VERT.ADJ = CONNECTING.LINES.VERT.ADJ,
                     CONNECTING.LINES.FLEX = CONNECTING.LINES.FLEX,
                     FONT.FACTOR = FONT.FACTOR,
                     SYMBOL.FACTOR = SYMBOL.FACTOR,
                     VERBOSE = VERBOSE)
  
  #-------------------CONFIGURATION VARIABLES BEGINS from config file---------
  
  gbl.var$verbose <- config.var$VERBOSE
  
  if(!is.null(config.file)) {
    config.var <- read.config(config.file, config.var)
  }
  
  if(!is.null(MYDATA.FILE) ){
    config.var$MYDATA.FILE <- MYDATA.FILE
  }
  if(!is.null(MYDATA.FORMAT) & length(MYDATA.FORMAT) == 1){
    config.var$MYDATA.FORMAT <- MYDATA.FORMAT
  }
  if(!is.null(MYDATA.LARGE.FORMAT) & length(MYDATA.LARGE.FORMAT) == 1){
    config.var$MYDATA.LARGE.FORMAT <- MYDATA.LARGE.FORMAT
  }
  if(!is.null(CORMATRIX.FILE)){
    config.var$CORMATRIX.FILE <- CORMATRIX.FILE
  }
  if(!is.null(CORMATRIX.FORMAT) & length(CORMATRIX.FORMAT) == 1){
    config.var$CORMATRIX.FORMAT <- CORMATRIX.FORMAT
  }
  if(!is.null(CORMATRIX.METHOD) & length(CORMATRIX.METHOD) == 1){
    config.var$CORMATRIX.METHOD <- CORMATRIX.METHOD
  }
  
  if(!is.null(BIOFEAT.USER.FILE)){
    config.var$BIOFEAT.USER.FILE <- BIOFEAT.USER.FILE
  }
  if(!is.null(PALETTE.FILE)){
    config.var$PALETTE.FILE <- PALETTE.FILE
  }
  if(!is.null(START)){
    config.var$START <- START
  }
  if(!is.null(END)){
    config.var$END <- END
  }
  if(!is.null(GENOME)){
    config.var$GENOME <- GENOME
  }
  if(is.null(config.var$MYDATA.FILE)) {
    stop("Invalid MYDATA data file: ", config.var$MYDATA.FILE, "\n")
  }
  
  #-------------- CHECK if all parameters have values -------------------
  check.configVar(config.var)
  
  #------------- connection to database
  if(!is.null(config.var$GENOME) & !is.null(config.var$BROWSER.SESSION)){
    mySession <- browserSession(config.var$BROWSER.SESSION)
    genome(mySession) <- config.var$GENOME
    gbl.var$mySession <- mySession
  }
  
  
  #------------- READ DATA for ANNOTATION TRACKS
  if(!is.null(config.var$BIOFEAT.USER.FILE)){
    split.biofeature.data.user.file <- strsplit(config.var$BIOFEAT.USER.FILE, ",")
    split.biofeature.data.user.type <- strsplit(config.var$BIOFEAT.USER.TYPE, ",")
    
    gbl.var$split.biofeature.data.user.file <- split.biofeature.data.user.file
    gbl.var$split.biofeature.data.user.type <- split.biofeature.data.user.type
    
    if(!is.null(config.var$BIOFEAT.USER.TYPE.PLOT)){
      split.biofeature.data.user.type.plot <- strsplit(config.var$BIOFEAT.USER.TYPE.PLOT, ",")
      gbl.var$split.biofeature.data.user.type.plot <- split.biofeature.data.user.type.plot
    }
  }
  
  #------------ SELECTION ANNOTATION TRACKS in comet.web
  if(!is.null(config.var$LIST.TRACKS)){
    split.list.tracks <- strsplit(config.var$LIST.TRACKS, ",")
    
    someenv.list.tracks<-hash()
    for(i in 1:length(split.list.tracks[[1]]))
    {
      info.track <- split.list.tracks[[1]][i]
      # if (config.var$VERBOSE) cat(" long",info.track,"\n")
      someenv.list.tracks[[ info.track ]]<- info.track
    } 
    gbl.var$split.list.tracks <- someenv.list.tracks
  }
  
  
  #------------- READ DATA and UPDATE VARIABLES
  if(!is.null(config.var$CORMATRIX.FILE)){
    split.cormatrix.file <- strsplit(config.var$CORMATRIX.FILE, ",")
    comatrix.file.length <- length(split.cormatrix.file[[1]])
    gbl.var$split.cormatrix.file <- split.cormatrix.file
  } else {
    warning("No visualisation of correlation matrice")
    DISP.MYDATA <- FALSE
  }
  
  if(!is.null(config.var$MYDATA.FILE)){
    split.mydata.file <- strsplit(config.var$MYDATA.FILE, ",")
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
    #split SAMPLE.LABELS on comma, then take the top 5 or less entries and create a list
    if(!is.null(config.var$SAMPLE.LABELS)) {
      split.tmp.sample.list <- strsplit(config.var$SAMPLE.LABELS, ",")
    } else {
      split.tmp.sample.list <- rep("sample", samples)
    }
    
    if(length(split.tmp.sample.list[[1]]) < samples) {
      warning("SAMPLE.LABELS contains less labels than the number of samples. Labels must be separated by commas without spaces.\n")
    }
    
    #--- FORMAT of DATA
    if(!is.null(config.var$MYDATA.FORMAT)){
      split.format <- strsplit(config.var$MYDATA.FORMAT, ",")
      split.format.length <- length(split.format[[1]])
      gbl.var$split.format <- split.format
    }
    
    #-- Visualisation of ASSOCIATION
    if(!is.null(config.var$DISP.ASSOCIATION)){
      split.association <- strsplit(config.var$DISP.ASSOCIATION, ",")
      split.association.length <- length(split.association[[1]])
      gbl.var$split.association <- split.association
    }
    
    #---VISUALISATION of REGION
    if(!is.null(config.var$DISP.REGION)){
      split.region <- strsplit(config.var$DISP.REGION, ",")
      split.region.length <- length(split.region[[1]])
      gbl.var$split.region <- split.region
    }
    
    split.tmp.sample.list.truncated <- substr(split.tmp.sample.list[[1]], 1, 25)
  } else {
    stop("No data visualise")
  }
  
  #DEBUG STATEMENT
  #if (config.var$VERBOSE) cat("samples ", samples, "\n")
  
  
  gbl.var$split.sample.labels <- list(head(split.tmp.sample.list.truncated, samples))
  
  #---LARGE REFERENCE DATA
  if(!is.null(config.var$MYDATA.LARGE.FILE)){
    large.split.mydata.file <- strsplit(config.var$MYDATA.LARGE.FILE, ",")
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
    #split SAMPLE.LABELS on comma, then take the top 5 or less entries and create a list
    large.split.tmp.sample.list <- strsplit(config.var$SAMPLE.LABELS.LARGE, ",")
    
    if(length(large.split.tmp.sample.list[[1]]) < large.samples) {
      warning("LARGE.SAMPLE.LABELS contains less labels than the number of samples. Labels must be separated by commas without spaces.\n")
    }
    
    large.split.tmp.sample.list.truncated <- substr(large.split.tmp.sample.list[[1]], 1, 25)
    
    #DEBUG STATEMENT
    #if (config.var$VERBOSE) cat("slarge.amples ", large.samples, "\n") 
    
    gbl.var$large.split.sample.labels <- list(head(large.split.tmp.sample.list.truncated, large.samples))
    
    #--- FORMAT of LARGE DATA
    if(!is.null(config.var$MYDATA.LARGE.FORMAT)){
      #if (config.var$VERBOSE) cat("ASSO", config.var$MYDATA.LARGE.FORMAT,"\n") 
      large.split.format <- strsplit(config.var$MYDATA.LARGE.FORMAT, ",")
      large.split.format.length <- length(large.split.format[[1]])
      gbl.var$large.split.format <- large.split.format
      # if (config.var$VERBOSE) cat("L",large.split.format.length,"  value ", gbl.var$large.split.format[[1]][1],"\n") 
    }
    
    #-- Visualisation of ASSOCIATION
    if(!is.null(config.var$DISP.ASSOCIATION.LARGE)){
      large.split.association <- strsplit(config.var$DISP.ASSOCIATION.LARGE, ",")
      large.split.association.length <- length(large.split.association[[1]])
      gbl.var$large.split.association <- large.split.association
    }
    
    #---VISUALISATION of REGION
    if(!is.null(config.var$DISP.REGION.LARGE)){
      large.split.region <- strsplit(config.var$DISP.REGION.LARGE, ",")
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
  if(PRINT.IMAGE == FALSE || is.null(config.var$IMAGE.NAME)){
    gbl.var <- draw.plot.comet.web(config.var, gbl.var,newpage=TRUE)
  } else{
    printPlot.comet.web(config.var, gbl.var)
  }
  
  #END FUNCTION LIST
  
  invisible(list(config.var, gbl.var))
  #DEBUG STATEMENT
  if (config.var$VERBOSE) cat("FINISH WEB COMET \n")
}


