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
# File: Analyse.R
# Author: Tiphaine Martin
# Email: tiphaine.martin@kcl.ac.uk
# Purpose: coMET allows the display of p-values from association
#           with a correlation heatmap.
# Version : 0.1
###########################################################################

#------------------READ in CONFiguration file------------------
read.config <- function(config.file, config.var) {
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)   cat("START READ.CONFIG\n")
  
  raw.config <- read.table(config.file, header=FALSE, as.is=TRUE, fill=TRUE, 
                           sep="\t", blank.lines.skip = TRUE)
  
  for(i in 1:length(raw.config$V1)) {
    tmp <- strsplit(raw.config$V1[i], "=")
    
    if(identical("MYDATA.FILE", tmp[[1]][1])) {
      config.var$MYDATA.FILE <- tmp[[1]][2]
    }
    
    if(identical("MYDATA.FORMAT", tmp[[1]][1])) {
      config.var$MYDATA.FORMAT <- tmp[[1]][2]
    }
    
    if(identical("MYDATA.LARGE.FILE", tmp[[1]][1])) {
      config.var$MYDATA.LARGE.FILE <- tmp[[1]][2]
    }
    
    if(identical("MYDATA.LARGE.FORMAT", tmp[[1]][1])) {
      config.var$MYDATA.LARGE.FORMAT <- tmp[[1]][2]
    }
    
    if(identical("CORMATRIX.FILE", tmp[[1]][1])) {
      config.var$CORMATRIX.FILE <- tmp[[1]][2]
    }
    
    if(identical("CORMATRIX.METHOD", tmp[[1]][1])) {
      config.var$CORMATRIX.METHOD <- tmp[[1]][2]
    }
    
    if(identical("CORMATRIX.FORMAT", tmp[[1]][1])) {
      config.var$CORMATRIX.FORMAT <- tmp[[1]][2]
    }
    
    if(identical("CORMATRIX.COLOR.SCHEME", tmp[[1]][1])) {
      config.var$CORMATRIX.COLOR.SCHEME <- tmp[[1]][2]
    }
    
    if(identical("MYDATA.REF", tmp[[1]][1])) {
      config.var$MYDATA.REF <- tmp[[1]][2]
    }
    
    if(identical("START", tmp[[1]][1])) {
      config.var$START <- as.numeric(tmp[[1]][2])
    }
    
    if(identical("END", tmp[[1]][1])) {
      config.var$END <- as.numeric(tmp[[1]][2])
    }
    
    if(identical("ZOOM", tmp[[1]][1])) {
      config.var$ZOOM <- as.logical(tmp[[1]][2])
    }
    
    if(identical("BIOFEAT.USER.FILE", tmp[[1]][1])) {
      config.var$BIOFEAT.USER.FILE <- tmp[[1]][2]
    }
    
    if(identical("BIOFEAT.USER.TYPE", tmp[[1]][1])) {
      config.var$BIOFEAT.USER.TYPE <- tmp[[1]][2]
    }
    
    if(identical("BIOFEAT.USER.TYPE.PLOT", tmp[[1]][1])) {
      config.var$BIOFEAT.USER.TYPE.PLOT <- tmp[[1]][2]
    }
    
    if(identical("LAB.Y", tmp[[1]][1])) {
      config.var$LAB.Y <- tmp[[1]][2]
    }
    
    if(identical("SYMBOLS", tmp[[1]][1])) {
      config.var$SYMBOLS <- tmp[[1]][2]
    }
    
    if(identical("SYMBOLS.LARGE", tmp[[1]][1])) {
      config.var$SYMBOLS.LARGE <- tmp[[1]][2]
    }
    
    if(identical("SAMPLE.LABELS", tmp[[1]][1])) {
      config.var$SAMPLE.LABELS <- tmp[[1]][2]
    }
    
    if(identical("SAMPLE.LABELS.LARGE", tmp[[1]][1])) {
      config.var$SAMPLE.LABELS.LARGE <- tmp[[1]][2]
    }
    
    if(identical("PVAL.THRESHOLD", tmp[[1]][1])) {
      config.var$PVAL.THRESHOLD <- as.numeric(tmp[[1]][2])
    }
    
    if(identical("DISP.PVAL.THRESHOLD", tmp[[1]][1])) {
      config.var$DISP.PVAL.THRESHOLD <- as.numeric(tmp[[1]][2])
    }
    
    if(identical("DISP.COLOR.REF", tmp[[1]][1])) {
      config.var$DISP.COLOR.REF <- as.logical(tmp[[1]][2])
    }
    
    if(identical("USE.COLORS", tmp[[1]][1])) {
      config.var$USE.COLORS <- as.logical(tmp[[1]][2])
    }
    
    if(identical("COLOR.LIST", tmp[[1]][1])) {
      config.var$COLOR.LIST <- tmp[[1]][2]
    }
    
    if(identical("COLOR.LIST.LARGE", tmp[[1]][1])) {
      config.var$COLOR.LIST.LARGE <- tmp[[1]][2]
    }
    if(identical("DISP.ASSOCIATION", tmp[[1]][1])) {
      config.var$DISP.ASSOCIATION <- tmp[[1]][2]
    }
    
    if(identical("DISP.ASSOCIATION.LARGE", tmp[[1]][1])) {
      config.var$DISP.ASSOCIATION.LARGE <- tmp[[1]][2]
    }
    
    if(identical("DISP.REGION", tmp[[1]][1])) {
      config.var$DISP.REGION <- tmp[[1]][2]
    }
    
    if(identical("DISP.REGION.LARGE", tmp[[1]][1])) {
      config.var$DISP.REGION.LARGE <- tmp[[1]][2]
    }
    
    if(identical("DISP.MYDATA.NAMES", tmp[[1]][1])) {
      config.var$DISP.MYDATA.NAMES <- as.logical(tmp[[1]][2])
    }
    
    if(identical("DISP.PVALUEPLOT", tmp[[1]][1])) {
      config.var$DISP.PVALUEPLOT <- tmp[[1]][2]
    }
    
    if(identical("LIST.TRACKS", tmp[[1]][1])) {
      config.var$LIST.TRACKS <- tmp[[1]][2]
    }
    
    if(identical("GENOME", tmp[[1]][1])) {
      config.var$GENOME <- tmp[[1]][2]
    }
    
    if(identical("DATASET.GENE", tmp[[1]][1])) {
      config.var$DATASET.GENE <- tmp[[1]][2]
    }
    
    if(identical("DATASET.SNP", tmp[[1]][1])) {
      config.var$DATASET.SNP <- tmp[[1]][2]
    }
    
    if(identical("VERSION.DBSNP", tmp[[1]][1])) {
      config.var$VERSION.DBSNP <- tmp[[1]][2]
    }
    
    if(identical("DATASET.SNP.STOMA", tmp[[1]][1])) {
      config.var$DATASET.SNP.STOMA <- tmp[[1]][2]
    }
    
    if(identical("DATASET.REGULATION", tmp[[1]][1])) {
      config.var$DATASET.REGULATION <- tmp[[1]][2]
    }
    
    if(identical("DATASET.STRU", tmp[[1]][1])) {
      config.var$DATASET.STRU <- tmp[[1]][2]
    }
    
    if(identical("DATASET.STRU.STOMA", tmp[[1]][1])) {
      config.var$DATASET.STRU.STOMA <- tmp[[1]][2]
    }
    
    if(identical("PATTERN.REGULATION", tmp[[1]][1])) {
      config.var$PATTERN.REGULATION <- tmp[[1]][2]
    }
    
    if(identical("BROWSER.SESSION", tmp[[1]][1])) {
      config.var$BROWSER.SESSION <- tmp[[1]][2]
    }
    
    if(identical("TRACKS.GVIZ", tmp[[1]][1])) {
      config.var$TRACKS.GVIZ <- tmp[[1]][2]
    }
    
    if(identical("TRACKS.GGBIO", tmp[[1]][1])) {
      config.var$TRACKS.GGBIO <- tmp[[1]][2]
    }
    
    if(identical("TRACKS.TRACKVIEWER", tmp[[1]][1])) {
      config.var$TRACKS.TRACKVIEWER <- tmp[[1]][2]
    }
    
    if(identical("PALETTE.FILE", tmp[[1]][1])) {
      config.var$PALETTE.FILE <- tmp[[1]][2]
    }
    
    if(identical("DISP.TYPE", tmp[[1]][1])) {
      config.var$DISP.TYPE <- tmp[[1]][2]
    }
    
    if(identical("DISP.CORMATRIXMAP", tmp[[1]][1])) {
      config.var$DISP.CORMATRIXMAP <- as.logical(tmp[[1]][2])
    }
    
    if(identical("DISP.MYDATA", tmp[[1]][1])) {
      config.var$DISP.MYDATA <- as.logical(tmp[[1]][2])
    }
    
    if(identical("DISP.COLOR.BAR", tmp[[1]][1])) {
      config.var$DISP.COLOR.BAR <- as.logical(tmp[[1]][2])
    }
    
    if(identical("DISP.PHYS.DIST", tmp[[1]][1])) {
      config.var$DISP.PHYS.DIST <- as.logical(tmp[[1]][2])
    }
    
    if(identical("DISP.CONNECTING.LINES", tmp[[1]][1])) {
      config.var$DISP.CONNECTING.LINES <- as.logical(tmp[[1]][2])
    }
    
    if(identical("IMAGE.TITLE", tmp[[1]][1])) {
      config.var$IMAGE.TITLE <- tmp[[1]][2]
    }
    
    if(identical("DISP.LEGEND", tmp[[1]][1])) {
      config.var$DISP.LEGEND <- as.logical(tmp[[1]][2])
    }
    
    if(identical("IMAGE.TYPE", tmp[[1]][1])) {
      config.var$IMAGE.TYPE <- tmp[[1]][2]
    }
    
    if(identical("IMAGE.SIZE", tmp[[1]][1])) {
      config.var$IMAGE.SIZE <- tmp[[1]][2]
    }
    
    if(identical("DISP.MULT.LAB.X", tmp[[1]][1])) {
      config.var$DISP.MULT.LAB.X <- as.logical(tmp[[1]][2])
    }
    
    if(identical("PRINT.IMAGE", tmp[[1]][1])) {
      config.var$PRINT.IMAGE <- as.logical(tmp[[1]][2])
    }
    
    if(identical("IMAGE.NAME", tmp[[1]][1])) {
      config.var$IMAGE.NAME <- tmp[[1]][2]
    }
    
    if(identical("CONNECTING.LINES.FACTOR", tmp[[1]][1])) {
      config.var$CONNECTING.LINES.FACTOR <- as.numeric(tmp[[1]][2])
    }
    
    if(identical("CONNECTING.LINES.ADJ", tmp[[1]][1])) {
      config.var$CONNECTING.LINES.ADJ <- as.numeric(tmp[[1]][2])
    }
    
    if(identical("CONNECTING.LINES.VERT.ADJ", tmp[[1]][1])) {
      config.var$CONNECTING.LINES.VERT.ADJ <- as.numeric(tmp[[1]][2])
    }
    
    if(identical("CONNECTING.LINES.FLEX", tmp[[1]][1])) {
      config.var$CONNECTING.LINES.FLEX <- as.numeric(tmp[[1]][2])
    }
    
    if(identical("FONT.FACTOR", tmp[[1]][1])) {
      config.var$FONT.FACTOR <- as.numeric(tmp[[1]][2])
    }
    
    if(identical("SYMBOL.FACTOR", tmp[[1]][1])) {
      config.var$SYMBOL.FACTOR <- as.numeric(tmp[[1]][2])
    }
    
    if(identical("VERBOSE", tmp[[1]][1])) {
      config.var$VERBOSE <- as.numeric(tmp[[1]][2])
    }
  }
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("FINISH READ.CONFIG\n")
  
  return(config.var)
}

#------------------Check the list of parameters in CONFiguration file------------------
check.configVar <- function(config.var) {
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("START CHECK.CONFIG\n")
  
  #######################INFO FILE 
  if(is.null(config.var$MYDATA.FILE)) {
    stop("No value for info file.\n 
          Parametre: MYDATA.FILE")
  }
  
  split.mydata.file <- strsplit(config.var$MYDATA.FILE, ",")
  if( length(split.mydata.file[[1]]) >1 ) {
    stop("Need only one file. Put other file via MYDATA.FILE \n")
  }
  
  if(is.null(config.var$MYDATA.FORMAT)) {
    stop("No value for the format of info file.\n 
          Parametre: MYDATA.FORMAT.\n 
          Potential value: SITE, REGION, SITE_ASSO, REGION_ASSO\n")
  }
  
  if(is.null(config.var$SYMBOLS)) {
    stop("No value for the format of info file.\n 
          Parametre: SYMBOLS\n 
          Potential value: circle-fill, square-fill diamont-fill\n")
  }
  
  if(is.null(config.var$SAMPLE.LABELS)) {
    stop("No value for the label of info file.\n 
          Parametre: SAMPLE.LABELS\n")
  }
  
  if(is.null(config.var$COLOR.LIST)) {
    if (config.var$VERBOSE)   cat("No value for the color of info file.\n 
          Parametre: COLOR.LIST\n")
  }
  
  
  if(grepl("ASSO", config.var$MYDATA.FORMAT)[1] & is.null(config.var$DISP.ASSOCIATION)) {
    stop("No value for the association of info file.\n 
          Parametre: DISP.ASSOCIATION\n")
  }
  
  if(grepl("REGION", config.var$MYDATA.FORMAT)[1] & is.null(config.var$DISP.REGION)) {
    stop("No value for the region of info file.\n 
          Parametre: DISP.REGION\n")
  }
  
  #---------------- Supplementary file
  if(is.null(config.var$MYDATA.LARGE.FILE)) {
    if (config.var$VERBOSE) cat("No value for the supplementary file.\n 
          Parametre: MYDATA.LARGE.FILE ")
  }
  
  if( ! is.null(config.var$MYDATA.LARGE.FILE)) {
    if(is.null(config.var$MYDATA.LARGE.FORMAT)) {
      stop("No value for the format of supplementary data.\n 
          Parametre: MYDATA.LARGE.FORMAT.\n 
          Potential value: SITE, REGION, SITE_ASSO, REGION_ASSO\n")
    }
    
    if(is.null(config.var$SYMBOLS.LARGE)) {
      stop("No value of symbole for supplementary data\n
           Parameter: SYMBOLS.LARGE\n")
    }
    
    if(is.null( config.var$SAMPLE.LABELS.LARGE)) {
      stop("No value of label for supplementary data\n
           Parameter: SAMPLE.LABELS.LARGE\n")
    }
    
    if(is.null(config.var$COLOR.LIST.LARGE)) {
      stop("No value of color for supplementary data\n
           Parameter: COLOR.LIST.LARGE\n")
    }
    
    if(grepl("ASSO", config.var$MYDATA.LARGE.FORMAT)[1] & 
         is.null(config.var$DISP.ASSOCIATION.LARGE)) {
      stop("No value of association for supplementary data\n
           Parameter: DISP.ASSOCIATION.LARGE\n")
    }
    
    if(grepl("REGION", config.var$MYDATA.LARGE.FORMAT)[1] &
         is.null(config.var$DISP.REGION.LARGE)) {
      stop("No value of region for supplementary data\n
           Parameter: DISP.REGION.LARGE\n")
    }
    
    split.mydata.large.file <- strsplit(config.var$MYDATA.LARGE.FILE, ",")
    split.mydata.large.format <- strsplit(config.var$MYDATA.LARGE.FORMAT, ",")
    split.mydata.large.region <- strsplit(config.var$DISP.REGION.LARGE, ",")
    split.mydata.large.asso <- strsplit(config.var$DISP.ASSOCIATION.LARGE, ",")
    split.mydata.large.color <- strsplit(config.var$COLOR.LIST.LARGE, ",")
    split.mydata.large.label <- strsplit(config.var$SAMPLE.LABELS.LARGE, ",")
    split.mydata.large.symbol <- strsplit(config.var$SYMBOLS.LARGE, ",")
    
    if( length(split.mydata.large.file[[1]]) !=  length(split.mydata.large.format[[1]]) ||
          length(split.mydata.large.file[[1]]) !=  length(split.mydata.large.region[[1]]) ||
          length(split.mydata.large.file[[1]]) !=  length(split.mydata.large.asso[[1]]) ||
          length(split.mydata.large.file[[1]]) !=  length(split.mydata.large.color[[1]]) ||
          length(split.mydata.large.file[[1]]) !=  length(split.mydata.large.label[[1]]) ||
          length(split.mydata.large.file[[1]]) !=  length(split.mydata.large.symbol[[1]])) {
      stop("Need to have the same number of element for all options related to supplementary data\n
           MYDATA.LARGE.FILE, MYDATA.LARGE.FORMAT, DISP.REGION.LARGE, DISP.ASSOCIATION.LARGE, COLOR.LIST.LARGE,
           SAMPLE.LABELS.LARGE, SYMBOLS.LARGE\n")
    }
  } 
  
  #-----------------------Correlation matrix
  if(is.null(config.var$CORMATRIX.FILE)) {
    stop("No value for the correlation file.\n 
          Parametre: CORMATRIX.FILE ")
  } else {
    if(is.null(config.var$DISP.MYDATA)) {
      stop("No value to show the correlation file.\n 
          Parametre: DISP.MYDATA ")
    }
  }
  
  if(is.null(config.var$CORMATRIX.FORMAT)) {
    stop("No value for the correlation file.\n 
          Parametre: CORMATRIX.METHOD \n
          Potential value: CORMATRIX, RAW")
  }
  
  if(config.var$CORMATRIX.FORMAT == "RAW" & 
       is.null(config.var$CORMATRIX.METHOD)) {
    stop("No value for the method of correlation file for raw.\n 
          Parametre: CORMATRIX.METHOD \n
          Potential value: spearman, pearson, kendall")
  }
  
  if(is.null(config.var$CORMATRIX.COLOR.SCHEME)) {
    stop("No value for the method of correlation file for raw.\n 
          Parametre: CORMATRIX.COLOR.SCHEME \n
          Potential value: spearman, pearson, kendall")
  }
  
  if(is.null(config.var$DISP.CORMATRIXMAP)) {
    stop("Show heatmap of correlation matrix.\n 
          Parametre: DISP.CORMATRIXMAP \n
          Potential value: TRUE or FALSE")
  }
  
  #-----------------------Parameters
  if(is.null(config.var$MYDATA.REF)) {
    if (config.var$VERBOSE)   cat("No value for the reference.\n
              Parameter: MYDATA.REF\n")
  }
  
  if(is.null(config.var$START)) {
    if (config.var$VERBOSE)   cat("No value for the beginning of genomic region\n
              Parameter: START\n")
  }
  
  if(is.null(config.var$END)) {
    if (config.var$VERBOSE)   cat("No value for the end of genomic region\n
            Parameter: END\n")
  }
  
  if(is.null(config.var$LAB.Y)) {
    stop("No value of axis\n
           Parameter: LAB.Y\n
           Potential value: log10, ln \n")
  }
  
  if(is.null(config.var$PVAL.THRESHOLD)) {
    stop("No value of threshold of the significance\n
           Parameter: PVAL.THRESHOLD\n
           Potential value: 10e-8\n")
  }
  
  if(is.null(config.var$DISP.PVAL.THRESHOLD)) {
    stop("No value of threshold to visualise the pvalue\n
           Parameter: DISP.PVAL.THRESHOLD\n
           Potential value: 0\n")
  }
  
  if(is.null(config.var$DISP.COLOR.REF)) {
    stop("No value of threshold to visualise the pvalue\n
           Parameter: DISP.COLOR.REF\n
           Potential value: TRUE or FALSE\n")
  }
  
  if(is.null(config.var$USE.COLORS)) {
    stop("No value to use\n
           Parameter: USE.COLORS\n
           Potential value: TRUE or FALSE\n")
  }
  
  if(is.null(config.var$DISP.MYDATA.NAMES)) {
    stop("No value to show the name of genes\n
           Parameter: DISP.MYDATA.NAMES\n
           Potential value: TRUE or FALSE\n")
  }
  
  if(is.null(config.var$DISP.PVALUEPLOT)) {
    stop("No value to show the pvalue plot \n
           Parameter: DISP.PVALUEPLOT\n
           Potential value: TRUE or FALSE\n")
  }
  
  if(is.null(config.var$GENOME)) {
    stop("No value of genome\n
           Parameter: GENOME\n
           Potential value: hg19\n")
  }
  
  
  if(is.null(config.var$BROWSER.SESSION)) {
    stop("No value of name of browser session\n
           Parameter: BROWSER.SESSION\n
           Potential value: UCSC\n")
  }
  
  if(is.null(config.var$LIST.TRACKS)) {
    if (config.var$VERBOSE)   cat("No value of annotation tracks\n
           Parameter: LIST.TRACKS\n
           Potential value: geneENSEMBL,CGI,ChromHMM,DNAse,RegENSEMBL,SNP\n")
  }else {
    
    if(is.null(config.var$TRACKS.GVIZ)) {
      if (config.var$VERBOSE)   cat("No list of Gviz track\n
           Parameter: TRACKS.GVIZ\n")
    }
    
    if(is.null(config.var$TRACKS.GGBIO)) {
      if (config.var$VERBOSE)   cat("No list of GGbio track\n
           Parameter: TRACKS.GGBIO\n")
    }
    
    if(is.null(config.var$TRACKS.TRACKVIEWER)) {
      if (config.var$VERBOSE)   cat("No list of TrackViewer track\n
           Parameter: TRACKS.TRACKVIEWER\n")
    }
  }
  
  if(is.null(config.var$DATASET.GENE)) {
    if (config.var$VERBOSE)   cat("No value of dataset for ENSEMBL\n
           Parameter: DATASET.GENE\n
          Potential value: hsapiens_gene\n")
  }
  
  if(is.null(config.var$DATASET.SNP)) {
    if (config.var$VERBOSE)   cat("No value of dataset for ENSEMBL\n
          Parameter: DATASET.SNP\n
          Potential value: hsapiens_snp\n")
  }else {
    if(is.null(config.var$VERSION.DBSNP)) {
      if (config.var$VERBOSE)   cat("No value of version of DBsnp\n
            Parameter: VERSION.DBSNP\n
            Potential value: snp138\n")
    }
  }
  
  if(is.null(config.var$DATASET.SNP.STOMA)) {
    if (config.var$VERBOSE)   cat("No value of dataset for ENSEMBL\n
          Parameter: DATASET.SNP.STOMA\n
          Potential value: hsapiens_snp_som\n")
  }
  
  if(is.null(config.var$DATASET.REGULATION)) {
    if (config.var$VERBOSE)   cat("No value of dataset for ENSEMBL\n
          Parameter: DATASET.REGULATION\n
          Potential value: hsapiens_feature_set\n")
  }
  
  if(is.null(config.var$DATASET.STRU)) {
    if (config.var$VERBOSE)   cat("No value of dataset for ENSEMBL\n
          Parameter: DATASET.STRU\n
          Potential value: hsapiens_structvar\n")
  }
  
  if(is.null(config.var$DATASET.STRU.STOMA)) {
    if (config.var$VERBOSE)   cat("No value of dataset for ENSEMBL\n
          Parameter: DATASET.STRU\n
          Potential value: hsapiens_structvar_som\n")
  }
  
  if(is.null(config.var$PATTERN.REGULATION)) {
    if (config.var$VERBOSE)   cat("No value of name of tissue or pattern of DBregulation for ENSEMBL\n
          Parameter: PATTERN.REGULATION\n
          Potential value: GM12878\n")
  }
  
  if(is.null(config.var$PALETTE.FILE)) {
    if (config.var$VERBOSE)   cat("No list of TrackViewer track\n
           Parameter: PALETTE.FILE\n")
  }
  
  if(is.null(config.var$BIOFEAT.USER.FILE)) {
    if (config.var$VERBOSE)   cat("No user file to annotation\n
           Parameter: BIOFEAT.USER.FILE\n")
  } else {
    if(is.null(config.var$BIOFEAT.USER.TYPE)) {
      stop("No type of visualise of user file\n
           Parameter: BIOFEAT.USER.TYPE\n")
    }
    
    if(is.null(config.var$BIOFEAT.USER.TYPE.PLOT)) {
      stop("No type of visualise of user file\n
           Parameter: BIOFEAT.USER.TYPE.PLOT\n")
    }
    
    split.biouser.file <- strsplit(config.var$BIOFEAT.USER.FILE, ",")
    split.biouser.type <- strsplit(config.var$BIOFEAT.USER.TYPE, ",")
    split.biouser.plot <- strsplit(config.var$BIOFEAT.USER.TYPE.PLOT, ",")
    
    
    if( length(split.biouser.file[[1]]) !=  length(split.biouser.type[[1]]) ||
          length(split.biouser.file[[1]]) !=  length(split.biouser.plot[[1]])  ) {
      stop("Need to have the same number of element for all options related to supplementary data\n
           BIOFEAT.USER.TYPE.PLOT, BIOFEAT.USER.FILE, BIOFEAT.USER.TYPE, BIOFEAT.USER.TYPE.PLOT\n")
    }
  }
  
  if(is.null(config.var$DISP.COLOR.BAR)) {
    stop("No value to show the color bar\n
           Parameter: DISP.COLOR.BAR\n
       Potential value; TRUE or FALSE\n")
  }
  
  if(is.null(config.var$DISP.PHYS.DIST)) {
    stop("No value to show the physique distance\n
           Parameter: DISP.PHYS.DIST\n
       Potential value; TRUE or FALSE\n")
  }
  
  if(is.null(config.var$DISP.CONNECTING.LINES)) {
    stop("No value to show the connecting lines\n
           Parameter: DISP.CONNECTING.LINES\n
       Potential value; TRUE or FALSE\n")
  }else {
    if(is.null(config.var$CONNECTING.LINES.FACTOR)) {
      stop("No value for the connecting lines factor \n
           Parameter: CONNECTING.LINES.FACTOR\n")
    }
    
    if(is.null(config.var$CONNECTING.LINES.ADJ)) {
      stop("No value for the connecting lines adj \n
           Parameter: CONNECTING.LINES.ADJ\n")
    }
    
    if(is.null(config.var$CONNECTING.LINES.VERT.ADJ)) {
      stop("No value for the connecting lines vertical adj \n
           Parameter: CONNECTING.LINES.VERT.ADJ\n")
    }
    
    if(is.null(config.var$CONNECTING.LINES.FLEX)) {
      stop("No value for the connecting lines flexible \n
           Parameter: CONNECTING.LINES.FLEX\n")
    }
  }
  
  if(is.null(config.var$DISP.LEGEND)) {
    stop("No value to show the legend\n
           Parameter: DISP.LEGEND\n
       Potential value; TRUE or FALSE\n")
  }
  
  if(is.null(config.var$DISP.MULT.LAB.X)) {
    if (config.var$VERBOSE)   cat("No value to show the multiple lable\n
           Parameter: DISP.MULT.LAB.X\n
       Potential value; TRUE or FALSE\n")
  }
  
  if(is.null(config.var$FONT.FACTOR)) {
    if (config.var$VERBOSE)   cat("No value for the font factor\n
           Parameter: FONT.FACTOR\n")
  }
  
  if(is.null(config.var$SYMBOL.FACTOR)) {
    if (config.var$VERBOSE)   cat("No value for the symbol factor\n
           Parameter: SYMBOL.FACTOR\n")
  }
  
  if(is.null(config.var$IMAGE.TITLE)) {
    stop("No title for the plot\n
           Parameter: IMAGE.TITLE\n")
  }
  
  if(is.null(config.var$IMAGE.TYPE)) {
    stop("No value for the type of image\n
           Parameter: IMAGE.TYPE\n
       Potential value; PDF or EPS\n")
  }
  
  if(is.null(config.var$IMAGE.SIZE)) {
    stop("No value for the size of image\n
           Parameter: IMAGE.SIZE\n
       Potential value; 3.5 or 7\n")
  }
  
  if(is.null(config.var$PRINT.IMAGE)) {
    stop("No value for the size of image\n
           Parameter: IMAGE.SIZE\n
       Potential value; 3.5 or 7\n")
  }
  
  if(is.null(config.var$IMAGE.NAME)) {
    stop("No value for the name of image\n
           Parameter: IMAGE.NAME\n")
  }
  
  #  if(is.null(config.var$ZOOM)) {
  #    if (config.var$VERBOSE)   cat("No value for the end of genomic region\n")
  #  }
  #  if(is.null(config.var$DISP.TYPE)) {
  #    if (config.var$VERBOSE)   cat("No value for the end of genomic region\n")
  #  }
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("FINISH CHECK.CONFIG\n")
  
  return(config.var)
}

#------------------UPDATE VARIABLE------------------
retrieve.data <- function(config.var, gbl.var) {
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("START RETRIEVE.DATA.MYDATA\n")
  
  #initialize distance variables
  min.dist <- NULL
  max.dist <- NULL
  
  gbl.var$mydata.gen <- config.var$GENOME
  split.mydata.file <- gbl.var$split.mydata.file
  split.cormatrix.file <- NULL
  split.cormatrix.file <- gbl.var$split.cormatrix.file
  #-------------------- READ FILES MYDATA
  for(i in 1:length(split.mydata.file[[1]])) {
    gbl.var$cur.sample <- i
    #-------------------- READ FILE CpG DATA
    gbl.var <-read.file.mydata(split.mydata.file,config.var, gbl.var,i)
    if (config.var$VERBOSE)  cat("READ ",i," MYDATA file \n")
    
    cur.mydata.num <- nrow(gbl.var$mydata.data)
    gbl.var$palette.size <- cur.mydata.num
    #create 2 hash tables containing MYDATA names and positions for current and previous sample sets
    #one hash has the name as the key and the other the location because R has no way of looking
    #up hash keys based on values
    format <- as.character(gbl.var$split.format[[1]][gbl.var$cur.sample])
    if (config.var$VERBOSE)  cat("format",format)
    if(grepl("REGION", format)[1] ) {
      mydata.hash.names.start <- new.env(hash=TRUE)
      mydata.hash.names.end <- new.env(hash=TRUE)
      mydata.hash.names.pos <- new.env(hash=TRUE)
      mydata.hash.pos.names <- new.env(hash=TRUE)
      mydata.hash.list.names.pos <- new.env(hash=TRUE)
      mydata.hash.list.pos.names <- new.env(hash=TRUE)
      
      for(p in 1:cur.mydata.num) {
        if(!is.na(gbl.var$mydata.data$CHROMOSOME[p]) & is.null(gbl.var$chr) ){
          gbl.var$mydata.chr <- gbl.var$mydata.data$CHROMOSOME[p] 
        }
        localiation <- 0
        if(!is.null(gbl.var$mydata.data$LOC.START[p]) & !is.null(gbl.var$mydata.data$LOC.END[p]) &  ( gbl.var$mydata.data$LOC.START[p] < gbl.var$mydata.data$LOC.END[p]  | gbl.var$mydata.data$LOC.START[p] == gbl.var$mydata.data$LOC.END[p]) ){
          localiation <- gbl.var$mydata.data$LOC.START[p] + ((gbl.var$mydata.data$LOC.END[p] -gbl.var$mydata.data$LOC.START[p]))/2
        }else {
          stop("Problem in the coordinates of region\n")
        }
        
        if(!is.na(gbl.var$mydata.data$MYDATA.NAME[p]) ) {
          assign(gbl.var$mydata.data$MYDATA.NAME[p], list(position=localiation), envir=mydata.hash.list.names.pos)
          assign(as.character(localiation), list(mydata.name=gbl.var$mydata.data$MYDATA.NAME[p]), envir=mydata.hash.list.pos.names)
          mydata.hash.names.pos[[gbl.var$mydata.data$MYDATA.NAME[p]]] <- localiation
          mydata.hash.pos.names[[as.character(localiation)]]<-gbl.var$mydata.data$MYDATA.NAME[p]
          mydata.hash.names.start[[gbl.var$mydata.data$MYDATA.NAME[p]]] <- gbl.var$mydata.data$LOC.START[p] 
          mydata.hash.names.end[[ gbl.var$mydata.data$MYDATA.NAME[p] ]] <- gbl.var$mydata.data$LOC.END[p]
          
        }
      }
      gbl.var$mydata.hash.list.names.pos <- mydata.hash.list.names.pos
      gbl.var$mydata.hash.names.pos <- mydata.hash.names.pos
      gbl.var$mydata.hash.names.start <- mydata.hash.names.start
      gbl.var$mydata.hash.names.end <- mydata.hash.names.end
      gbl.var$mydata.hash.list.pos.names <- mydata.hash.list.pos.names
      gbl.var$mydata.hash.pos.names <- mydata.hash.pos.names
      
    } else {
      mydata.hash.list.names.pos <- new.env(hash=TRUE)
      mydata.hash.list.pos.names <- new.env(hash=TRUE)
      mydata.hash.names.pos <- new.env(hash=TRUE)
      mydata.hash.pos.names <- new.env(hash=TRUE)
      for(p in 1:cur.mydata.num) {
        if(!is.na(gbl.var$mydata.data$CHROMOSOME[p]) & is.null(gbl.var$chr) ){
          gbl.var$mydata.chr <- gbl.var$mydata.data$CHROMOSOME[p] 
        }
        
        if(!is.na(gbl.var$mydata.data$MYDATA.NAME[p]) | !is.na(gbl.var$mydata.data$LOC[p])) {
          assign(gbl.var$mydata.data$MYDATA.NAME[p], list(position=gbl.var$mydata.data$LOC[p]), envir=mydata.hash.list.names.pos)
          assign(as.character(gbl.var$mydata.data$LOC[p]), list(mydata.name=gbl.var$mydata.data$MYDATA.NAME[p]), envir=mydata.hash.list.pos.names)
          mydata.hash.names.pos[[ gbl.var$mydata.data$MYDATA.NAME[p] ]] <- gbl.var$mydata.data$LOC[p]
          mydata.hash.pos.names[[as.character(gbl.var$mydata.data$LOC[p]) ]] <- gbl.var$mydata.data$MYDATA.NAME[p]
        }
      }
      gbl.var$mydata.hash.list.names.pos <- mydata.hash.list.names.pos
      gbl.var$mydata.hash.list.pos.names <- mydata.hash.list.pos.names
      gbl.var$mydata.hash.names.pos <- mydata.hash.names.pos
      gbl.var$mydata.hash.pos.names <- mydata.hash.pos.names
    }
    
    #create a variable with only the MYDATA names
    mydata.hash.names <- ls(envir=gbl.var$mydata.hash.list.names.pos)
    #create a variable with only the MYDATA positions
    mydata.hash.pos <- ls(envir=gbl.var$mydata.hash.list.pos.names)
    
    
    #DEBUG STATEMENT
    ## if (config.var$VERBOSE)  cat("length(mydata.hash.names) ", length(mydata.hash.names), "\n")
    ## if (config.var$VERBOSE)  cat("length(mydata.hash.pos) ", length(mydata.hash.pos), "\n")
    
    #create a variable with only the sorted MYDATA positions
    gbl.var$sorted.mydata.pos <- sort(as.numeric(mydata.hash.pos))
    
    #sort the MYDATA names
    for(n in 1:length(mydata.hash.pos)) {
      #  if (config.var$VERBOSE)  cat("POS: ", mydata.hash.pos[i], "NAME: ", get(mydata.hash.pos[i], envir=mydata.hash.pos.name)$mydata.name, "\n")
      gbl.var$sorted.mydata.names[n] <- get(mydata.hash.pos[n], envir=gbl.var$mydata.hash.list.pos.names)$mydata.name
    }
    
    #get the number of unique MYDATA in all current sample sets
    gbl.var$mydata.num <- length(gbl.var$sorted.mydata.names)
    
    #DEBUG STATEMENT
     if (config.var$VERBOSE)  cat("gbl.var$mydata.num ", gbl.var$mydata.num, "\n")
    
    #------------- Define Format visualisation CpG data
    
    #Set the gbl.var$cex.factor, which determines font size
    formatting.var.list <- set.image.parameters(config.var, gbl.var)
    
    gbl.var$font.size <- formatting.var.list$font.size
    gbl.var$line.width <- formatting.var.list$line.width
    gbl.var$cex.factor <- formatting.var.list$cex.factor
    gbl.var$cex.factor.symbol <- formatting.var.list$cex.factor.symbol
    
    #DEBUG STATEMENT
    # if (config.var$VERBOSE)  cat("gbl.var$cex.factor ", gbl.var$cex.factor, "\n")
    
    if(config.var$IMAGE.SIZE == 3.5) {
      gbl.var$mydata.names <- substr(gbl.var$mydata.data$MYDATA.NAME, 1, 9)
    } else if(config.var$IMAGE.SIZE == 7) {
      gbl.var$mydata.names <- substr(gbl.var$mydata.data$MYDATA.NAME, 1, 9)
    } else {
      stop("Invalid image size: ", config.var$IMAGE.SIZE, "\n")
    }
    
    #----------------- DEFINE VALUES OF X,Y and DISTANCE to position data
    if(!is.null(config.var$START)) {
        min.dist <- config.var$START
    } else {
      if(is.null(min.dist)){
        if(!is.null(gbl.var$mydata.data$LOC)) {
          stop("No value for LOC (min.dist)",gbl.var$mydata.data$LOC)
        }
        min.dist <- min(gbl.var$mydata.data$LOC, na.rm=TRUE)
      }else {
        min.dist <- min(c(min.dist, gbl.var$mydata.data$LOC), na.rm = TRUE)
      }
      min.dist <- min.dist - 500
    }
    
    
    gbl.var$min.user.x <- min.dist
    
    if(!is.null(config.var$END)){
      max.dist <- config.var$END 
      #max.dist <- config.var$END
    } else {
      
    if(is.null(max.dist)){
      if(is.null(gbl.var$mydata.data$LOC)){
        stop("NO value for LOC (max.dist)", gbl.var$mydata.data$LOC)
      }
        max.dist <-   max(gbl.var$mydata.data$LOC, na.rm = TRUE)
    }else {
      max.dist <- max(c(max.dist, gbl.var$mydata.data$LOC), na.rm = TRUE)
    }
      max.dist <- max.dist + 500
    }
    
    gbl.var$max.user.x <- max.dist
    
    gbl.var$total.dist <- max.dist - min.dist
    gbl.var$min.dist <- min.dist
    gbl.var$max.dist <- max.dist
    
    
    # FIX VALUES and change scale of Pvalue
    fix.var <- fix.values(config.var, gbl.var)
    config.var <- fix.var$config.var
    gbl.var <- fix.var$gbl.var
    
    if(!is.null(config.var$START)) {
      gbl.var$min.x <- config.var$START 
    } else {
      if(is.null(gbl.var$mydata.data$LOC)){
        stop("NULL value for min.x", gbl.var$mydata.data$LOC)
      }
      if(is.null(gbl.var$min.x)){
        gbl.var$min.x <-  min(gbl.var$mydata.data$LOC, na.rm = TRUE)
      }else {
        gbl.var$min.x <- min(c(gbl.var$min.x, gbl.var$mydata.data$LOC), na.rm = TRUE)
      }
      gbl.var$min.x <- gbl.var$min.x - 500 
    }
    #  if (config.var$VERBOSE)  cat("gbl.var$min.x", gbl.var$min.x, "\n")
    
    
    if(!is.null(config.var$END)) {
      gbl.var$max.x <- config.var$END 
    } else {
      if(is.null(gbl.var$mydata.data$LOC)){
        stop("NULL value for max.x", gbl.var$mydata.data$LOC)
      }
      if(is.null(gbl.var$max.x)){
      gbl.var$max.x <- max(gbl.var$mydata.data$LOC, na.rm = TRUE)
    }else {
      gbl.var$max.x <- max(c(gbl.var$max.x, gbl.var$mydata.data$LOC), na.rm = TRUE)
    }
      gbl.var$max.x <- gbl.var$max.x + 500
    }
    #   if (config.var$VERBOSE)  cat("gbl.var$max.x", gbl.var$max.x, "\n")
    
    #------------------- DEFINE DATA FOR ONLY ZOOM REGION -----------------------
    ## if (config.var$VERBOSE)  cat("Sort data", gbl.var$sorted.mydata.pos, "\n")
    ## if (config.var$VERBOSE)  cat("Sort  length data", is.vector(gbl.var$sorted.mydata.pos), "\n")
    remove.sorted.mydata.pos <-which(gbl.var$sorted.mydata.pos < gbl.var$min.x | gbl.var$sorted.mydata.pos > gbl.var$max.x)
    # if (config.var$VERBOSE)  cat("Sort oouutt", length(remove.sorted.mydata.pos), "\n")
    if(length(remove.sorted.mydata.pos) > 0){
      sorted.mydata.pos <-  gbl.var$sorted.mydata.pos
      #  if (config.var$VERBOSE)  cat("List gb  ",sorted.mydata.pos,"\n")
      sorted.mydata.pos.zoom <- sorted.mydata.pos[-remove.sorted.mydata.pos]
      #  if (config.var$VERBOSE)  cat("List gb  ",sorted.mydata.pos.zoom,"\n")
    }else {
      sorted.mydata.pos.zoom <- gbl.var$sorted.mydata.pos
    }
    gbl.var$sorted.mydata.pos.zoom <- sorted.mydata.pos.zoom  
    # if (config.var$VERBOSE)  cat("gbl.var$sorted.mydata.pos.zoom ", gbl.var$sorted.mydata.pos.zoom, "\n")
    
    
    #----------------- DEFINE THE HIGHEST VALUE OF Y
    gbl.var$min.y <- min(c(gbl.var$min.y, min(gbl.var$mydata.data$MYDATA.PVAL, na.rm = TRUE)))
    # Ensure that the y scale extends beyond the data region, round to a whole number
    gbl.var$max.y <- ceiling(max(c(gbl.var$max.y, max(gbl.var$mydata.data$MYDATA.PVAL, na.rm = TRUE))))
    
    #DEBUG STATEMENT
    # if (config.var$VERBOSE)  cat("gbl.var$min.y ", gbl.var$min.y, "\n")
    #  if (config.var$VERBOSE)  cat("gbl.var$max.y ", gbl.var$max.y, "\n")
    
    #--------------- DEFINE THE VALUE OF X-AXIS
    gbl.var$cur.exp <- 0
    
    max.min.diff <- gbl.var$max.x - gbl.var$min.x
    exp.test <- TRUE
    
    while(exp.test) {
      if(max.min.diff > 10) {
        gbl.var$cur.exp <- gbl.var$cur.exp + 1
        max.min.diff <- max.min.diff / 10
      }
      else {
        exp.test <- FALSE
      }
    }
    
    gbl.var$r.x <- c(gbl.var$min.x, gbl.var$max.x) # 0 - max(LOC)
    #cat(gbl.var$r.x, "\n")
    
    gbl.var$axis.x <- seq(gbl.var$min.x , gbl.var$max.x, 10^(gbl.var$cur.exp - 1)) # BY Variable
    gbl.var$axis.x[length(gbl.var$axis.x)] <- gbl.var$max.x
    
    
    #DEBUG STATEMENT
    # if (config.var$VERBOSE)  cat("config.var$DISP.LEGEND ", config.var$DISP.LEGEND, "\n")
    # if (config.var$VERBOSE)  cat("gbl.var$max.y ", gbl.var$max.y, "\n")
    
    if(config.var$DISP.LEGEND & !config.var$DISP.CORMATRIXMAP) {
      
      gbl.var$axis.y <- seq(as.integer(gbl.var$min.y), gbl.var$max.y, 1) # BY Variable
      gbl.var$r.y <- c(gbl.var$min.y, gbl.var$max.y + 0.5) # 0 - max(PVAL) for MYDATA
    }
    else {
      gbl.var$axis.y <- seq(as.integer(gbl.var$min.y), gbl.var$max.y, 1) # BY Variable
      gbl.var$r.y <- c(gbl.var$min.y, gbl.var$max.y) # 0 - max(PVAL) for MYDATA
    }
    
    if(gbl.var$cur.sample == 1){
      best.position <- 0
      format <- as.character(gbl.var$split.format[[1]][gbl.var$cur.sample])
      if(is.null(config.var$MYDATA.REF)){
        best.position <- which.max(gbl.var$mydata.data$MYDATA.PVAL)
      }else {
        best.position <- which(gbl.var$mydata.data$MYDATA.NAME == config.var$MYDATA.REF)
      }
      if(grepl("REGION", format)[1] ) {
        gbl.var$mydata.ref.pos <- gbl.var$mydata.data$LOC.START[best.position] + ((gbl.var$mydata.data$LOC.END[best.position] - gbl.var$mydata.data$LOC.START[best.position])/2)
        gbl.var$mydata.ref.pos <- gbl.var$mydata.data$LOC.START[best.position] + ((gbl.var$mydata.data$LOC.END[best.position] - gbl.var$mydata.data$LOC.START[best.position])/2)
      } else {
        gbl.var$mydata.ref.pos <- gbl.var$mydata.data$LOC[best.position]
        gbl.var$mydata.ref.pos <- gbl.var$mydata.data$LOC[best.position]
      }
      gbl.var$mydata.best.position <- best.position
      if (config.var$VERBOSE)  cat("Ref",gbl.var$mydata.ref.pos,"\n")
      
      #-----------------READ MATRICE for CORMATRIX MAP
      # if (config.var$VERBOSE)  cat("config.var$DISP.CORMATRIXMAP ", config.var$DISP.CORMATRIXMAP, "\n")
      #Need to read even if we don't show data
      gbl.var <-read.file.cormatrix(config.var, gbl.var,split.cormatrix.file)
    }
  }
  
  
  #------------------- READ LARGE DATA
  if(!is.null(config.var$MYDATA.LARGE.FILE)){
    large.split.mydata.file=gbl.var$large.split.mydata.file
    
    for(i in 1:length(large.split.mydata.file[[1]])) {
      #-------------------- READ FILE CpG DATA
      gbl.var <-read.file.mydata.large(large.split.mydata.file,config.var, gbl.var,i)
      
      cur.mydata.large.num <- nrow(gbl.var$mydata.large.data)
      #create 1 hash table containing MYDATA names and positions for current and previous sample sets
      #one hash has the name as the key and the other the location because R has no way of looking
      #up hash keys based on values
      format <- as.character(gbl.var$large.split.format[[1]][i])
      if(grepl("REGION", format)[1] ) {
        mydata.large.hash.names.start <- new.env(hash=TRUE)
        mydata.large.hash.names.end <- new.env(hash=TRUE)
        mydata.large.hash.names.pos <-new.env(hash=TRUE)
        mydata.large.hash.pos.names <- hash()
        
        for(p in 1:cur.mydata.large.num) {
          
          localisation <- 0
          
          if( (!is.na(gbl.var$mydata.large.data$LOC.START[p]) & !is.na(gbl.var$mydata.large.data$LOC.END[p])) &  ( gbl.var$mydata.large.data$LOC.START[p] < gbl.var$mydata.large.data$LOC.END[p]  | gbl.var$mydata.large.data$LOC.START[p] == gbl.var$mydata.large.data$LOC.END[p])){
            localisation <- gbl.var$mydata.large.data$LOC.START[p] + ((gbl.var$mydata.large.data$LOC.END[p] -gbl.var$mydata.large.data$LOC.START[p]))/2
          }else {
            stop("Problem in the coordinates of region\n")
          }
          
          if(!is.na(gbl.var$mydata.large.data$MYDATA.NAME[p]) ) {
            mydata.large.hash.names.pos[[ as.character(gbl.var$mydata.large.data$MYDATA.NAME[p]) ]]<- localisation
            #assign(gbl.var$mydata.large.data$MYDATA.NAME[p], localisation, envir=mydata.large.hash.names.pos)
            mydata.large.hash.names.start[[gbl.var$mydata.large.data$MYDATA.NAME[p] ]] <- gbl.var$mydata.large.data$LOC.START[p]
            mydata.large.hash.names.end[[gbl.var$mydata.large.data$MYDATA.NAME[p] ]] <- gbl.var$mydata.large.data$LOC.END[p]
            mydata.large.hash.pos.names[[as.character(localisation)]] <- gbl.var$mydata.large.data$MYDATA.NAME[p]
          }
        }
        gbl.var$mydata.large.hash.names.pos <- mydata.large.hash.names.pos
        gbl.var$mydata.large.hash.names.start <- mydata.large.hash.names.start
        gbl.var$mydata.large.hash.names.end <- mydata.large.hash.names.end
        gbl.var$mydata.large.hash.pos.names <- mydata.large.hash.pos.names
        #  if (config.var$VERBOSE)  cat("test hash",values(gbl.var$mydata.large.hash.names.pos),"\n")
        
      } else {
        mydata.large.hash.names.pos <- new.env(hash=TRUE)
        mydata.large.hash.pos.names <- new.env(hash=TRUE)
        for(p in 1:cur.mydata.large.num) {
          localisation <- gbl.var$mydata.large.data$LOC[p]
          if(!is.na(gbl.var$mydata.large.data$MYDATA.NAME[p]) | !is.na(gbl.var$mydata.large.data$LOC[p])) {
            mydata.large.hash.names.pos[[ as.character(gbl.var$mydata.large.data$MYDATA.NAME[p]) ]]<- localisation
            mydata.large.hash.pos.names[[as.character(localisation)]] <- gbl.var$mydata.large.data$MYDATA.NAME[p]
          }
        }
        gbl.var$mydata.large.hash.names.pos <- mydata.large.hash.names.pos
        gbl.var$mydata.large.hash.pos.names <- mydata.large.hash.pos.names
      }
      #create a variable with only the MYDATA positions
      mydata.large.hash.pos <- ls(envir=gbl.var$mydata.large.hash.pos.names)
      
      #DEBUG STATEMENT
      ## if (config.var$VERBOSE)  cat("length(mydata.hash.pos) ", length(mydata.hash.pos), "\n")
      
      #create a variable with only the sorted MYDATA positions
      gbl.var$sorted.mydata.large.pos <- sort(as.numeric(mydata.large.hash.pos))
      
      #------------------- ZOOM on LARGE DATA -----------------------
      ## if (config.var$VERBOSE)  cat("Sort data", gbl.var$sorted.mydata.large.pos, "\n")
      ## if (config.var$VERBOSE)  cat("Sort  length data", is.vector(gbl.var$sorted.mydata.large.pos), "\n")
      remove.sorted.mydata.large.pos <-which(gbl.var$sorted.mydata.large.pos < gbl.var$min.x | gbl.var$sorted.mydata.large.pos > gbl.var$max.x)
      # if (config.var$VERBOSE)  cat("Sort oouutt", length(remove.sorted.mydata.large.pos), "\n")
      if(length(remove.sorted.mydata.large.pos) > 0){
        sorted.mydata.large.pos <-  gbl.var$sorted.mydata.large.pos
        #  if (config.var$VERBOSE)  cat("List gb  ",sorted.mydata.large.pos,"\n")
        sorted.mydata.large.pos.zoom <- sorted.mydata.large.pos[-remove.sorted.mydata.large.pos]
        #  if (config.var$VERBOSE)  cat("List gb  ",sorted.mydata.large.pos.zoom,"\n")
      }else {
        sorted.mydata.large.pos.zoom <- gbl.var$sorted.mydata.large.pos
      }
      gbl.var$sorted.mydata.large.pos.zoom <- sorted.mydata.large.pos.zoom  
      # if (config.var$VERBOSE)  cat("gbl.var$sorted.mydata.large.pos.zoom ", gbl.var$sorted.mydata.large.pos.zoom, "\n")
      
    }
    
    #----------------- DEFINE THE HIGHEST VALUE OF Y
    #DEBUG STATEMENT
    #  if (config.var$VERBOSE)  cat("gbl.var$min.y ", gbl.var$min.y, "\n")
    # if (config.var$VERBOSE)  cat("gbl.var$max.y ", gbl.var$max.y, "\n")
    if(config.var$LAB.Y == "ln") {
      gbl.var$mydata.large.data$MYDATA.PVAL <- -log(gbl.var$mydata.large.data$MYDATA.PVAL) # LOG Variable
    }
    if(config.var$LAB.Y == "log") {
      gbl.var$mydata.large.data$MYDATA.PVAL <- -log10(gbl.var$mydata.large.data$MYDATA.PVAL) # LOG Variable
    }
    
    gbl.var$min.y <- min(c(gbl.var$min.y, min(gbl.var$mydata.large.data$MYDATA.PVAL, na.rm = TRUE)))
    # Ensure that the y scale extends beyond the data region, round to a whole number
    gbl.var$max.y <- ceiling(max(c(gbl.var$max.y, max(gbl.var$mydata.large.data$MYDATA.PVAL, na.rm = TRUE))))
    
    #DEBUG STATEMENT
    # if (config.var$VERBOSE)  cat("gbl.var$min.y ", gbl.var$min.y, "\n")
    # if (config.var$VERBOSE)  cat("gbl.var$max.y ", gbl.var$max.y, "\n")
    if(config.var$DISP.LEGEND & !config.var$DISP.CORMATRIXMAP) {
      
      gbl.var$axis.y <- seq(as.integer(gbl.var$min.y), gbl.var$max.y, 1) # BY Variable
      gbl.var$r.y <- c(gbl.var$min.y, gbl.var$max.y + 0.5) # 0 - max(PVAL) for MYDATA
    }
    else {
      gbl.var$axis.y <- seq(as.integer(gbl.var$min.y), gbl.var$max.y, 1) # BY Variable
      gbl.var$r.y <- c(gbl.var$min.y, gbl.var$max.y) # 0 - max(PVAL) for MYDATA
    }
    
  }
  
  #----- CREATE LIST of BIOFEATURES -------------------------
  
  if(!is.null(config.var$BIOFEAT.USER.FILE)) {
    gbl.var <- create.tracks.user(config.var, gbl.var)
  }
  
  
  #----- FIX VALUES -------------------------
  
  fix.var.generic <- fix.values.generic(config.var, gbl.var)
  
  config.var <- fix.var.generic$config.var
  gbl.var <- fix.var.generic$gbl.var
  
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("FINISH RETRIEVE.DATA\n")
  
  return(list(config.var = config.var, gbl.var = gbl.var))
}


#-------------------FIX VALUES GENERIC----------------------------
fix.values.generic <- function(config.var, gbl.var) {
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("START FIX.VALUES.GENERIC\n")
  # if (config.var$VERBOSE)  cat("PVALUE",gbl.var$mydata.data[1,4],"\n")
  
  if(config.var$LAB.Y == "ln") {
    
    config.var$PVAL.THRESHOLD <- -log(config.var$PVAL.THRESHOLD)
    
    #DEBUG STATEMENT
    # if (config.var$VERBOSE)  cat("PVAL.THRESHOLD ", config.var$DISP.PVAL.THRESHOLD, "\n")
    # if (config.var$VERBOSE)  cat("PVAL.THRESHOLD ", config.var$PVAL.THRESHOLD, "\n")
    
    if(gbl.var$pval.flag) {
      if(config.var$DISP.PVAL.THRESHOLD != 1) {
        config.var$DISP.PVAL.THRESHOLD <- -log(config.var$DISP.PVAL.THRESHOLD)
        gbl.var$pval.flag <- FALSE
      } else {
        config.var$DISP.PVAL.THRESHOLD <- -1
        gbl.var$pval.flag <- FALSE
      }
    }
  } else if (config.var$LAB.Y == "log") {
    
    config.var$PVAL.THRESHOLD <- -log10(config.var$PVAL.THRESHOLD)
    
    #DEBUG STATEMENT
    # if (config.var$VERBOSE)  cat("DISP.PVAL.THRESHOLD ", config.var$DISP.PVAL.THRESHOLD, "\n")
    # if (config.var$VERBOSE)  cat("PVAL.THRESHOLD ", config.var$PVAL.THRESHOLD, "\n")
    
    if(gbl.var$pval.flag) {
      if(config.var$DISP.PVAL.THRESHOLD != 1) {
        config.var$DISP.PVAL.THRESHOLD <- -log10(config.var$DISP.PVAL.THRESHOLD)
        gbl.var$pval.flag <- FALSE
      } else {
        config.var$DISP.PVAL.THRESHOLD <- -1
        gbl.var$pval.flag <- FALSE
      }
    }
  } else {
    stop("Invalid LAB.Y: ", config.var$LAB.Y, ". Specify either 'log' or 'ln'.\n")
  }
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("FINISH FIX.VALUES.GENERIC\n")
  
  return(list(config.var = config.var, gbl.var = gbl.var))
}

#-------------------FIX VALUES for MYDATA----------------------------
fix.values <- function(config.var, gbl.var) {
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("START FIX.VALUES\n")
  # if (config.var$VERBOSE)  cat("PVALUE",gbl.var$mydata.data[1,4],"\n")
  
  if(config.var$LAB.Y == "ln") {
    
    gbl.var$mydata.data$MYDATA.PVAL <- -log(gbl.var$mydata.data$MYDATA.PVAL) # LOG Variable
    
  } else if (config.var$LAB.Y == "log") {
    
    gbl.var$mydata.data$MYDATA.PVAL <- -log10(gbl.var$mydata.data$MYDATA.PVAL) # LOG Variable
    
  } else {
    stop("Invalid LAB.Y: ", config.var$LAB.Y, ". Specify either 'log' or 'ln'.\n")
  }
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("FINISH FIX.VALUES\n")
  
  return(list(config.var = config.var, gbl.var = gbl.var))
}

#-------------------FIX VALUES for MYDATA----------------------------
fix.values.large <- function(config.var, gbl.var) {
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("START FIX.VALUES\n")
  # if (config.var$VERBOSE)  cat("PVALUE",gbl.var$mydata.data[1,4],"\n")
  
  if(config.var$LAB.Y == "ln") {
    
    gbl.var$mydata.large.data$MYDATA.PVAL <- -log(gbl.var$mydata.large.data$MYDATA.PVAL) # LOG Variable
    
  } else if (config.var$LAB.Y == "log") {
    
    gbl.var$mydata.large.data$MYDATA.PVAL <- -log10(gbl.var$mydata.large.data$MYDATA.PVAL) # LOG Variable
    
  } else {
    stop("Invalid LAB.Y: ", config.var$LAB.Y, ". Specify either 'log' or 'ln'.\n")
  }
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("FINISH FIX.VALUES\n")
  
  return(list(config.var = config.var, gbl.var = gbl.var))
}



#------------------- CHECK format CpG data -------
# option 1 is for reference data, option 2 is for large data
check.format.mydata <- function(gbl.var,option,numfile){
  #DEBUG STATEMENT
  
  #DEBUG STATEMENT
  
  if(option == 1) {
    if (gbl.var$verbose)  cat("START CHECK.FORMAT.MYDATA\n")
    mydata.test <- gbl.var$mydata.data
    format <- as.character(gbl.var$split.format[[1]][numfile])
  }else {
    if (gbl.var$verbose)  cat("START CHECK.FORMAT.MYDATA.LARGE\n")
    mydata.test <- gbl.var$mydata.large.data
    format <- as.character(gbl.var$large.split.format[[1]][numfile])
  }
  
  if (gbl.var$verbose)  cat("mydata.test ", dim(mydata.test), "\n")
  if (gbl.var$verbose)  cat("mydata.test ", colnames(mydata.test), "\n")
  if (gbl.var$verbose)  cat("mydata.format ", format, "\n")
  
  #check for data columns and compatibility with display type option
  if (is.null(mydata.test)) {
    stop("Missing MYDATA data file \n")
  } 
  if (is.null(mydata.test$MYDATA.NAME)) {
    stop("Missing MYDATA data file column MYDATA.NAME\n")
  } 
  if (length(grep("^[0-9]", mydata.test$MYDATA.NAME)) > 0) {
    stop("MYDATA names cannot start with numbers.\n")
  } 
  if ((format == "DTR" | format == "DTR_ASSO" | format == "SITE" | format == "SITE_ASSO")) {
    if(is.null(mydata.test$LOC)) {
      stop("Missing MYDATA data file column LOC\n")
    }
    if ((length(grep("^[a-zA-Z]",mydata.test$LOC)) > 0) ){
      stop("Missing column LOC has to be number\n")
    }
  } else if (format == "REGION" | format == "REGION_ASSO") {
    if(is.null(mydata.test$LOC.START)) {
      stop("Missing MYDATA data file column LOC.START\n")
    }
    if ((length(grep("^[a-zA-Z]",mydata.test$LOC.START)) > 0) ){
      stop("Missing column LOC.START has to be number\n")
    }
    if (is.null(mydata.test$LOC.END)) {
      stop("Missing MYDATA data file column LOC.END\n")
    } 
    if ((length(grep("^[a-zA-Z]",mydata.test$LOC.END)) > 0) ){
      stop("Missing column LOC.END has to be number\n")
    }
  } 
  if (is.null(mydata.test$MYDATA.PVAL)) {
    stop("Missing MYDATA data file column MYDATA.PVAL\n")
  }  else if (length(grep("^[a-zA-Z]",mydata.test$MYDATA.PVAL)) > 0){
    stop("Missing column MYDATA.PVAL has to be number\n")
  } 
  if (format == "REGION_ASSO" | format == "DTR_ASSO" |  format == "SITE_ASSO") {
    if(is.null(mydata.test$MYDATA.ASSO) ) {
      stop("Missing MYDATA data file column MYDATA.ASSO\n")
    }
    #  if (gbl.var$verbose)  cat("BUG Association ",length(mydata.test$MYDATA.ASSO)," \n")
    # if (gbl.var$verbose)  cat("BUG Association ",length(grep("[0-9]",mydata.test$MYDATA.ASSO))," \n")
    if ((length(grep("[a-zA-Z]",mydata.test$MYDATA.ASSO)) > 0)){
      stop("Missing column MYDATA.ASSO has to be number\n")
    } else if ((length(grep("[0-9|-|+]",mydata.test$MYDATA.ASSO)) == length(mydata.test$MYDATA.ASSO))){
      pos <- which(as.numeric(mydata.test$MYDATA.ASSO) > 0)
      neg <- which(as.numeric(mydata.test$MYDATA.ASSO) < 0)
      noth <- which(as.numeric(mydata.test$MYDATA.ASSO) == 0)
      mydata.test$MYDATA.ASSO[pos] <- "+"
      mydata.test$MYDATA.ASSO[neg] <- "-"
      # mydata.test$MYDATA.ASSO[noth] <- "0"
      
      if(option == 1) {
        gbl.var$mydata.data <- mydata.test
      }else {
        gbl.var$mydata.large.data <- mydata.test
      }
    } else if ((length(grep("[+|-]",mydata.test$MYDATA.ASSO)) == length(mydata.test$MYDATA.ASSO))){  
    } else {
      dd <- length(grep("[0-9|+|-]",mydata.test$MYDATA.ASSO))
      da <- length(mydata.test$MYDATA.ASSO)
      stop("Do not recognize column MYDATA.ASSO ", dd, " vs ",da," \n")
    }
    
  } 
  if (is.null(mydata.test$CHROMOSOME)) {
    stop("Missing MYDATA data file column CHROMOSOME\n")
  } else if (!grepl( "^chr",mydata.test$CHROMOSOME)[1]) {
    if (gbl.var$verbose)  cat("Missing column CHROMOSOME should be at UCSC format\n")
    mydata.test$CHROMOSOME <- paste("chr",mydata.test$CHROMOSOME,sep="")
    if(option == 1) {
      gbl.var$mydata.data <- mydata.test
    }else {
      gbl.var$mydata.large.data <- mydata.test
    }
    
  }  
  
  # if (gbl.var$verbose)  cat("format",is.character(format),"\n")
  # if (format == "REGION_ASSO" | format == "DTR_ASSO" | format == "SITE_ASSO") {
  #    if (gbl.var$verbose)  cat("Association ",gbl.var$mydata.large.data$Association[1]," \n") 
  # }
  
  #Sort DATA from localisation
  if(option == 1) {
    mydata.test.nosort <- gbl.var$mydata.data 
    mydata.test.sort <- mydata.test.nosort
    # if (gbl.var$verbose)  cat("test",grepl( "REGION",format)[1],"\n")
    if (grepl( "REGION",format)[1]) {
      mydata.test.sort <- mydata.test.nosort[order(mydata.test.nosort$LOC.START),]
    } else {
      mydata.test.sort <- mydata.test.nosort[order(mydata.test.nosort$LOC),]
    }
    gbl.var$mydata.data <- mydata.test.sort
  }else {
    mydata.test.nosort <- gbl.var$mydata.large.data 
    mydata.test.sort <- mydata.test.nosort
    if (grepl( "REGION",format)[1]) {
      mydata.test.sort <- mydata.test.nosort[order(mydata.test.nosort$LOC.START),]
    } else {
      mydata.test.sort <- mydata.test.nosort[order(mydata.test.nosort$LOC),]
    }
    gbl.var$mydata.large.data <- mydata.test.sort
  }
  
  #DEBUG STATEMENT
  if (gbl.var$verbose)  cat("FINISH CHECK.FORMAT.MYDATA\n")
  return(gbl.var)
}

#------------------- READ file CpG data -------
read.file.mydata <- function(split.mydata.file,config.var, gbl.var,numfile){
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("START READ.FILE.MYDATA ",numfile,"\n")
  #  if (config.var$VERBOSE)  cat("Format",config.var$MYDATA.FORMAT,"\n")
  
  #  split.mydata.file <-gbl.var$split.mydata.file
  # numfile <- gbl.var$cur.sample
  general.data.raw <- read.delim(split.mydata.file[[1]][numfile], header=TRUE, sep="\t", as.is=TRUE, blank.lines.skip = TRUE, fill=TRUE)
  gbl.var$general.data <- general.data.raw
  
  #----------- MYDATA table with pvalue
  if (config.var$MYDATA.FORMAT == "DTR") {
    gbl.var$mydata.data <- gbl.var$general.data[,c(1,2,7,14)]
    colnames(gbl.var$mydata.data)<-c("MYDATA.NAME","CHROMOSOME","LOC","MYDATA.PVAL")
  } else if (config.var$MYDATA.FORMAT == "DTR_ASSO") {
    gbl.var$mydata.data <- gbl.var$general.data[,c(1,2,7,14,15)]
    colnames(gbl.var$mydata.data)<-c("MYDATA.NAME","CHROMOSOME","LOC","MYDATA.PVAL","MYDATA.ASSO")
  } else if (config.var$MYDATA.FORMAT == "SITE") {
    gbl.var$mydata.data <- gbl.var$general.data
    colnames(gbl.var$mydata.data)<-c("MYDATA.NAME","CHROMOSOME","LOC","MYDATA.PVAL")
  } else if (config.var$MYDATA.FORMAT == "SITE_ASSO") {
    gbl.var$mydata.data <- gbl.var$general.data
    colnames(gbl.var$mydata.data)<-c("MYDATA.NAME","CHROMOSOME","LOC","MYDATA.PVAL","MYDATA.ASSO")
  } else if (config.var$MYDATA.FORMAT == "REGION") {
    gbl.var$mydata.data <- gbl.var$general.data
    colnames(gbl.var$mydata.data)<-c("MYDATA.NAME","CHROMOSOME","LOC.START","LOC.END","MYDATA.PVAL")
  } else if (config.var$MYDATA.FORMAT == "REGION_ASSO") {
    gbl.var$mydata.data <- gbl.var$general.data
    colnames(gbl.var$mydata.data)<-c("MYDATA.NAME","CHROMOSOME","LOC.START","LOC.END","MYDATA.PVAL","MYDATA.ASSO")
  } else {
    stop("This format does not exist\n")
  }
  
  #------------------- CHECK FILE CpG DATA
  gbl.var<-check.format.mydata(gbl.var,1,numfile)
  
  if(!is.null(gbl.var$mydata.data$MYDATA.PVAL[which(is.na(gbl.var$mydata.data$MYDATA.PVAL))])){
    gbl.var$mydata.data$MYDATA.PVAL[which(is.na(gbl.var$mydata.data$MYDATA.PVAL))] <- 1
  }
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("FINISH READ.FILE.MYDATA ",numfile,"\n")
  return(gbl.var)
}

#------------------- READ LARGE FILE OF POSITION and PVALUE-------
read.file.mydata.large <- function(large.split.mydata.file, config.var, gbl.var,numfile.large){
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("START READ.FILE.MYDATA.LARGE ",numfile.large," \n")
  # if (config.var$VERBOSE)  cat("Format",config.var$MYDATA.LARGE.FORMAT,"...", numfile.large,"\n")
  
  general.data.large.raw.tmp <- read.delim(large.split.mydata.file[[1]][numfile.large], header=TRUE, sep="\t", as.is=TRUE, blank.lines.skip = TRUE, fill=TRUE)
  
  large.format <- gbl.var$large.split.format[[1]][numfile.large]
  
  #----------- MYDATA table with pvalue
  if (large.format == "DTR") {
    general.data.large.raw <- general.data.large.raw.tmp[,c(1,2,7,14)]
    colnames(general.data.large.raw)<-c("MYDATA.NAME","CHROMOSOME","LOC","MYDATA.PVAL")
  } else if (large.format == "DTR_ASSO") {
    general.data.large.raw <- general.data.large.raw.tmp[,c(1,2,7,14,15)]
    colnames(general.data.large.raw)<-c("MYDATA.NAME","CHROMOSOME","LOC","MYDATA.PVAL","MYDATA.ASSO")
  } else if (large.format == "SITE") {
    general.data.large.raw <- general.data.large.raw.tmp
    colnames(general.data.large.raw)<-c("MYDATA.NAME","CHROMOSOME","LOC","MYDATA.PVAL")
  } else if (large.format == "SITE_ASSO") {
    general.data.large.raw <- general.data.large.raw.tmp
    colnames(general.data.large.raw)<-c("MYDATA.NAME","CHROMOSOME","LOC","MYDATA.PVAL","MYDATA.ASSO")
  }else if (large.format == "REGION") {
    general.data.large.raw <- general.data.large.raw.tmp
    colnames(general.data.large.raw)<-c("MYDATA.NAME","CHROMOSOME","LOC.START","LOC.END","MYDATA.PVAL")
  } else if (large.format == "REGION_ASSO") {
    general.data.large.raw <- general.data.large.raw.tmp
    colnames(general.data.large.raw)<-c("MYDATA.NAME","CHROMOSOME","LOC.START","LOC.END","MYDATA.PVAL","MYDATA.ASSO")
  }else{
    stop("This format does not exist\n")
  }
  
  
  if (!grepl( "^chr",general.data.large.raw$CHROMOSOME)[1]) {
    if (config.var$VERBOSE)  cat("Missing column CHROMOSOME should be at UCSC format\n")
    gbl.var$mydata.large.data <- general.data.large.raw[which(general.data.large.raw$CHROMOSOME == gsub("^chr", "", gbl.var$mydata.chr)),]
  }else {
    gbl.var$mydata.large.data <- general.data.large.raw[which(general.data.large.raw$CHROMOSOME == gbl.var$mydata.chr),]
    
  }
  
  #------------------- CHECK FILE CpG DATA
  gbl.var<-check.format.mydata(gbl.var,2,numfile.large)
  
  if(!is.null(gbl.var$mydata.large.data$MYDATA.PVAL[which(is.na(gbl.var$mydata.large.data$MYDATA.PVAL))])){
    gbl.var$mydata.large.data$MYDATA.PVAL[which(is.na(gbl.var$mydata.large.data$MYDATA.PVAL))] <- 1
  }
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("FINISH READ.FILE.MYDATA.LARGE \n")
  return(gbl.var)
}


#------------------- READ CORMATRIX FILE -------
read.file.cormatrix <- function(config.var, gbl.var,split.cormatrix.file=NULL){
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("START READ.FILE.COMATRIX\n")
  
  if(is.null(split.cormatrix.file)){
    #  if (config.var$VERBOSE)  cat("CORMATRIX empty \n")
    stop("EMPTY CORMATRIX data file \n")
  } else {
      if (config.var$VERBOSE)  cat("Format data",config.var$CORMATRIX.FORMAT,"\n")
    
    #----------- MYDATA table with pvalue
    if (config.var$CORMATRIX.FORMAT == "DTR_CORMATRIX") {
      gbl.var<-compute.cormatrix(config.var, gbl.var)
      
    } else if (config.var$CORMATRIX.FORMAT == "CORMATRIX") {
      cormatrix.data.raw<- read.delim(split.cormatrix.file[[1]][1], sep="\t", header=TRUE, as.is=TRUE, blank.lines.skip = TRUE, fill=TRUE)
      #--- Check the same size cormatrix and data
      if((nrow(gbl.var$mydata.data) == nrow(cormatrix.data.raw)) & (nrow(gbl.var$mydata.data) == ncol(cormatrix.data.raw))){
       # cormatrix.data.raw[match(gbl.var$mydata.data$MYDATA.NAME, colnames(cormatrix.data.raw)),]
       #  cormatrix.data.raw[gbl.var$mydata.data$MYDATA.NAME,]
        if(identical(as.vector(gbl.var$mydata.data$MYDATA.NAME),as.vector(colnames(cormatrix.data.raw)))) {
          gbl.var$cormatrix.data <-as.matrix(cormatrix.data.raw)
        } else {
          stop("Invalid MYDATA data file: ", config.var$MYDATA.FILE, " and CORMATRIX ",split.cormatrix.file[[1]][1] ," do not have the same omic sites/regions or not in the same order, need to put the ascendant order of their start position in the two files.\n Error: ",all.equal(as.vector(gbl.var$mydata.data$MYDATA.NAME),as.vector(colnames(cormatrix.data.raw))),"\n")
        }
      }else {
        stop("Invalid MYDATA data file: ", config.var$MYDATA.FILE, " and CORMATRIX ",split.cormatrix.file[[1]][1] ," do not have the same size\n")
      }
      
    } else if (config.var$CORMATRIX.FORMAT == "RAW_REV"){
      # if (config.var$VERBOSE)  cat("Format data",config.var$CORMATRIX.FORMAT,"\n")
      matrix.data.raw_rot<- read.delim(split.cormatrix.file[[1]][1], header=TRUE, sep="\t", as.is=TRUE, blank.lines.skip = TRUE, fill=TRUE)
      matrix.data.raw <- matrix.data.raw_rot
      matrix.data.raw[match(gbl.var$mydata.data$MYDATA.NAME, matrix.data.raw[1,]),]
      matrix.data.raw <-matrix.data.raw[,-1]
      gbl.var$matrix.data <-matrix.data.raw
      nr=nrow(matrix.data.raw)
      nr1=nrow(gbl.var$mydata.data)
      #if (config.var$VERBOSE)  cat("Matrice size",nr,"\n")
      #if (config.var$VERBOSE)  cat("Matrice size",nr1,"\n")
      #--- Check the same size cormatrix and data
      if(nrow(gbl.var$mydata.data) == nrow(matrix.data.raw) ){
        gbl.var$matrix.data <-matrix.data.raw
        gbl.var<-compute.cormatrix(config.var, gbl.var)
      }else {
        stop("Invalid MYDATA data file: ", config.var$MYDATA.FILE," value (" ,nr1,") and CORMATRIX ",split.cormatrix.file[[1]][1] ," value (" ,nr,") do not have the same size\n")
      }
      
    }else if (config.var$CORMATRIX.FORMAT == "RAW"){
      matrix.data.raw_rot<- read.delim(split.cormatrix.file[[1]][1], header=TRUE, sep="\t", as.is=TRUE, blank.lines.skip = TRUE, fill=TRUE)
      matrix.data.raw <- t(matrix.data.raw_rot)
      matrix.data.raw[match(gbl.var$mydata.data$MYDATA.NAME, matrix.data.raw[1,]),]
      gbl.var$matrix.data <-matrix.data.raw
      nr=nrow(matrix.data.raw)
      nr1=nrow(gbl.var$mydata.data)
       if (config.var$VERBOSE)  cat("Matrice size",nr,"\n")
       if (config.var$VERBOSE)  cat("Matrice size",nr1,"\n")
      #--- Check the same size cormatrix and data
      if(nrow(gbl.var$mydata.data) == nrow(matrix.data.raw) ){
        gbl.var$matrix.data <-matrix.data.raw
        gbl.var<-compute.cormatrix(config.var, gbl.var)
      }else {
        stop("Invalid MYDATA data file: ", config.var$MYDATA.FILE," value (" ,nr1,") and CORMATRIX ",split.cormatrix.file[[1]][1] ," value (" ,nr,") do not have the same size\n")
      }
      
    } else if (config.var$CORMATRIX.FORMAT == "DTR_RAW"){
      matrix.data.raw.tmp <- gbl.var$general.data[,-c(1:19)]
      nr=nrow(matrix.data.raw.tmp)
      nr1=nrow(gbl.var$mydata.data)
      # if (config.var$VERBOSE)  cat("Matrice size",nr,"\n")
      # if (config.var$VERBOSE)  cat("Matrice size",nr1,"\n")
      #--- Check the same size cormatrix and data
      if(nrow(gbl.var$mydata.data) == nrow(matrix.data.raw) ){
        gbl.var$matrix.data <-matrix.data.raw.tmp
        gbl.var<-compute.cormatrix(config.var, gbl.var)
      }else {
        stop("Invalid MYDATA data file: ", config.var$MYDATA.FILE, " and MATRIX do not have the same size\n")
      }
    }
    
    # checking size of data
    num.mydata<-dim(gbl.var$mydata.data)[1]
    row.LD<-dim(gbl.var$cormatrix.data)[1]
    col.LD<-dim(gbl.var$cormatrix.data)[2]
    
    if (config.var$VERBOSE)  cat("Matrice size num.mydat a",num.mydata,"\n")
    if (config.var$VERBOSE)  cat("Matrice size row.LD ",row.LD,"\n")
    if (config.var$VERBOSE)  cat("Matrice size col.LD ",col.LD,"\n")
    
    if (row.LD != num.mydata | col.LD != num.mydata ) {
      stop("Not the same number of row or column between the correlation matrice and MYDATA data\n")
    }  
    
    rownames(gbl.var$cormatrix.data) <- gbl.var$mydata.data$MYDATA.NAME
    #--------- REDUCE the matrix to triangle ---
    gbl.var$cormatrix.data.full <- gbl.var$cormatrix.data
    gbl.var$cormatrix.data[lower.tri(gbl.var$cormatrix.data)] <- NA
    diag(gbl.var$cormatrix.data) <- NA
    
    #--------- VECTOR for the reference gene ---
    best.position <- gbl.var$mydata.best.position
    gbl.var$ref <-gbl.var$cormatrix.data.full[,best.position]
    if (config.var$VERBOSE)  cat("Ref position",best.position,"\n")
  }
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("FINISH READ.FILE.COMET\n")
  return(gbl.var)
}

#------------------- Compute correlation between methylation  -------
compute.cormatrix <- function(config.var, gbl.var){
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("START COMPUTE.CORMATRIX\n")
  
  matrix.data.raw <- gbl.var$matrix.data
  if(!(nrow(gbl.var$mydata.data) == nrow(matrix.data.raw) )){
    stop("Invalid MYDATA data file: ", config.var$MYDATA.FILE, " and MATRIX do not have the same size\n")
  }
  
  if (config.var$CORMATRIX.METHOD == "spearman") {
    gbl.var$cormatrix.data <- cor(t(matrix.data.raw),method="spearman")
  } else if (config.var$CORMATRIX.METHOD == "pearson") {
    gbl.var$cormatrix.data <- cor(t(matrix.data.raw),method="pearson")
  } else if (config.var$CORMATRIX.METHOD == "kendall") {
    gbl.var$cormatrix.data <- cor(t(matrix.data.raw),method="kendall")
  }else {
    stop("Invalid CORMATRIX method : ", config.var$CORMATRIX.METHOD, "\n")
  }
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("FINISH COMPUTE.CORMATRIX\n")
  return(gbl.var)
}

#------------------- Compute list of tracks at GVIZ format  -------
createList.trackUser <- function(config.var, gbl.var){
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("START CREATELIST.TRACKUSER\n")
  split.biofeature.data.user.file <- gbl.var$split.biofeature.data.user.file
  split.biofeature.data.user.type <- gbl.var$split.biofeature.data.user.type
  split.biofeature.data.user.type.plot <- gbl.var$split.biofeature.data.user.type.plot
  
  listtracks_gvizuser <-list()
  for(i in 1:length(split.biofeature.data.user.file[[1]])) {
    cur.sample.biofeat <- i
    biofile <- split.biofeature.data.user.file[[1]]
    if (split.biofeature.data.user.type[[1]][cur.sample.biofeat] == "GeneRegion") {
      gTrack <- GeneRegionTrack(range = biofile, genome = gbl.var$mydata.gen,
                                name = "GeneRegion", chromosome = gbl.var$mydata.chr)
      if(length(listtracks_gvizuser) == 0) {
        listtracks_gvizuser <- list(gTrack)
      } else {
        listtracks_gvizuser <- c(listtracks_gvizuser,gTrack)
      }
    }
    
    if (split.biofeature.data.user.type[[1]][cur.sample.biofeat] == "Annotation") {
      aTrack <- AnnotationTrack(range = biofile, genome = gbl.var$mydata.gen,
                                name = "Annotation", chromosome = gbl.var$mydata.chr)
      if(length(listtracks_gvizuser) == 0) {
        listtracks_gvizuser <- list(aTrack)
      } else {
        listtracks_gvizuser <- c(listtracks_gvizuser,aTrack)
      }
    }
    
    if (split.biofeature.data.user.type[[1]][cur.sample.biofeat] == "Data"){
      typlePlot <- split.biofeature.data.user.type.plot[[1]][cur.sample.biofeat]
      dTrack <- DataTrack(range = biofile, genome = gbl.var$mydata.gen, name = "Plot Data",
                          type = typlePlot, chromosome = gbl.var$mydata.chr)
      
      if(length(listtracks_gvizuser) == 0) {
        listtracks_gvizuser <- list(dTrack)
      } else {
        listtracks_gvizuser <- c(listtracks_gvizuser,dTrack)
      }
    }
    
  }
  
  gbl.var$listtracks_user <- listtracks_gvizuser
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("FINISH CREATELIST.TRACKUSER\n")
  return(gbl.var)
}

#-------------------SET PARAMETERS----------------------------

set.image.parameters <- function(config.var, gbl.var) {
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("START SET.IMAGE.PARAMETERS\n")
  
  sec.cex.factor <- 1
  
  #DEBUG STATEMENT
  # if (config.var$VERBOSE)  cat("IMAGE.SIZE ", config.var$IMAGE.SIZE, "\n")
  
  if(config.var$IMAGE.SIZE == 3.5) {
    font.size <- 4
    line.width <- 0.5
    sec.cex.factor <- 1.5
    cex.factor.symbol <- 0.25
  } else if(config.var$IMAGE.SIZE == 7) {
    font.size <- 9
    line.width <- 1
    sec.cex.factor <- 2
    cex.factor.symbol <- 0.5
  } else {
    stop("Invalid image size: ", config.var$IMAGE.SIZE, "\n")
  }
  
  if(gbl.var$mydata.num > 0 & gbl.var$mydata.num <= 20 ) {
    cex.factor <- 0.54 * sec.cex.factor
  } else if(gbl.var$mydata.num > 20 & gbl.var$mydata.num <= 60 ) {
    cex.factor <- 0.54 * sec.cex.factor
  } else if(gbl.var$mydata.num > 60 & gbl.var$mydata.num <= 90) {
    cex.factor <- 0.4 * sec.cex.factor
  } else if(gbl.var$mydata.num > 90 & gbl.var$mydata.num <= 120) {
    cex.factor <- 0.2 * sec.cex.factor
  } else if(gbl.var$mydata.num > 120) {
    cex.factor <- 0.1 * sec.cex.factor
    line.width <- 0.5
  } else {
    stop("Invalid image size: ", gbl.var$mydata.num, "\n")
  }
  
  #DEBUG STATEMENT
  # if (config.var$VERBOSE)  cat("gbl.var$mydata.num ", gbl.var$mydata.num, "\n")
  # if (config.var$VERBOSE)  cat("font.size ", font.size, "\n")
  # if (config.var$VERBOSE)  cat("line.width ", line.width, "\n")
  # if (config.var$VERBOSE)  cat("cex.factor.symbol ", cex.factor.symbol, "\n")
  # if (config.var$VERBOSE)  cat("cex.factor ", cex.factor, "\n")
  # if (config.var$VERBOSE)  cat("sec.cex.factor ", sec.cex.factor, "\n")
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("FINISH SET.IMAGE.PARAMETERS\n")
  
  return(list(font.size = font.size, line.width = line.width, cex.factor.symbol = cex.factor.symbol, cex.factor = cex.factor))
}
