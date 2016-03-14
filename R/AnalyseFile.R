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
# Version : 0.99.9
###########################################################################

#------------------READ in CONFiguration file------------------
read.config <- function(config.file, config.var) {
  
  #DEBUG STATEMENT
  if (config.var$verbose)   cat("START READ.CONFIG\n")
  
  raw.config <- read.table(config.file, header=FALSE, as.is=TRUE, fill=TRUE, 
                           sep="\t", blank.lines.skip = TRUE)
  
  for(i in 1:length(raw.config$V1)) {
    tmp <- strsplit(raw.config$V1[i], "=")
    
    if(identical("mydata.file", tmp[[1]][1])) {
      config.var$mydata.file <- tmp[[1]][2]
    }
    
    if(identical("mydata.type", tmp[[1]][1])) {
      config.var$mydata.type <- tmp[[1]][2]
    }
    
    if(identical("mydata.format", tmp[[1]][1])) {
      config.var$mydata.format <- tmp[[1]][2]
    }
    
    if(identical("mydata.large.file", tmp[[1]][1])) {
      config.var$mydata.large.file <- tmp[[1]][2]
    }
    
    if(identical("mydata.large.format", tmp[[1]][1])) {
      config.var$mydata.large.format <- tmp[[1]][2]
    }
    
    if(identical("mydata.large.type", tmp[[1]][1])) {
      config.var$mydata.large.type <- tmp[[1]][2]
    }
    
    if(identical("cormatrix.file", tmp[[1]][1])) {
      config.var$cormatrix.file <- tmp[[1]][2]
    }
    
    if(identical("cormatrix.method", tmp[[1]][1])) {
      config.var$cormatrix.method <- tmp[[1]][2]
    }
    
    if(identical("cormatrix.format", tmp[[1]][1])) {
      config.var$cormatrix.format <- tmp[[1]][2]
    }
    
    if(identical("cormatrix.type", tmp[[1]][1])) {
      config.var$cormatrix.type <- tmp[[1]][2]
    }
    
    if(identical("cormatrix.conf.level", tmp[[1]][1])) {
      config.var$cormatrix.conf.level <- as.numeric(tmp[[1]][2])
    }
    
    if(identical("cormatrix.sig.level", tmp[[1]][1])) {
      config.var$cormatrix.sig.level <- as.numeric(tmp[[1]][2])
    }
    
    if(identical("cormatrix.adjust", tmp[[1]][1])) {
      config.var$cormatrix.adjust <- tmp[[1]][2]
    }
    
    if(identical("cormatrix.output", tmp[[1]][1])) {
      config.var$cormatrix.output <- tmp[[1]][2]
    }
    
    if(identical("cormatrix.color.scheme", tmp[[1]][1])) {
      config.var$cormatrix.color.scheme <- tmp[[1]][2]
    }
    
    if(identical("mydata.ref", tmp[[1]][1])) {
      config.var$mydata.ref <- tmp[[1]][2]
    }
    
    if(identical("start", tmp[[1]][1])) {
      config.var$start <- as.numeric(tmp[[1]][2])
    }
    
    if(identical("end", tmp[[1]][1])) {
      config.var$end <- as.numeric(tmp[[1]][2])
    }
    
    if(identical("zoom", tmp[[1]][1])) {
      config.var$zoom <- as.logical(tmp[[1]][2])
    }
    
    if(identical("biofeat.user.file", tmp[[1]][1])) {
      config.var$biofeat.user.file <- tmp[[1]][2]
    }
    
    if(identical("biofeat.user.type", tmp[[1]][1])) {
      config.var$biofeat.user.type <- tmp[[1]][2]
    }
    
    if(identical("biofeat.user.type.plot", tmp[[1]][1])) {
      config.var$biofeat.user.type.plot <- tmp[[1]][2]
    }
    
    if(identical("lab.Y", tmp[[1]][1])) {
      config.var$lab.Y <- tmp[[1]][2]
    }
    
    if(identical("symbols", tmp[[1]][1])) {
      config.var$symbols <- tmp[[1]][2]
    }
    
    if(identical("symbols.large", tmp[[1]][1])) {
      config.var$symbols.large <- tmp[[1]][2]
    }
    
    if(identical("sample.labels", tmp[[1]][1])) {
      config.var$sample.labels <- tmp[[1]][2]
    }
    
    if(identical("sample.labels.large", tmp[[1]][1])) {
      config.var$sample.labels.large <- tmp[[1]][2]
    }
    
    if(identical("pval.threshold", tmp[[1]][1])) {
      config.var$pval.threshold <- as.numeric(tmp[[1]][2])
    }
    
    if(identical("pval.threshold.2", tmp[[1]][1])) {
      config.var$pval.threshold.2 <- as.numeric(tmp[[1]][2])
    }
    
    if(identical("disp.pval.threshold", tmp[[1]][1])) {
      config.var$disp.pval.threshold <- as.logical(tmp[[1]][2])
    }
    
    if(identical("disp.color.ref", tmp[[1]][1])) {
      config.var$disp.color.ref <- as.logical(tmp[[1]][2])
    }
    
    if(identical("use.colors", tmp[[1]][1])) {
      config.var$use.colors <- as.logical(tmp[[1]][2])
    }
    
    if(identical("color.list", tmp[[1]][1])) {
      config.var$color.list <- tmp[[1]][2]
    }
    
    if(identical("color.list.large", tmp[[1]][1])) {
      config.var$color.list.large <- tmp[[1]][2]
    }
    
    if(identical("disp.association", tmp[[1]][1])) {
      #no format logical because of reading as string with comma
      config.var$disp.association <- tmp[[1]][2]
    }
    
    if(identical("disp.association.large", tmp[[1]][1])) {
      #no format logical because of multiple cases
      config.var$disp.association.large <- tmp[[1]][2]
    }
    
    if(identical("disp.beta.association", tmp[[1]][1])) {
      #no format logical because of reading as string with comma
      config.var$disp.beta.association <- tmp[[1]][2]
    }
    
    if(identical("disp.beta.association.large", tmp[[1]][1])) {
      #no format logical because of multiple cases
      config.var$disp.beta.association.large <- tmp[[1]][2]
    }
    
    if(identical("factor.beta", tmp[[1]][1])) {
      #no format logical because of multiple cases
      config.var$factor.beta <- as.numeric(tmp[[1]][2])
    }
    
    if(identical("disp.region", tmp[[1]][1])) {
      #no format logical because of reading as string with comma
      config.var$disp.region <- tmp[[1]][2]
    }
    
    if(identical("disp.region.large", tmp[[1]][1])) {
      #no logical because of multiple case
      config.var$disp.region.large <- tmp[[1]][2]
    }
    
    if(identical("disp.mydata.names", tmp[[1]][1])) {
      config.var$disp.mydata.names <- as.logical(tmp[[1]][2])
    }
    
    if(identical("disp.pvalueplot", tmp[[1]][1])) {
      config.var$disp.pvalueplot <- as.logical(tmp[[1]][2])
    }
    
    if(identical("fontsize.gviz", tmp[[1]][1])) {
      config.var$fontsize.gviz <- tmp[[1]][2]
    }
    
    if(identical("list.tracks", tmp[[1]][1])) {
      config.var$list.tracks <- tmp[[1]][2]
    }
    
    if(identical("genome", tmp[[1]][1])) {
      config.var$genome <- tmp[[1]][2]
    }
    
    if(identical("dataset.gene", tmp[[1]][1])) {
      config.var$dataset.gene <- tmp[[1]][2]
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
    
    if(identical("pattern.regulation", tmp[[1]][1])) {
      config.var$pattern.regulation <- tmp[[1]][2]
    }
    
    if(identical("BROWSER.SESSION", tmp[[1]][1])) {
      config.var$BROWSER.SESSION <- tmp[[1]][2]
    }
    
    if(identical("tracks.gviz", tmp[[1]][1])) {
      config.var$tracks.gviz <- tmp[[1]][2]
    }
    
    if(identical("tracks.ggbio", tmp[[1]][1])) {
      config.var$tracks.ggbio <- tmp[[1]][2]
    }
    
    if(identical("tracks.trackviewer", tmp[[1]][1])) {
      config.var$tracks.trackviewer <- tmp[[1]][2]
    }
    
    if(identical("palette.file", tmp[[1]][1])) {
      config.var$palette.file <- tmp[[1]][2]
    }
    
    if(identical("disp.type", tmp[[1]][1])) {
      config.var$disp.type <-tmp[[1]][2]
    }
    
    if(identical("disp.cormatrixmap", tmp[[1]][1])) {
      config.var$disp.cormatrixmap <- as.logical(tmp[[1]][2])
    }
    
    if(identical("disp.mydata", tmp[[1]][1])) {
      config.var$disp.mydata <- as.logical(tmp[[1]][2])
    }
    
    if(identical("disp.color.bar", tmp[[1]][1])) {
      config.var$disp.color.bar <- as.logical(tmp[[1]][2])
    }
    
    if(identical("disp.phys.dist", tmp[[1]][1])) {
      config.var$disp.phys.dist <- as.logical(tmp[[1]][2])
    }
    
    if(identical("disp.connecting.lines", tmp[[1]][1])) {
      config.var$disp.connecting.lines <- as.logical(tmp[[1]][2])
    }
    
    if(identical("image.title", tmp[[1]][1])) {
      config.var$image.title <- tmp[[1]][2]
    }
    
    if(identical("disp.legend", tmp[[1]][1])) {
      config.var$disp.legend <- as.logical(tmp[[1]][2])
    }
    
    if(identical("image.type", tmp[[1]][1])) {
      config.var$image.type <- tmp[[1]][2]
    }
    
    if(identical("image.size", tmp[[1]][1])) {
      config.var$image.size <- tmp[[1]][2]
    }
    
    if(identical("disp.mult.lab.X", tmp[[1]][1])) {
      config.var$disp.mult.lab.X <- as.logical(tmp[[1]][2])
    }
    
    if(identical("print.image", tmp[[1]][1])) {
      config.var$print.image <- as.logical(tmp[[1]][2])
    }
    
    if(identical("image.name", tmp[[1]][1])) {
      config.var$image.name <- tmp[[1]][2]
    }
    
    if(identical("connecting.lines.factor", tmp[[1]][1])) {
      config.var$connecting.lines.factor <- as.numeric(tmp[[1]][2])
    }
    
    if(identical("connecting.lines.adj", tmp[[1]][1])) {
      config.var$connecting.lines.adj <- as.numeric(tmp[[1]][2])
    }
    
    if(identical("connecting.lines.vert.adj", tmp[[1]][1])) {
      config.var$connecting.lines.vert.adj <- as.numeric(tmp[[1]][2])
    }
    
    if(identical("connecting.lines.flex", tmp[[1]][1])) {
      config.var$connecting.lines.flex <- as.numeric(tmp[[1]][2])
    }
    
    if(identical("font.factor", tmp[[1]][1])) {
      config.var$font.factor <- as.numeric(tmp[[1]][2])
    }
    
    if(identical("symbol.factor", tmp[[1]][1])) {
      config.var$symbol.factor <- as.numeric(tmp[[1]][2])
    }
    
    if(identical("verbose", tmp[[1]][1])) {
      config.var$verbose <- as.logical(tmp[[1]][2])
    }
  }
  
  #DEBUG STATEMENT
  if (config.var$verbose)  cat("FINISH READ.CONFIG\n")
  
  return(config.var)
}

#------------------Check the list of parameters in CONFiguration file------------------
## for comet and comet.web
check.configVar <- function(config.var) {
  
  #DEBUG STATEMENT
  if (config.var$verbose)  cat("START CHECK.CONFIG\n")
  
  #######################INFO FILE 
  if(is.null(config.var$mydata.file)) {
    stop("No value for info file.\n 
          Parametre: mydata.file")
  }
  
  if(is.null(config.var$mydata.type)) {
    stop("No value for the format of info file.\n 
          Parametre: mydata.type.\n 
          Potential value: FILE, MATRIX\n")
  } else {
    
    if(config.var$mydata.type == "file") {
      split.mydata.file <- strsplit(config.var$mydata.file, ",")
      if( length(split.mydata.file[[1]]) >1 ) {
        
        stop("Need only one file. Put other files via mydata.large.file \n")
      }
    } else {
      split.mydata.file <- list(config.var$mydata.file)
      if( length(split.mydata.file) >1 ) {
        
        stop("Need only one file. Put other files via mydata.large.file \n")
      }
    }
    
  }
  
  if(is.null(config.var$mydata.format)) {
    
    stop("No value for the format of info file.\n 
          Parametre: mydata.format.\n 
          Potential value: site, region, site_asso, region_asso\n")
  }
  
  if(is.null(config.var$symbols)) {
    
    stop("No value for the format of info file.\n 
          Parametre: symbols\n 
          Potential value: circle-fill, square-fill diamont-fill\n")
  }
  
  if(is.null(config.var$sample.labels)) {
    
    stop("No value for the label of info file.\n 
          Parametre: sample.labels\n")
  }
  
  if(is.null(config.var$color.list)) {
    
    if (config.var$verbose)   cat("No value for the color of info file.\n 
          Parametre: color.list\n")
  }
  
  
  if(grepl("asso", config.var$mydata.format)[1] & is.null(config.var$disp.association)) {
    
    stop("No value for the association of info file.\n 
          Parametre: disp.association\n")
  }
  
  if(grepl("asso", config.var$mydata.format)[1] & 
     (is.null(config.var$disp.beta.association) | is.null(config.var$factor.beta)))  {
    
    stop("No value for the association of info file.\n 
          Parametre: disp.beta.association or factor.beta\n")
  }
  
  if(grepl("region", config.var$mydata.format)[1] & is.null(config.var$disp.region)) {
    
    stop("No value for the region of info file.\n 
          Parametre: disp.region\n")
  }
  
  #---------------- Supplementary file
  if(is.null(config.var$mydata.large.file)) {
    if (config.var$verbose) cat("No value for the supplementary file.\n 
          Parametre: mydata.large.file ")
  }
  
  if( ! is.null(config.var$mydata.large.file)) {
    if(is.null(config.var$mydata.large.format)) {
      
      stop("No value for the format of supplementary data.\n 
          Parametre: mydata.large.format.\n 
          Potential value: site, region, site_asso, region_asso\n")
    }
    if(is.null(config.var$mydata.large.type)) {
      
      stop("No value for the format of supplementary data.\n 
          Parametre: mydata.large.type.\n 
          Potential value: LISTFILE or listdataframe \n")
    }
    
    if(is.null(config.var$symbols.large)) {
      
      stop("No value of symbole for supplementary data\n
           Parameter: symbols.large\n")
    }
    
    if(is.null( config.var$sample.labels.large)) {
      
      stop("No value of label for supplementary data\n
           Parameter: sample.labels.large\n")
    }
    
    if(is.null(config.var$color.list.large)) {
      
      stop("No value of color for supplementary data\n
           Parameter: color.list.large\n")
    }
    
    if(grepl("asso", config.var$mydata.large.format)[1] & 
       is.null(config.var$disp.association.large)) {
      
      stop("No value of association for supplementary data\n
           Parameter: disp.association.large\n")
    }
    
    if(grepl("asso", config.var$mydata.large.format)[1] & 
       (is.null(config.var$disp.beta.association.large) | is.null(config.var$factor.beta)))  {
      
      stop("No value for the association of info file.\n 
          Parametre: disp.beta.association.large or factor.beta\n")
    }
    
    if(grepl("region", config.var$mydata.large.format)[1] &
       is.null(config.var$disp.region.large)) {
      
      stop("No value of region for supplementary data\n
           Parameter: disp.region.large\n")
    }
    
    if(config.var$mydata.large.format == "file") {
      split.mydata.large.file <- strsplit(config.var$mydata.large.file, ",")
      length.datalarge <- length(split.mydata.large.file[[1]])
    } else {
      split.mydata.large.file <- config.var$mydata.large.file
      length.datalarge <- length(split.mydata.large.file)
    }
    
    split.mydata.large.format <- strsplit(config.var$mydata.large.format, ",")
    split.mydata.large.region <- strsplit(config.var$disp.region.large, ",")
    split.mydata.large.asso <- strsplit(config.var$disp.association.large, ",")
    split.mydata.large.beta.asso <- strsplit(config.var$disp.beta.association.large, ",")
    split.mydata.large.color <- strsplit(config.var$color.list.large, ",")
    split.mydata.large.label <- strsplit(config.var$sample.labels.large, ",")
    split.mydata.large.symbol <- strsplit(config.var$symbols.large, ",")
    
    if( length.datalarge !=  length(split.mydata.large.format[[1]]) ||
        length.datalarge !=  length(split.mydata.large.region[[1]]) ||
        length.datalarge !=  length(split.mydata.large.asso[[1]]) ||
        length.datalarge !=  length(split.mydata.large.beta.asso[[1]]) ||
        length.datalarge !=  length(split.mydata.large.color[[1]]) ||
        length.datalarge !=  length(split.mydata.large.label[[1]]) ||
        length.datalarge !=  length(split.mydata.large.symbol[[1]])) {
      
      stop("Need to have the same number of element for all options related to supplementary data\n
           mydata.large.file, mydata.large.format, disp.region.large, disp.association.large, 
           disp.beta.association.large, color.list.large,
           sample.labels.large, symbols.large: ",length.datalarge,"\n")
    }
  } 
  
  #-----------------------Correlation matrix
  if(is.null(config.var$cormatrix.file)) {
    
    stop("No value for the correlation file.\n 
          Parametre: cormatrix.file ")
  } else {
    if(is.null(config.var$disp.mydata)) {
      
      stop("No value to show the correlation file.\n 
          Parametre: disp.mydata ")
    }
    
    if(is.null(config.var$cormatrix.type)) {
      
      stop("No value for the correlation file.\n 
          Parametre: cormatrix.type \n
          Potential value: listfile or listdataframe")
    }
    
    if(is.null(config.var$cormatrix.conf.level)) {
      
      stop("No value for the confiance level in the visualisation of correlation matrix.\n 
          Parametre: cormatrix.conf.level \n
          Potential value:  under or equal 1")
    }
    
    if(is.null(config.var$cormatrix.sig.level)) {
      
      stop("No value for the pvalue in the visualisation of correlation matrix.\n 
          Parametre: cormatrix.sig.level \n
          Potential value:  under or equal 1")
    }
    
    if(is.null(config.var$cormatrix.adjust)) {
      
      stop("No value for the adjust hypothesis.\n 
          Parametre: cormatrix.adjust \n
          Potential value: holm, hochberg, hommel, bonferroni, BH, BY, fdr, none ")
    }
    
    if(is.null(config.var$cormatrix.format)) {
      
      stop("No value for the correlation file.\n 
          Parametre: cormatrix.format \n
          Potential value: CORMATRIX, RAW")
    }
    
    if(config.var$cormatrix.format == "raw" & 
       is.null(config.var$cormatrix.method)) {
      
      stop("No value for the method of correlation file for raw.\n 
          Parametre: cormatrix.method \n
          Potential value: spearman, pearson, kendall")
    }
    
    if(is.null(config.var$cormatrix.color.scheme)) {
      
      stop("No value for the method of correlation file for raw.\n 
          Parametre: cormatrix.color.scheme \n
          Potential value: spearman, pearson, kendall")
    }
    
    if(is.null(config.var$disp.cormatrixmap)) {
      
      stop("Show heatmap of correlation matrix.\n 
          Parametre: disp.cormatrixmap \n
          Potential value: TRUE or FALSE")
    }
  }
  
  #-----------------------Parameters
  if(is.null(config.var$mydata.ref)) {
    if (config.var$verbose)   cat("No value for the reference.\n
              Parameter: mydata.ref\n")
  }
  
  if(is.null(config.var$start)) {
    if (config.var$verbose)   cat("No value for the beginning of genomic region\n
              Parameter: START\n")
  }
  
  if(is.null(config.var$end)) {
    if (config.var$verbose)   cat("No value for the end of genomic region\n
            Parameter: END\n")
  }
  
  if(is.null(config.var$lab.Y)) {
    
    stop("No value of axis\n
           Parameter: lab.Y\n
           Potential value: log10, ln \n")
  }
  
  if(is.null(config.var$pval.threshold)) {
    
    stop("No value of threshold of the significance\n
           Parameter: pval.threshold\n
           Potential value: 10e-8\n")
  }
  
  if(is.null(config.var$pval.threshold.2)) {
    
    stop("No value of threshold 2 of the significance\n
           Parameter: pval.threshold. 2\n
           Potential value: 10e-5\n")
  }
  
  if(is.null(config.var$disp.pval.threshold)) {
    
    stop("No value of threshold to visualise the pvalue\n
           Parameter: disp.pval.threshold\n
           Potential value: 0\n")
  }
  
  if(is.null(config.var$disp.color.ref)) {
    
    stop("No value of threshold to visualise the pvalue\n
           Parameter: disp.color.ref\n
           Potential value: TRUE or FALSE\n")
  }
  
  if(is.null(config.var$use.colors)) {
    
    stop("No value to use\n
           Parameter: use.colors\n
           Potential value: TRUE or FALSE\n")
  }
  
  if(is.null(config.var$disp.mydata.names)) {
    
    stop("No value to show the name of genes\n
           Parameter: disp.mydata.names\n
           Potential value: TRUE or FALSE\n")
  }
  
  if(is.null(config.var$disp.pvalueplot)) {
    
    stop("No value to show the pvalue plot \n
           Parameter: disp.pvalueplot\n
           Potential value: TRUE or FALSE\n")
  }
  
  if(is.null(config.var$genome)) {
    
    stop("No value of genome\n
           Parameter: genome\n
           Potential value: hg19\n")
  }
  
  
  if(is.null(config.var$BROWSER.SESSION)) {
    
    stop("No value of name of browser session\n
           Parameter: BROWSER.SESSION\n
           Potential value: UCSC\n")
  }
  
  if(is.null(config.var$list.tracks)) {
    if (config.var$verbose)   cat("No value of annotation tracks\n
           Parameter: list.tracks\n
           Potential value: geneENSEMBL,CGI,ChromHMM,DNAse,RegENSEMBL,SNP\n")
  }else {
    
    if(is.null(config.var$tracks.gviz)) {
      if (config.var$verbose)   cat("No list of Gviz track\n
           Parameter: tracks.gviz\n")
    }
    
    if(is.null(config.var$tracks.ggbio)) {
      if (config.var$verbose)   cat("No list of GGbio track\n
           Parameter: tracks.ggbio\n")
    }
    
    if(is.null(config.var$tracks.trackviewer)) {
      if (config.var$verbose)   cat("No list of TrackViewer track\n
           Parameter: tracks.trackviewer\n")
    }
  }
  
  if(is.null(config.var$dataset.gene)) {
    if (config.var$verbose)   cat("No value of dataset for ENSEMBL\n
           Parameter: dataset.gene\n
          Potential value: hsapiens_gene\n")
  }
  
  if(is.null(config.var$DATASET.SNP)) {
    if (config.var$verbose)   cat("No value of dataset for ENSEMBL\n
          Parameter: DATASET.SNP\n
          Potential value: hsapiens_snp\n")
  }else {
    if(is.null(config.var$VERSION.DBSNP)) {
      if (config.var$verbose)   cat("No value of version of DBsnp\n
            Parameter: VERSION.DBSNP\n
            Potential value: snp138\n")
    }
  }
  
  if(is.null(config.var$DATASET.SNP.STOMA)) {
    if (config.var$verbose)   cat("No value of dataset for ENSEMBL\n
          Parameter: DATASET.SNP.STOMA\n
          Potential value: hsapiens_snp_som\n")
  }
  
  if(is.null(config.var$DATASET.REGULATION)) {
    if (config.var$verbose)   cat("No value of dataset for ENSEMBL\n
          Parameter: DATASET.REGULATION\n
          Potential value: hsapiens_feature_set\n")
  }
  
  if(is.null(config.var$DATASET.STRU)) {
    if (config.var$verbose)   cat("No value of dataset for ENSEMBL\n
          Parameter: DATASET.STRU\n
          Potential value: hsapiens_structvar\n")
  }
  
  if(is.null(config.var$DATASET.STRU.STOMA)) {
    if (config.var$verbose)   cat("No value of dataset for ENSEMBL\n
          Parameter: DATASET.STRU\n
          Potential value: hsapiens_structvar_som\n")
  }
  
  if(is.null(config.var$pattern.regulation)) {
    if (config.var$verbose)   cat("No value of name of tissue or pattern of DBregulation for ENSEMBL\n
          Parameter: pattern.regulation\n
          Potential value: GM12878\n")
  }
  
  if(is.null(config.var$palette.file)) {
    if (config.var$verbose)   cat("No list of TrackViewer track\n
           Parameter: palette.file\n")
  }
  
  if(is.null(config.var$biofeat.user.file)) {
    if (config.var$verbose)   cat("No user file to annotation\n
           Parameter: biofeat.user.file\n")
  } else {
    if(is.null(config.var$biofeat.user.type)) {
      
      stop("No type of visualise of user file\n
           Parameter: biofeat.user.type\n")
    }
    
    if(is.null(config.var$biofeat.user.type.plot)) {
      
      stop("No type of visualise of user file\n
           Parameter: biofeat.user.type.plot\n")
    }
    
    split.biouser.file <- strsplit(config.var$biofeat.user.file, ",")
    split.biouser.type <- strsplit(config.var$biofeat.user.type, ",")
    split.biouser.plot <- strsplit(config.var$biofeat.user.type.plot, ",")
    
    
    if( length(split.biouser.file[[1]]) !=  length(split.biouser.type[[1]]) ||
        length(split.biouser.file[[1]]) !=  length(split.biouser.plot[[1]])  ) {
      
      stop("Need to have the same number of element for all options related to supplementary data\n
           biofeat.user.type.plot, biofeat.user.file, biofeat.user.type, biofeat.user.type.plot\n")
    }
  }
  
  if(is.null(config.var$disp.color.bar)) {
    
    stop("No value to show the color bar\n
           Parameter: disp.color.bar\n
       Potential value; TRUE or FALSE\n")
  }
  
  if(is.null(config.var$disp.phys.dist)) {
    
    stop("No value to show the physique distance\n
           Parameter: disp.phys.dist\n
       Potential value; TRUE or FALSE\n")
  }
  
  if(is.null(config.var$disp.connecting.lines)) {
    
    stop("No value to show the connecting lines\n
           Parameter: disp.connecting.lines\n
       Potential value; TRUE or FALSE\n")
  }else {
    if(is.null(config.var$connecting.lines.factor)) {
      
      stop("No value for the connecting lines factor \n
           Parameter: connecting.lines.factor\n")
    }
    
    if(is.null(config.var$connecting.lines.adj)) {
      
      stop("No value for the connecting lines adj \n
           Parameter: connecting.lines.adj\n")
    }
    
    if(is.null(config.var$connecting.lines.vert.adj)) {
      
      stop("No value for the connecting lines vertical adj \n
           Parameter: connecting.lines.vert.adj\n")
    }
    
    if(is.null(config.var$connecting.lines.flex)) {
      
      stop("No value for the connecting lines flexible \n
           Parameter: connecting.lines.flex\n")
    }
  }
  
  if(is.null(config.var$disp.legend)) {
    
    stop("No value to show the legend\n
           Parameter: disp.legend\n
       Potential value; TRUE or FALSE\n")
  }
  
  if(is.null(config.var$disp.mult.lab.X)) {
    if (config.var$verbose)   cat("No value to show the multiple lable\n
           Parameter: disp.mult.lab.X\n
       Potential value; TRUE or FALSE\n")
  }
  
  if(is.null(config.var$font.factor)) {
    if (config.var$verbose)   cat("No value for the font factor\n
           Parameter: font.factor\n")
  }
  
  if(is.null(config.var$symbol.factor)) {
    if (config.var$verbose)   cat("No value for the symbol factor\n
           Parameter: symbol.factor\n")
  }
  
  if(is.null(config.var$image.title)) {
    
    stop("No title for the plot\n
           Parameter: image.title\n")
  }
  
  if(is.null(config.var$image.type)) {
    
    stop("No value for the type of image\n
           Parameter: image.type\n
       Potential value; PDF or EPS\n")
  }
  
  if(is.null(config.var$image.size)) {
    
    stop("No value for the size of image\n
           Parameter: image.size\n
       Potential value; 3.5 or 7\n")
  }
  
  if(is.null(config.var$print.image)) {
    
    stop("No value for the size of image\n
           Parameter: image.size\n
       Potential value; 3.5 or 7\n")
  }
  
  if(is.null(config.var$image.name)) {
    
    stop("No value for the name of image\n
           Parameter: image.name\n")
  }
  
  #  if(is.null(config.var$zoom)) {
  #    if (config.var$verbose)   cat("No value for the end of genomic region\n")
  #  }
  #  if(is.null(config.var$disp.type)) {
  #    if (config.var$verbose)   cat("No value for the end of genomic region\n")
  #  }
  
  #DEBUG STATEMENT
  if (config.var$verbose)  cat("FINISH CHECK.CONFIG\n")
  
  return(config.var)
}

## for comet.list
check.configVar.cometlist <- function(config.var) {
  
  #DEBUG STATEMENT
  if (config.var$verbose)  cat("START CHECK.CONFIG.COMETLIST\n")
  
  #-----------------------Correlation matrix
  if(is.null(config.var$cormatrix.file)) {
    
    stop("No value for the correlation file.\n 
         Parametre: cormatrix.file ")
  } else {
    
    if(is.null(config.var$cormatrix.type)) {
      
      stop("No value for the correlation file.\n 
           Parametre: cormatrix.type \n
           Potential value: listfile or listdataframe")
    }
    
    if(is.null(config.var$cormatrix.conf.level)) {
      
      stop("No value for the confiance level in the visualisation of correlation matrix.\n 
           Parametre: cormatrix.conf.level \n
           Potential value:  under or equal 1")
    }
    
    if(is.null(config.var$cormatrix.sig.level)) {
      
      stop("No value for the pvalue in the visualisation of correlation matrix.\n 
           Parametre: cormatrix.sig.level \n
           Potential value:  under or equal 1")
    }
    
    if(is.null(config.var$cormatrix.adjust)) {
      
      stop("No value for the adjust hypothesis.\n 
           Parametre: cormatrix.adjust \n
           Potential value: holm, hochberg, hommel, bonferroni, BH, BY, fdr, none")
    }
    
    if(is.null(config.var$cormatrix.output)) {
      
      stop("No name for the output file.\n 
           Parametre: cormatrix.output \n
           No need the extention, need absolue path")
    }
    
    if(is.null(config.var$cormatrix.format)) {
      
      stop("No value for the correlation file.\n 
           Parametre: cormatrix.format \n
           Potential value: CORMATRIX, RAW")
    }
    
    if(config.var$cormatrix.format == "raw" & 
       is.null(config.var$cormatrix.method)) {
      
      stop("No value for the method of correlation file for raw.\n 
           Parametre: cormatrix.method \n
           Potential value: spearman, pearson, kendall")
    }
  }
  
  #-----------------------Parameters
  if(is.null(config.var$mydata.ref)) {
    if (config.var$verbose)   cat("No value for the reference.\n
                                  Parameter: mydata.ref\n")
  }
  
  #DEBUG STATEMENT
  if (config.var$verbose)  cat("FINISH CHECK.CONFIG.COMETLIST\n")
  
  return(config.var)
}

#------------------UPDATE VARIABLE------------------
## For comet and comet.web
retrieve.data <- function(config.var, gbl.var) {
  
  #DEBUG STATEMENT
  if (config.var$verbose)  cat("START RETRIEVE.DATA\n")
  
  #initialize distance variables
  min.dist <- NULL
  max.dist <- NULL
  all.beta <- NULL
  
  gbl.var$mydata.gen <- config.var$genome
  split.mydata.file <- gbl.var$split.mydata.file
  split.cormatrix.file <- NULL
  split.cormatrix.file <- gbl.var$split.cormatrix.file
  
  numdata <- 0
  if(gbl.var$split.type == "file"){
    numdata <- length(split.mydata.file[[1]])
  } else {
    numdata <- length(split.mydata.file)
  }
  #-------------------- READ FILES MYDATA
  for(i in 1:numdata) {
    gbl.var$cur.sample <- i
    #-------------------- READ FILE CpG DATA
    gbl.var <-read.file.mydata(split.mydata.file,config.var, gbl.var,i)
    if (config.var$verbose)  cat("READ ",i," MYDATA file \n")
    
    cur.mydata.num <- nrow(gbl.var$mydata.data)
    #gbl.var$palette.size <- cur.mydata.num
    gbl.var$palette.size <- 100
    #create 2 hash tables containing MYDATA names and positions for current and previous sample sets
    #one hash has the name as the key and the other the location because R has no way of looking
    #up hash keys based on values
    format <- as.character(gbl.var$split.format[[1]][gbl.var$cur.sample])
    if (config.var$verbose)  cat("format",format,"\n")
    if(grepl("region", format)[1] ) {
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
    ## if (config.var$verbose)  cat("length(mydata.hash.names) ", length(mydata.hash.names), "\n")
    ## if (config.var$verbose)  cat("length(mydata.hash.pos) ", length(mydata.hash.pos), "\n")
    
    #create a variable with only the sorted MYDATA positions
    gbl.var$sorted.mydata.pos <- sort(as.numeric(mydata.hash.pos))
    
    #sort the MYDATA names
    for(n in 1:length(mydata.hash.pos)) {
      #  if (config.var$verbose)  cat("POS: ", mydata.hash.pos[i], "NAME: ", get(mydata.hash.pos[i], envir=mydata.hash.pos.name)$mydata.name, "\n")
      gbl.var$sorted.mydata.names[n] <- get(mydata.hash.pos[n], envir=gbl.var$mydata.hash.list.pos.names)$mydata.name
    }
    
    #get the number of unique MYDATA in all current sample sets
    gbl.var$mydata.num <- length(gbl.var$sorted.mydata.names)
    
    #DEBUG STATEMENT
    if (config.var$verbose)  cat("gbl.var$mydata.num ", gbl.var$mydata.num, "\n")
    
    #------------- Define Format visualisation CpG data
    
    #Set the gbl.var$cex.factor, which determines font size
    formatting.var.list <- set.image.parameters(config.var, gbl.var)
    
    gbl.var$font.size <- formatting.var.list$font.size
    gbl.var$line.width <- formatting.var.list$line.width
    gbl.var$cex.factor <- formatting.var.list$cex.factor
    gbl.var$cex.factor.symbol <- formatting.var.list$cex.factor.symbol
    
    #DEBUG STATEMENT
    # if (config.var$verbose)  cat("gbl.var$cex.factor ", gbl.var$cex.factor, "\n")
    
    if(config.var$image.size == 3.5) {
      gbl.var$mydata.names <- substr(gbl.var$mydata.data$MYDATA.NAME, 1, 9)
    } else if(config.var$image.size == 7) {
      gbl.var$mydata.names <- substr(gbl.var$mydata.data$MYDATA.NAME, 1, 9)
    } else {
      
      stop("Invalid image size: ", config.var$image.size, "\n")
    }
    
    #----------------- DEFINE VALUES OF X,Y and DISTANCE to position data
    if(!is.null(config.var$start)) {
      min.dist <- config.var$start
    } else {
      if(is.null(min.dist)){
        if(is.null(gbl.var$mydata.data$LOC)) {
          
          stop("No value for LOC (min.dist)",gbl.var$mydata.data$LOC)
        }
        min.dist <- min(gbl.var$mydata.data$LOC, na.rm=TRUE)
      }else {
        min.dist <- min(c(min.dist, gbl.var$mydata.data$LOC), na.rm = TRUE)
      }
      min.dist <- min.dist - 500
    }
    
    
    gbl.var$min.user.x <- min.dist
    
    if(!is.null(config.var$end)){
      max.dist <- config.var$end 
      #max.dist <- config.var$end
    } else {
      
      if(is.null(max.dist)){
        if(is.null(gbl.var$mydata.data$LOC)){
          
          stop("NO value for LOC (max.dist)", gbl.var$mydata.data$LOC)
        }
        max.dist <- max(gbl.var$mydata.data$LOC, na.rm = TRUE)
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
    
    if(!is.null(config.var$start)) {
      gbl.var$min.x <- config.var$start 
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
    #  if (config.var$verbose)  cat("gbl.var$min.x", gbl.var$min.x, "\n")
    
    
    if(!is.null(config.var$end)) {
      gbl.var$max.x <- config.var$end 
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
    #   if (config.var$verbose)  cat("gbl.var$max.x", gbl.var$max.x, "\n")
    
    #------------------- DEFINE DATA FOR ONLY zoom REGION -----------------------
    ## if (config.var$verbose)  cat("Sort data", gbl.var$sorted.mydata.pos, "\n")
    ## if (config.var$verbose)  cat("Sort  length data", is.vector(gbl.var$sorted.mydata.pos), "\n")
    remove.sorted.mydata.pos <-which(gbl.var$sorted.mydata.pos < gbl.var$min.x | gbl.var$sorted.mydata.pos > gbl.var$max.x)
    # if (config.var$verbose)  cat("Sort oouutt", length(remove.sorted.mydata.pos), "\n")
    if(length(remove.sorted.mydata.pos) > 0){
      sorted.mydata.pos <-  gbl.var$sorted.mydata.pos
      #  if (config.var$verbose)  cat("List gb  ",sorted.mydata.pos,"\n")
      sorted.mydata.pos.zoom <- sorted.mydata.pos[-remove.sorted.mydata.pos]
      #  if (config.var$verbose)  cat("List gb  ",sorted.mydata.pos.zoom,"\n")
    }else {
      sorted.mydata.pos.zoom <- gbl.var$sorted.mydata.pos
    }
    gbl.var$sorted.mydata.pos.zoom <- sorted.mydata.pos.zoom  
    # if (config.var$verbose)  cat("gbl.var$sorted.mydata.pos.zoom ", gbl.var$sorted.mydata.pos.zoom, "\n")
    
    
    #----------------- DEFINE THE HIGHEST VALUE OF Y
    gbl.var$min.y <- min(c(gbl.var$min.y, min(gbl.var$mydata.data$MYDATA.PVAL, na.rm = TRUE)))
    # Ensure that the y scale extends beyond the data region, round to a whole number
    gbl.var$max.y <- ceiling(max(c(gbl.var$max.y, max(gbl.var$mydata.data$MYDATA.PVAL, na.rm = TRUE))))
    
    #DEBUG STATEMENT
    # if (config.var$verbose)  cat("gbl.var$min.y ", gbl.var$min.y, "\n")
    #  if (config.var$verbose)  cat("gbl.var$max.y ", gbl.var$max.y, "\n")
    
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
    # if (config.var$verbose)  cat("config.var$disp.legend ", config.var$disp.legend, "\n")
    # if (config.var$verbose)  cat("gbl.var$max.y ", gbl.var$max.y, "\n")
    
    if(config.var$disp.legend & !config.var$disp.cormatrixmap) {
      
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
      if(is.null(config.var$mydata.ref)){
        best.position <- which.max(gbl.var$mydata.data$MYDATA.PVAL)
      }else {
        best.position <- which(gbl.var$mydata.data$MYDATA.NAME == config.var$mydata.ref)
      }
      if(grepl("region", format)[1] ) {
        gbl.var$mydata.ref.pos <- gbl.var$mydata.data$LOC.START[best.position] + ((gbl.var$mydata.data$LOC.END[best.position] - gbl.var$mydata.data$LOC.START[best.position])/2)
        gbl.var$mydata.ref.pos <- gbl.var$mydata.data$LOC.START[best.position] + ((gbl.var$mydata.data$LOC.END[best.position] - gbl.var$mydata.data$LOC.START[best.position])/2)
      } else {
        gbl.var$mydata.ref.pos <- gbl.var$mydata.data$LOC[best.position]
        gbl.var$mydata.ref.pos <- gbl.var$mydata.data$LOC[best.position]
      }
      gbl.var$mydata.best.position <- best.position
      if (config.var$verbose)  cat("Ref",gbl.var$mydata.ref.pos,"\n")
      
      #---------- EXTRACT BETA VALUES
      if(grepl("asso", format)[1] & is.integer(gbl.var$mydata.data$MYDATA.ASSO.ORIGINAL[1])) {
        all.beta <- c(all.beta,gbl.var$mydata.data$MYDATA.ASSO.ORIGINAL)
      }
      
      #-----------------READ MATRICE for CORMATRIX MAP
      # if (config.var$verbose)  cat("config.var$disp.cormatrixmap ", config.var$disp.cormatrixmap, "\n")
      #Need to read even if we don't show data
      gbl.var <-read.file.cormatrix(config.var, gbl.var,split.cormatrix.file)
    }
    
    
  }
  
  
  #------------------- READ LARGE DATA
  if(!is.null(config.var$mydata.large.file)){
    large.split.mydata.file=gbl.var$large.split.mydata.file
    
    numlargedata <- 0
    if(gbl.var$large.split.type == "listfile"){
      numlargedata <- length(large.split.mydata.file[[1]])
    } else {
      numlargedata <- length(large.split.mydata.file)
    }
    for(i in 1:numlargedata) {
      #-------------------- READ FILE CpG DATA
      gbl.var <-read.file.mydata.large(large.split.mydata.file,config.var, gbl.var,i)
      
      cur.mydata.large.num <- nrow(gbl.var$mydata.large.data)
      #create 1 hash table containing MYDATA names and positions for current and previous sample sets
      #one hash has the name as the key and the other the location because R has no way of looking
      #up hash keys based on values
      format <- as.character(gbl.var$large.split.format[[1]][i])
      
      #---------- EXTRACT BETA VALUES
      if(grepl("asso", format)[1] & is.integer(gbl.var$mydata.large.data$MYDATA.ASSO.ORIGINAL[1])) {
        all.beta <- c(all.beta,gbl.var$mydata.large.data$MYDATA.ASSO.ORIGINAL)
      }
      
      if(grepl("region", format)[1] ) {
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
        #  if (config.var$verbose)  cat("test hash",values(gbl.var$mydata.large.hash.names.pos),"\n")
        
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
      ## if (config.var$verbose)  cat("length(mydata.hash.pos) ", length(mydata.hash.pos), "\n")
      
      #create a variable with only the sorted MYDATA positions
      gbl.var$sorted.mydata.large.pos <- sort(as.numeric(mydata.large.hash.pos))
      
      #------------------- zoom on LARGE DATA -----------------------
      ## if (config.var$verbose)  cat("Sort data", gbl.var$sorted.mydata.large.pos, "\n")
      ## if (config.var$verbose)  cat("Sort  length data", is.vector(gbl.var$sorted.mydata.large.pos), "\n")
      remove.sorted.mydata.large.pos <-which(gbl.var$sorted.mydata.large.pos < gbl.var$min.x | gbl.var$sorted.mydata.large.pos > gbl.var$max.x)
      # if (config.var$verbose)  cat("Sort oouutt", length(remove.sorted.mydata.large.pos), "\n")
      if(length(remove.sorted.mydata.large.pos) > 0){
        sorted.mydata.large.pos <-  gbl.var$sorted.mydata.large.pos
        #  if (config.var$verbose)  cat("List gb  ",sorted.mydata.large.pos,"\n")
        sorted.mydata.large.pos.zoom <- sorted.mydata.large.pos[-remove.sorted.mydata.large.pos]
        #  if (config.var$verbose)  cat("List gb  ",sorted.mydata.large.pos.zoom,"\n")
      }else {
        sorted.mydata.large.pos.zoom <- gbl.var$sorted.mydata.large.pos
      }
      gbl.var$sorted.mydata.large.pos.zoom <- sorted.mydata.large.pos.zoom  
      # if (config.var$verbose)  cat("gbl.var$sorted.mydata.large.pos.zoom ", gbl.var$sorted.mydata.large.pos.zoom, "\n")
      
      
      
    }
    
    #----------------- DEFINE THE HIGHEST VALUE OF Y
    #DEBUG STATEMENT
    #  if (config.var$verbose)  cat("gbl.var$min.y ", gbl.var$min.y, "\n")
    # if (config.var$verbose)  cat("gbl.var$max.y ", gbl.var$max.y, "\n")
    if(config.var$lab.Y == "ln") {
      gbl.var$mydata.large.data$MYDATA.PVAL <- -log(gbl.var$mydata.large.data$MYDATA.PVAL) # LOG Variable
    }
    if(config.var$lab.Y == "log") {
      gbl.var$mydata.large.data$MYDATA.PVAL <- -log10(gbl.var$mydata.large.data$MYDATA.PVAL) # LOG Variable
    }
    
    gbl.var$min.y <- min(c(gbl.var$min.y, min(gbl.var$mydata.large.data$MYDATA.PVAL, na.rm = TRUE)))
    # Ensure that the y scale extends beyond the data region, round to a whole number
    gbl.var$max.y <- ceiling(max(c(gbl.var$max.y, max(gbl.var$mydata.large.data$MYDATA.PVAL, na.rm = TRUE))))
    
    #DEBUG STATEMENT
    # if (config.var$verbose)  cat("gbl.var$min.y ", gbl.var$min.y, "\n")
    # if (config.var$verbose)  cat("gbl.var$max.y ", gbl.var$max.y, "\n")
    if(config.var$disp.legend & !config.var$disp.cormatrixmap) {
      
      gbl.var$axis.y <- seq(as.integer(gbl.var$min.y), gbl.var$max.y, 1) # BY Variable
      gbl.var$r.y <- c(gbl.var$min.y, gbl.var$max.y + 0.5) # 0 - max(PVAL) for MYDATA
    }
    else {
      gbl.var$axis.y <- seq(as.integer(gbl.var$min.y), gbl.var$max.y, 1) # BY Variable
      gbl.var$r.y <- c(gbl.var$min.y, gbl.var$max.y) # 0 - max(PVAL) for MYDATA
    }
    
  }
  
  #------ DEFINE MEAN and SD of ALL BETA
  if (config.var$verbose)  cat("Longer of all.beta",length(all.beta),"\n")
  if(is.null(all.beta)){
    gbl.var$mean.beta <- 0
    gbl.var$sd.beta <- 0
  } else {
    gbl.var$mean.beta <- mean(all.beta)
    gbl.var$sd.beta <- sd(all.beta)
  }
  
  if (config.var$verbose)   cat("Mean beta",gbl.var$mean.beta,"\n")
  if (config.var$verbose)   cat("SD beta",gbl.var$sd.beta,"\n")
  
  
  #----- CREATE LIST of BIOFEATURES -------------------------
  
  if(!is.null(config.var$biofeat.user.file)) {
    gbl.var <- create.tracks.user(config.var, gbl.var)
  }
  
  
  #----- FIX VALUES -------------------------
  
  fix.var.generic <- fix.values.generic(config.var, gbl.var)
  
  config.var <- fix.var.generic$config.var
  gbl.var <- fix.var.generic$gbl.var
  
  
  #DEBUG STATEMENT
  if (config.var$verbose)  cat("FINISH RETRIEVE.DATA\n")
  
  return(list(config.var = config.var, gbl.var = gbl.var))
}

## For comet.list
retrieve.data.cometlist <- function(config.var, gbl.var) {
  #DEBUG STATEMENT
  if (config.var$verbose)  cat("START RETRIEVE.DATA.COMETLIST\n")
  
  #initialize distance variables
  split.cormatrix.file <- NULL
  split.cormatrix.file <- gbl.var$split.cormatrix.file
  
  #-----------------READ MATRICE for CORMATRIX MAP
  # if (config.var$verbose)  cat("config.var$disp.cormatrixmap ", config.var$disp.cormatrixmap, "\n")
  #Need to read even if we don't show data
  gbl.var <-read.file.cormatrix(config.var, gbl.var,split.cormatrix.file)
  
  
  #DEBUG STATEMENT
  if (config.var$verbose)  cat("FINISH RETRIEVE.DATA.COMET.LIST\n")
  
  return(list(config.var = config.var, gbl.var = gbl.var))
}

#-------------------FIX VALUES GENERIC----------------------------
fix.values.generic <- function(config.var, gbl.var) {
  
  #DEBUG STATEMENT
  if (config.var$verbose)  cat("START FIX.VALUES.GENERIC\n")
  # if (config.var$verbose)  cat("PVALUE",gbl.var$mydata.data[1,4],"\n")
  
  if(config.var$lab.Y == "ln") {
    
    config.var$pval.threshold <- -log(config.var$pval.threshold)
    config.var$pval.threshold.2 <- -log(config.var$pval.threshold.2)
    
    #DEBUG STATEMENT
    # if (config.var$verbose)  cat("pval.threshold ", config.var$disp.pval.threshold, "\n")
    # if (config.var$verbose)  cat("pval.threshold ", config.var$pval.threshold, "\n")
    
    if(gbl.var$pval.flag) {
      if(config.var$disp.pval.threshold != 1) {
        config.var$disp.pval.threshold <- -log(config.var$disp.pval.threshold)
        gbl.var$pval.flag <- FALSE
      } else {
        config.var$disp.pval.threshold <- -1
        gbl.var$pval.flag <- FALSE
      }
    }
  } else if (config.var$lab.Y == "log") {
    
    config.var$pval.threshold <- -log10(config.var$pval.threshold)
    config.var$pval.threshold.2 <- -log10(config.var$pval.threshold.2)
    
    #DEBUG STATEMENT
    # if (config.var$verbose)  cat("disp.pval.threshold ", config.var$disp.pval.threshold, "\n")
    # if (config.var$verbose)  cat("pval.threshold ", config.var$pval.threshold, "\n")
    
    if(gbl.var$pval.flag) {
      if(config.var$disp.pval.threshold != 1) {
        config.var$disp.pval.threshold <- -log10(config.var$disp.pval.threshold)
        gbl.var$pval.flag <- FALSE
      } else {
        config.var$disp.pval.threshold <- -1
        gbl.var$pval.flag <- FALSE
      }
    }
  } else {
    
    stop("Invalid lab.Y: ", config.var$lab.Y, ". Specify either 'log' or 'ln'.\n")
  }
  
  #DEBUG STATEMENT
  if (config.var$verbose)  cat("FINISH FIX.VALUES.GENERIC\n")
  
  return(list(config.var = config.var, gbl.var = gbl.var))
}

#-------------------FIX VALUES for MYDATA----------------------------
fix.values <- function(config.var, gbl.var) {
  
  #DEBUG STATEMENT
  if (config.var$verbose)  cat("START FIX.VALUES\n")
  # if (config.var$verbose)  cat("PVALUE",gbl.var$mydata.data[1,4],"\n")
  
  if(config.var$lab.Y == "ln") {
    
    gbl.var$mydata.data$MYDATA.PVAL <- -log(gbl.var$mydata.data$MYDATA.PVAL) # LOG Variable
    
  } else if (config.var$lab.Y == "log") {
    
    gbl.var$mydata.data$MYDATA.PVAL <- -log10(gbl.var$mydata.data$MYDATA.PVAL) # LOG Variable
    
  } else {
    
    stop("Invalid lab.Y: ", config.var$lab.Y, ". Specify either 'log' or 'ln'.\n")
  }
  
  #DEBUG STATEMENT
  if (config.var$verbose)  cat("FINISH FIX.VALUES\n")
  
  return(list(config.var = config.var, gbl.var = gbl.var))
}

#-------------------FIX VALUES for MYDATA----------------------------
fix.values.large <- function(config.var, gbl.var) {
  
  #DEBUG STATEMENT
  if (config.var$verbose)  cat("START FIX.VALUES\n")
  # if (config.var$verbose)  cat("PVALUE",gbl.var$mydata.data[1,4],"\n")
  
  if(config.var$lab.Y == "ln") {
    
    gbl.var$mydata.large.data$MYDATA.PVAL <- -log(gbl.var$mydata.large.data$MYDATA.PVAL) # LOG Variable
    
  } else if (config.var$lab.Y == "log") {
    
    gbl.var$mydata.large.data$MYDATA.PVAL <- -log10(gbl.var$mydata.large.data$MYDATA.PVAL) # LOG Variable
    
  } else {
    
    stop("Invalid lab.Y: ", config.var$lab.Y, ". Specify either 'log' or 'ln'.\n")
  }
  
  #DEBUG STATEMENT
  if (config.var$verbose)  cat("FINISH FIX.VALUES\n")
  
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
    if (gbl.var$verbose)  cat("START CHECK.FORMAT.MYDATA.large\n")
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
  if ((format == "DTR" | format == "dtr_asso" | format == "site" | format == "site_asso")) {
    if(is.null(mydata.test$LOC)) {
      
      stop("Missing MYDATA data file column LOC\n")
    }
    if ((length(grep("^[a-zA-Z]",mydata.test$LOC)) > 0) ){
      
      stop("Missing column LOC has to be number\n")
    }
  } else if (format == "region" | format == "region_asso") {
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
  if (format == "region_asso" | format == "dtr_asso" |  format == "site_asso") {
    if(is.null(mydata.test$MYDATA.ASSO) ) {
      
      stop("Missing MYDATA data file column MYDATA.ASSO\n")
    }
    mydata.test$MYDATA.ASSO.ORIGINAL <- mydata.test$MYDATA.ASSO
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
  # if (format == "region_asso" | format == "DTR_ASSO" | format == "site_asso") {
  #    if (gbl.var$verbose)  cat("Association ",gbl.var$mydata.large.data$Association[1]," \n") 
  # }
  
  #Sort DATA from localisation
  if(option == 1) {
    mydata.test.nosort <- gbl.var$mydata.data 
    mydata.test.sort <- mydata.test.nosort
    # if (gbl.var$verbose)  cat("test",grepl( "region",format)[1],"\n")
    if (grepl( "region",format)[1]) {
      mydata.test.sort <- mydata.test.nosort[order(mydata.test.nosort$LOC.START),]
    } else {
      mydata.test.sort <- mydata.test.nosort[order(mydata.test.nosort$LOC),]
    }
    gbl.var$mydata.data <- mydata.test.sort
  }else {
    mydata.test.nosort <- gbl.var$mydata.large.data 
    mydata.test.sort <- mydata.test.nosort
    if (grepl( "region",format)[1]) {
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
  if (config.var$verbose)  cat("START READ.FILE.MYDATA ",numfile,"\n")
  #  if (config.var$verbose)  cat("Format",config.var$mydata.format,"\n")
  
  #  split.mydata.file <-gbl.var$split.mydata.file
  # numfile <- gbl.var$cur.sample
  if(gbl.var$split.type[[1]][numfile] == "file" | gbl.var$split.type[[1]][numfile] == "listfile" ) {
    general.data.raw <- read.delim(split.mydata.file[[1]][numfile], header=TRUE, sep="\t", as.is=TRUE, blank.lines.skip = TRUE, fill=TRUE)
  } else {
    general.data.raw <- split.mydata.file[[numfile]]
  }
  gbl.var$general.data <- general.data.raw
  
  #----------- MYDATA table with pvalue
  if (config.var$mydata.format == "DTR") {
    gbl.var$mydata.data <- gbl.var$general.data[,c(1,2,7,14)]
    colnames(gbl.var$mydata.data)<-c("MYDATA.NAME","CHROMOSOME","LOC","MYDATA.PVAL")
  } else if (config.var$mydata.format == "dtr_asso") {
    gbl.var$mydata.data <- gbl.var$general.data[,c(1,2,7,14,15)]
    colnames(gbl.var$mydata.data)<-c("MYDATA.NAME","CHROMOSOME","LOC","MYDATA.PVAL","MYDATA.ASSO")
  } else if (config.var$mydata.format == "site") {
    gbl.var$mydata.data <- gbl.var$general.data
    colnames(gbl.var$mydata.data)<-c("MYDATA.NAME","CHROMOSOME","LOC","MYDATA.PVAL")
  } else if (config.var$mydata.format == "site_asso") {
    gbl.var$mydata.data <- gbl.var$general.data
    colnames(gbl.var$mydata.data)<-c("MYDATA.NAME","CHROMOSOME","LOC","MYDATA.PVAL","MYDATA.ASSO")
  } else if (config.var$mydata.format == "region") {
    gbl.var$mydata.data <- gbl.var$general.data
    colnames(gbl.var$mydata.data)<-c("MYDATA.NAME","CHROMOSOME","LOC.START","LOC.END","MYDATA.PVAL")
  } else if (config.var$mydata.format == "region_asso") {
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
  if (config.var$verbose)  cat("FINISH READ.FILE.MYDATA ",numfile,"\n")
  return(gbl.var)
}

#------------------- READ LARGE FILE OF POSITION and PVALUE-------
read.file.mydata.large <- function(large.split.mydata.file, config.var, gbl.var,numfile.large){
  #DEBUG STATEMENT
  if (config.var$verbose)  cat("START READ.FILE.MYDATA.large ",numfile.large," \n")
  # if (config.var$verbose)  cat("Format",config.var$mydata.large.format,"...", numfile.large,"\n")
  
  if(gbl.var$large.split.type == "listfile"){
    general.data.large.raw.tmp <- read.delim(large.split.mydata.file[[1]][numfile.large], header=TRUE, sep="\t", as.is=TRUE, blank.lines.skip = TRUE, fill=TRUE)
  } else {
    general.data.large.raw.tmp <- large.split.mydata.file[[numfile.large]]
  }
  large.format <- gbl.var$large.split.format[[1]][numfile.large]
  
  #----------- MYDATA table with pvalue
  if (large.format == "DTR") {
    general.data.large.raw <- general.data.large.raw.tmp[,c(1,2,7,14)]
    colnames(general.data.large.raw)<-c("MYDATA.NAME","CHROMOSOME","LOC","MYDATA.PVAL")
  } else if (large.format == "dtr_asso") {
    general.data.large.raw <- general.data.large.raw.tmp[,c(1,2,7,14,15)]
    colnames(general.data.large.raw)<-c("MYDATA.NAME","CHROMOSOME","LOC","MYDATA.PVAL","MYDATA.ASSO")
  } else if (large.format == "site") {
    general.data.large.raw <- general.data.large.raw.tmp
    colnames(general.data.large.raw)<-c("MYDATA.NAME","CHROMOSOME","LOC","MYDATA.PVAL")
  } else if (large.format == "site_asso") {
    general.data.large.raw <- general.data.large.raw.tmp
    colnames(general.data.large.raw)<-c("MYDATA.NAME","CHROMOSOME","LOC","MYDATA.PVAL","MYDATA.ASSO")
  }else if (large.format == "region") {
    general.data.large.raw <- general.data.large.raw.tmp
    colnames(general.data.large.raw)<-c("MYDATA.NAME","CHROMOSOME","LOC.START","LOC.END","MYDATA.PVAL")
  } else if (large.format == "region_asso") {
    general.data.large.raw <- general.data.large.raw.tmp
    colnames(general.data.large.raw)<-c("MYDATA.NAME","CHROMOSOME","LOC.START","LOC.END","MYDATA.PVAL","MYDATA.ASSO")
  }else{
    
    stop("This format does not exist\n")
  }
  
  
  if (!grepl( "^chr",general.data.large.raw$CHROMOSOME)[1]) {
    if (config.var$verbose)  cat("Missing column CHROMOSOME should be at UCSC format\n")
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
  if (config.var$verbose)  cat("FINISH READ.FILE.MYDATA.large \n")
  return(gbl.var)
}


#------------------- READ CORMATRIX FILE -------
read.file.cormatrix <- function(config.var, gbl.var,split.cormatrix.file=NULL){
  #DEBUG STATEMENT
  if (config.var$verbose)  cat("START READ.FILE.COMATRIX\n")
  
  if(is.null(split.cormatrix.file)){
    #  if (config.var$verbose)  cat("CORMATRIX empty \n")
    
    stop("EMPTY CORMATRIX data file \n")
  } else {
    if (config.var$verbose)  cat("Format data",config.var$cormatrix.format,"\n")
    
    #----------- MYDATA table with pvalue
    if (config.var$cormatrix.format == "DTR_CORMATRIX") {
      cormatrix.data.raw <- gbl.var$general.data[,-c(1:19)]
      #This format is only possible for comet and comet.web function
      if(gbl.var$presence.mydata == 0){
        #--- Check the same size cormatrix and data
        if((nrow(gbl.var$mydata.data) == nrow(cormatrix.data.raw)) & (nrow(gbl.var$mydata.data) == ncol(cormatrix.data.raw))){
          if(config.var$verbose){
            cat("Same number of data in the primary data and correlation matrix\n")
          }
          # cormatrix.data.raw[match(gbl.var$mydata.data$MYDATA.NAME, colnames(cormatrix.data.raw)),]
          #  cormatrix.data.raw[gbl.var$mydata.data$MYDATA.NAME,]
          if(identical(as.vector(gbl.var$mydata.data$MYDATA.NAME),as.vector(colnames(cormatrix.data.raw)))) {
            gbl.var$cormatrix.data <-as.matrix(cormatrix.data.raw)
          } else {
            
            stop("Invalid MYDATA data file: ", config.var$mydata.file, " and CORMATRIX ",split.cormatrix.file[[1]][1] ," do not have the same omic sites/regions or not in the same order, need to put the ascendant order of their start position in the two files.\n Error: ",all.equal(as.vector(gbl.var$mydata.data$MYDATA.NAME),as.vector(colnames(cormatrix.data.raw))),"\n")
          }
        }else {
          
          stop("Invalid MYDATA data file: ", config.var$mydata.file, " and CORMATRIX ",split.cormatrix.file[[1]][1] ," do not have the same size\n")
        }
      }
    } else if (config.var$cormatrix.format == "cormatrix") {
      #This format is only possible for comet and comet.web function
      if(gbl.var$split.cormatrix.type == "listfile"){
        cormatrix.data.raw<- read.delim(split.cormatrix.file[[1]][1], sep="\t", header=TRUE, as.is=TRUE, blank.lines.skip = TRUE, fill=TRUE)
      } else{
        cormatrix.data.raw<-split.cormatrix.file[[1]]
      }
      
      if(gbl.var$presence.mydata == 0){
        ## For comet and comet.web
        #--- Check the same size cormatrix and data
        if((nrow(gbl.var$mydata.data) == nrow(cormatrix.data.raw)) & (nrow(gbl.var$mydata.data) == ncol(cormatrix.data.raw))){
          # cormatrix.data.raw[match(gbl.var$mydata.data$MYDATA.NAME, colnames(cormatrix.data.raw)),]
          #  cormatrix.data.raw[gbl.var$mydata.data$MYDATA.NAME,]
          if(identical(as.vector(gbl.var$mydata.data$MYDATA.NAME),as.vector(colnames(cormatrix.data.raw)))) {
            gbl.var$cormatrix.data <-as.matrix(cormatrix.data.raw)
          } else {
            
            stop("Invalid MYDATA data file: ", config.var$mydata.file, " and CORMATRIX ",split.cormatrix.file[[1]][1] ," do not have the same omic sites/regions or not in the same order, need to put the ascendant order of their start position in the two files.\n Error: ",all.equal(as.vector(gbl.var$mydata.data$MYDATA.NAME),as.vector(colnames(cormatrix.data.raw))),"\n")
          }
        }else {
          
          stop("Invalid MYDATA data file: ", config.var$mydata.file, " and CORMATRIX ",split.cormatrix.file[[1]][1] ," do not have the same size\n")
        }
      } else {
        ## For comet.list
        gbl.var$cormatrix.data <-as.matrix(cormatrix.data.raw)
      }
      
    } else if (config.var$cormatrix.format == "raw_rev"){
      # if (config.var$verbose)  cat("Format data",config.var$cormatrix.format,"\n")
      if(gbl.var$split.cormatrix.type == "listfile"){
        matrix.data.raw_rot<- read.delim(split.cormatrix.file[[1]][1], header=TRUE, sep="\t", as.is=TRUE, blank.lines.skip = TRUE, fill=TRUE)
      } else {
        matrix.data.raw_rot <- split.cormatrix.file[[1]][1]
      }
      matrix.data.raw <- matrix.data.raw_rot
      if(gbl.var$presence.mydata == 0){
        ## For comet and comet.web
        matrix.data.raw <- matrix.data.raw[match(gbl.var$mydata.data$MYDATA.NAME, rownames(matrix.data.raw)),]
      }
      #Update the rownames
      rownames(matrix.data.raw) <-matrix.data.raw[,1]
      matrix.data.raw <-matrix.data.raw[,-1]
      gbl.var$matrix.data <-matrix.data.raw
      nr=nrow(matrix.data.raw)
      if(gbl.var$presence.mydata == 0){
        ## For comet and comet.web
        nr1=nrow(gbl.var$mydata.data)
        #if (config.var$verbose)  cat("Matrice size",nr,"\n")
        #if (config.var$verbose)  cat("Matrice size",nr1,"\n")
        #--- Check the same size cormatrix and data
        if(nrow(gbl.var$mydata.data) == nrow(matrix.data.raw) ){
          gbl.var$matrix.data <-matrix.data.raw
          gbl.var<-compute.cormatrix(config.var, gbl.var)
          # gbl.var <-compute.pvalue.cormatrix(config.var, gbl.var)
        }else {
          
          stop("Invalid MYDATA data file: ", config.var$mydata.file," value (" ,nr1,") and CORMATRIX ",split.cormatrix.file[[1]][1] ," value (" ,nr,") do not have the same size\n")
        }
      } else {
        ## For comet.list
        gbl.var$matrix.data <-matrix.data.raw
        gbl.var<-compute.cormatrix(config.var, gbl.var)
        #gbl.var <-compute.pvalue.cormatrix(config.var, gbl.var)
      }
      
    }else if (config.var$cormatrix.format == "raw"){
      if(gbl.var$split.cormatrix.type == "listfile"){
        matrix.data.raw_rot<- read.delim(split.cormatrix.file[[1]][1], header=TRUE, sep="\t", as.is=TRUE, blank.lines.skip = TRUE, fill=TRUE)
      } else {
        matrix.data.raw_rot <- split.cormatrix.file[[1]]
      }
      matrix.data.raw <- t(matrix.data.raw_rot)
      if(gbl.var$presence.mydata == 0){
        ## For comet and comet.web
        matrix.data.raw <- matrix.data.raw[match(gbl.var$mydata.data$MYDATA.NAME, rownames(matrix.data.raw)),]
      }
      gbl.var$matrix.data <-matrix.data.raw
      nr=nrow(matrix.data.raw)
      if (config.var$verbose)  cat("Matrice size",nr,"\n")
      if(gbl.var$presence.mydata == 0){
        ## For comet and comet.web
        nr1=nrow(gbl.var$mydata.data)
        if (config.var$verbose)  cat("Matrice size",nr1,"\n")
        #--- Check the same size cormatrix and data
        if(nrow(gbl.var$mydata.data) == nrow(matrix.data.raw) ){
          gbl.var$matrix.data <-matrix.data.raw
          gbl.var<-compute.cormatrix(config.var, gbl.var)
          #   gbl.var <-compute.pvalue.cormatrix(config.var, gbl.var)
        }else {
          
          stop("Invalid MYDATA data file: ", config.var$mydata.file," value (" ,nr1,") and CORMATRIX ",split.cormatrix.file[[1]][1] ," value (" ,nr,") do not have the same size\n")
        }
      } else {
        ## For comet.list
        gbl.var$matrix.data <-matrix.data.raw
        gbl.var<-compute.cormatrix(config.var, gbl.var)
        #  gbl.var <-compute.pvalue.cormatrix(config.var, gbl.var)
      }
    } else if (config.var$cormatrix.format == "DTR_RAW"){
      matrix.data.raw.tmp <- gbl.var$general.data[,-c(1:19)]
      rownames(matrix.data.raw.tmp) <- gbl.var$general.data[,1]
      nr=nrow(matrix.data.raw.tmp)
      if(gbl.var$presence.mydata == 0){
        ## For comet and comet.web
        nr1=nrow(gbl.var$mydata.data)
        # if (config.var$verbose)  cat("Matrice size",nr,"\n")
        # if (config.var$verbose)  cat("Matrice size",nr1,"\n")
        #--- Check the same size cormatrix and data
        if(nrow(gbl.var$mydata.data) == nrow(matrix.data.raw) ){
          gbl.var$matrix.data <-matrix.data.raw.tmp
          gbl.var<-compute.cormatrix(config.var, gbl.var)
          #  gbl.var <-compute.pvalue.cormatrix(config.var, gbl.var)
        }else {
          
          stop("Invalid MYDATA data file: ", config.var$mydata.file, " and MATRIX do not have the same size\n")
        }
      } else {
        ## For comet.list
        gbl.var$matrix.data <-matrix.data.raw.tmp
        gbl.var<-compute.cormatrix(config.var, gbl.var)
        #  gbl.var <-compute.pvalue.cormatrix(config.var, gbl.var)
      }
    }
    
    # checking size of data
    num.mydata<-0
    if(gbl.var$presence.mydata == 0){
      ## For comet and comet.web
      num.mydata<-dim(gbl.var$mydata.data)[1]
    }
    row.LD<-dim(gbl.var$cormatrix.data)[1]
    col.LD<-dim(gbl.var$cormatrix.data)[2]
    
    if (config.var$verbose)  cat("Matrice size num.mydat a",num.mydata,"\n")
    if (config.var$verbose)  cat("Matrice size row.LD ",row.LD,"\n")
    if (config.var$verbose)  cat("Matrice size col.LD ",col.LD,"\n")
    
    if(gbl.var$presence.mydata == 0){
      ## For comet and comet.web
      
      if (row.LD != num.mydata | col.LD != num.mydata ) {
        
        stop("Not the same number of row or column between the correlation matrice and MYDATA data\n")
      }  
      
      rownames(gbl.var$cormatrix.data) <- gbl.var$mydata.data$MYDATA.NAME
    }
    #--------- REDUCE the matrix to triangle ---
    gbl.var$cormatrix.data.full <- gbl.var$cormatrix.data
    gbl.var$cormatrix.data[lower.tri(gbl.var$cormatrix.data)] <- NA
    diag(gbl.var$cormatrix.data) <- NA
    
    gbl.var$cormatrix.pvalue.data.full <- gbl.var$cormatrix.pvalue.data
    #keep only the pvalue adjusted 
    #(upper matrix=adjusted pvalue; lower matrix=raw pvalue because of psych library)
    gbl.var$cormatrix.pvalue.data[lower.tri(gbl.var$cormatrix.pvalue.data)] <- NA
    diag(gbl.var$cormatrix.pvalue.data) <- NA
    
    
    #--------- VECTOR for the reference gene ---
    if(gbl.var$presence.mydata == 0){
      ## For comet and comet.web
      best.position <- gbl.var$mydata.best.position
      gbl.var$ref <-gbl.var$cormatrix.data.full[,best.position]
      if (config.var$verbose)  cat("Ref position",best.position,"\n")
    }
  }
  #DEBUG STATEMENT
  if (config.var$verbose)  cat("FINISH READ.FILE.COMET\n")
  return(gbl.var)
}

#------------------- Compute correlation between methylation  -------
compute.cormatrix <- function(config.var, gbl.var){
  #DEBUG STATEMENT
  if (config.var$verbose)  cat("START COMPUTE.CORMATRIX\n")
  
  matrix.data.raw <- gbl.var$matrix.data
  if(gbl.var$presence.mydata == 0){
    ## For comet and comet.web
    if(!(nrow(gbl.var$mydata.data) == nrow(matrix.data.raw) )){
      
      stop("Invalid MYDATA data file: ", config.var$mydata.file, " and MATRIX do not have the same size\n")
    }
  }
  
  mat <- t(matrix.data.raw)
  numrow<- nrow(mat)
  # if (config.var$verbose) cat("Number row",numrow,"\n")
  if (config.var$cormatrix.method == "spearman" | 
      config.var$cormatrix.method == "pearson" |
      config.var$cormatrix.method == "kendall") {
    tmp <- corr.test(mat, adjust=config.var$cormatrix.adjust,
                     method=config.var$cormatrix.method,
                     alpha=config.var$cormatrix.conf.level)
    
    cor.mat <- as.matrix(tmp$r)
    p.mat <-as.matrix(tmp$p)
    CI.mat <-tmp$ci
    
    numrowcor<- nrow(cor.mat)
    #if (config.var$verbose) cat("Number row correlation",numrowcor,"\n")
    
    colnames(cor.mat) <- rownames(cor.mat) <- rownames(matrix.data.raw)
    colnames(p.mat) <- rownames(p.mat) <- rownames(matrix.data.raw)
    
    gbl.var$cormatrix.data <- cor.mat
    gbl.var$cormatrix.pvalue.data <- p.mat
    gbl.var$cormatrix.CI.data <- CI.mat
    
  }else {
    
    stop("Invalid CORMATRIX method : ", config.var$cormatrix.method, "\n")
  }
  
  #DEBUG STATEMENT
  if (config.var$verbose)  cat("FINISH COMPUTE.CORMATRIX\n")
  return(gbl.var)
}


#------------------- Compute list of tracks at GVIZ format  -------
createList.trackUser <- function(config.var, gbl.var){
  #DEBUG STATEMENT
  if (config.var$verbose)  cat("START CREATELIST.TRACKUSER\n")
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
  if (config.var$verbose)  cat("FINISH CREATELIST.TRACKUSER\n")
  return(gbl.var)
}

#-------------------SET PARAMETERS----------------------------

set.image.parameters <- function(config.var, gbl.var) {
  
  #DEBUG STATEMENT
  if (config.var$verbose)  cat("START SET.IMAGE.PARAMETERS\n")
  
  sec.cex.factor <- 1
  
  #DEBUG STATEMENT
  # if (config.var$verbose)  cat("image.size ", config.var$image.size, "\n")
  
  if(config.var$image.size == 3.5) {
    font.size <- 4
    line.width <- 0.5
    sec.cex.factor <- 1.5
    cex.factor.symbol <- 0.25
  } else if(config.var$image.size == 7) {
    font.size <- 9
    line.width <- 1
    sec.cex.factor <- 2
    cex.factor.symbol <- 0.5
  } else {
    
    stop("Invalid image size: ", config.var$image.size, "\n")
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
  # if (config.var$verbose)  cat("gbl.var$mydata.num ", gbl.var$mydata.num, "\n")
  if (config.var$verbose)  cat("font.size ", font.size, "\n")
  # if (config.var$verbose)  cat("line.width ", line.width, "\n")
  # if (config.var$verbose)  cat("cex.factor.symbol ", cex.factor.symbol, "\n")
  # if (config.var$verbose)  cat("cex.factor ", cex.factor, "\n")
  # if (config.var$verbose)  cat("sec.cex.factor ", sec.cex.factor, "\n")
  
  #DEBUG STATEMENT
  if (config.var$verbose)  cat("FINISH SET.IMAGE.PARAMETERS\n")
  
  return(list(font.size = font.size, line.width = line.width, cex.factor.symbol = cex.factor.symbol, cex.factor = cex.factor))
}
