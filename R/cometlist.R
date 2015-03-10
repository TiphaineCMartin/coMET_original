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
# File: cometlist.R
# Author: Tiphaine Martin
# Email: Tiphaine.Martin@kcl.ac.uk
# Purpose: coMET.list allows the obtention of co-omic features in a region of interest.
# Version : 0.99.9
###########################################################################

#coMET allows the display of p-values from association with a correlation heatmap. 
#This version coMET.web is used for web site.

# The values given such as attributes of function erase ones of configuration file.
comet.list <- function(cormatrix.file = NULL,
                       cormatrix.method = "spearman",
                       cormatrix.format= "raw",
                       cormatrix.conf.level=0.05,
                       cormatrix.sig.level= 1,
                       cormatrix.adjust="none",
                       cormatrix.type = "listdataframe",
                       cormatrix.output="cormatrix_list",
                       config.file = NULL,
                       verbose = FALSE) {
  
  #-------------------MAIN FUNCTION----------------------------
  
  #DEBUG STATEMENT
  if (verbose) cat("START COMET.LIST\n")
  
  #-------------------GLOBAL VARIABLES BEGINS------------------
  
  #create three variables (global variables), which are lists of all other variables
  #used throughout the program to be passed to functions with values changed as
  #necessary.
  gbl.var <- NULL
  
  #DATA VARIABLES
  presence.mydata<- 1 # 0 presence of omic-WAS results (object mydata), 1 no presence of data
  split.cormatrix.file <- NULL
  
  #FORMATTING VARIABLES
  cormatrix.data <- NULL
  split.cormatrix.type <- NULL
  matrix.data <- NULL
  cormatrix.data.full <- NULL
  cormatrix.pvalue.data <- NULL
  cormatrix.pvalue.data.full <- NULL
  cormatrix.CI.data <- NULL
  
  verbose <- FALSE
  
  gbl.var <- list(
    presence.mydata=presence.mydata,
    cormatrix.data = cormatrix.data,
    split.cormatrix.type = split.cormatrix.type,
    split.cormatrix.file = split.cormatrix.file,
    cormatrix.pvalue.data = cormatrix.pvalue.data,
    cormatrix.pvalue.data.full = cormatrix.pvalue.data.full,
    matrix.data = matrix.data,
    cormatrix.data.full = cormatrix.data.full,
    cormatrix.data.full = cormatrix.data.full,
    cormatrix.CI.data = cormatrix.CI.data,
    verbose = verbose
  )
  
  #-------------------GLOBAL VARIABLES ENDS------------------
  
  
  #-------------------UPDATE CONFIGURATION VARIABLES---------
  config.var <- list(cormatrix.file = cormatrix.file,
                     cormatrix.method = cormatrix.method,
                     cormatrix.format = cormatrix.format,
                     cormatrix.conf.level = cormatrix.conf.level,
                     cormatrix.sig.level = cormatrix.sig.level,
                     cormatrix.adjust = cormatrix.adjust, 
                     cormatrix.type = cormatrix.type,
                     cormatrix.output = cormatrix.output,
                     verbose = verbose)
  
  #-------------------CONFIGURATION VARIABLES BEGINS from config file---------
  gbl.var$verbose <- config.var$verbose
  
  if(!is.null(config.file)) {
    config.var <- read.config(config.file, config.var)
  }
  
  #-------------- CHECK if all parameters have values -----------------
  check.configVar.cometlist(config.var)
  
  #------------- READ DATA and UPDATE VARIABLES
  if(!is.null(config.var$cormatrix.file)){
    if(!is.null(config.var$cormatrix.type)){
      split.cormatrix.type <- strsplit(config.var$cormatrix.type, ",")
      comatrix.type.length <- length(split.cormatrix.type[[1]])
      gbl.var$split.cormatrix.type <- split.cormatrix.type
    } else {
      stop("Need to define the format of your correlation matrix from FILE or from MATRIX")
    }
    if(gbl.var$split.cormatrix.type == "listfile"){
      split.cormatrix.file <- strsplit(config.var$cormatrix.file, ",")
      comatrix.file.length <- length(split.cormatrix.file[[1]])
      gbl.var$split.cormatrix.file <- split.cormatrix.file
    }else {
      split.cormatrix.file <- config.var$cormatrix.file
      comatrix.file.length <- length(split.cormatrix.file)
      gbl.var$split.cormatrix.file <- split.cormatrix.file
    }
  } else {
    stop("No correlation matrice. Impossible to create the list of correlation between omic features")
  }
  
  #---- FIX GENERAL VARIABLES
  fix.var <- retrieve.data.cometlist(config.var, gbl.var)
  config.var <- fix.var$config.var
  gbl.var <- fix.var$gbl.var
  
  
  #############RUN the method
  write.comet.list(config.var, gbl.var)
  
  #END FUNCTION LIST
  
  invisible(list(config.var, gbl.var))
  #DEBUG STATEMENT
  if (config.var$verbose)  cat("FINISH COMET.LIST \n")
}


