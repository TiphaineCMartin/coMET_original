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
# File: DrawPlot.R
# Author: Tiphaine Martin
# Email: tiphaine.martin@kcl.ac.uk
# Purpose: coMET allows the display of p-values from association
#           with a correlation heatmap.
# Version : 0.99.9
###########################################################################

#------------------- DRAW Structure + all tracks ---------------
draw.plot.comet.web <- function(config.var, gbl.var,newpage=TRUE){
  #DEBUG STATEMENT
  if (config.var$verbose) cat("START DRAW.PLOT.COMET.WEB\n")
  
  if(newpage == "TRUE"){
    grid.newpage()
  }
  
  #---------------- DRAW DIFFERENT TRACKS ---------
  #-------------DRAW the grid for the windows
  
  gbl.var <- draw.plot.grid.setup(config.var, gbl.var)
  
  #------------ DRAW LARGE MYDATA DOT in plot
  if(! is.null(config.var$mydata.large.file)){
    for(i in 1:gbl.var$large.mydata.samples) {
      
      gbl.var$cur.sample.large <- i
      large.split.mydata.file=gbl.var$large.split.mydata.file
      #-------------------- READ FILE CpG DATA
      gbl.var <-read.file.mydata.large(large.split.mydata.file,config.var, gbl.var,i)
      
      #----- FIX VALUES -------------------------
      
      fix.var <- fix.values.large(config.var, gbl.var)
      
      config.var <- fix.var$config.var
      gbl.var <- fix.var$gbl.var
      
      #------------ DRAW PVALUE of DATA
      if(config.var$disp.mydata) {
        draw.plot.grid.mydata.large(config.var, gbl.var)
      }
    }
  }
  
  #------------ DRAW MYDATA DOT in plot
  for(i in 1:gbl.var$mydata.samples) {
    
    gbl.var$cur.sample <- i
    split.mydata.file=gbl.var$split.mydata.file
    #-------------------- READ FILE CpG DATA
    gbl.var <-read.file.mydata(split.mydata.file,config.var, gbl.var,i)
    
    #----- FIX VALUES -------------------------
    
    fix.var <- fix.values(config.var, gbl.var)
    
    config.var <- fix.var$config.var
    gbl.var <- fix.var$gbl.var
    
    #------------ DRAW PVALUE of DATA
    if(config.var$disp.mydata) {
      draw.plot.grid.mydata(config.var, gbl.var)
    }
  }
  
  #------------ DRAW LEGEND 
  if(config.var$disp.legend) {
    draw.legend(config.var, gbl.var)
  }
  
  #------------ DRAW CORMATRIXMAP
  if(config.var$disp.cormatrixmap) {
    draw.plot.cormatrix.plot(config.var, gbl.var)
  }
  
  #------------ DRAW BIOFEATURES : BIOFEAT and OTHER ANNOTATION 
  # draw annotation tracks
  draw.plot.annotation(config.var, gbl.var)
  
  if(has.key("geneENSEMBL", gbl.var$split.list.tracks) | 
     (has.key("transcriptENSEMBL", gbl.var$split.list.tracks))
     | (has.key("genesUCSC",gbl.var$split.list.tracks) )) {
    # draw name of genes
    draw.name.genes.web(config.var, gbl.var)
  }
  # draw name of tracks
  draw.name.tracks.web(config.var,gbl.var)
  
  
  
  #DEBUG STATEMENT
  if (config.var$verbose) cat("FINISH DRAW.PLOT\n")
  return(gbl.var)
}

#------------------- DRAW all tracks + Generic---------------
draw.plot.comet <- function(config.var, gbl.var,newpage=TRUE){
  #DEBUG STATEMENT
  if (config.var$verbose) cat("START DRAW.PLOT.COMET\n")
  
  if(newpage == "TRUE"){
    grid.newpage()
  }
  
  #---------------- DRAW DIFFERENT TRACKS ---------
  #-------------DRAW the grid for the windows
  
  gbl.var <- draw.plot.grid.setup(config.var, gbl.var)
  
  #------------ DRAW LARGE MYDATA DOT in plot
  if(! is.null(config.var$mydata.large.file)){
    for(i in 1:gbl.var$large.mydata.samples) {
      
      gbl.var$cur.sample.large <- i
      large.split.mydata.file=gbl.var$large.split.mydata.file
      #-------------------- READ FILE CpG DATA
      gbl.var <-read.file.mydata.large(large.split.mydata.file,config.var, gbl.var,i)
      
      #----- FIX VALUES -------------------------
      
      fix.var <- fix.values.large(config.var, gbl.var)
      
      config.var <- fix.var$config.var
      gbl.var <- fix.var$gbl.var
      
      #------------ DRAW PVALUE of DATA
      if(config.var$disp.mydata) {
        draw.plot.grid.mydata.large(config.var, gbl.var)
      }
    }
  }
  
  #------------ DRAW MYDATA DOT in plot
  for(i in 1:gbl.var$mydata.samples) {
    
    gbl.var$cur.sample <- i
    split.mydata.file=gbl.var$split.mydata.file
    #-------------------- READ FILE CpG DATA
    gbl.var <-read.file.mydata(split.mydata.file,config.var, gbl.var,i)
    
    #----- FIX VALUES -------------------------
    
    fix.var <- fix.values(config.var, gbl.var)
    
    config.var <- fix.var$config.var
    gbl.var <- fix.var$gbl.var
    
    #------------ DRAW PVALUE of DATA
    if(config.var$disp.mydata) {
      draw.plot.grid.mydata(config.var, gbl.var)
    }
  }
  
  #------------ DRAW LEGEND 
  if(config.var$disp.legend) {
    draw.legend(config.var, gbl.var)
  }
  
  #------------ DRAW CORMATRIXMAP
  if(config.var$disp.cormatrixmap) {
    draw.plot.cormatrix.plot(config.var, gbl.var)
  }
  
  #------------ DRAW BIOFEATURES : BIOFEAT and OTHER ANNOTATION 
  # draw annotation tracks
  draw.plot.annotation(config.var, gbl.var)
  
  
  config.var <- fix.var$config.var
  gbl.var <- fix.var$gbl.var
  
  #DEBUG STATEMENT
  if (config.var$verbose) cat("FINISH DRAW.PLOT.COMET\n")
  
  return(gbl.var)
}

#------------------- DRAW all tracks + Generic + No PVALUE PLOT---------------
draw.plot.comet.nopval <- function(config.var, gbl.var,newpage=TRUE){
  #DEBUG STATEMENT
  if (config.var$verbose) cat("START DRAW.PLOT.COMET.NOPVAL\n")
  
  if(newpage == "TRUE"){
    grid.newpage()
  }
  
  #---------------- DRAW DIFFERENT TRACKS ---------
  #-------------DRAW the grid for the windows
  
  gbl.var <- draw.plot.grid.setup(config.var, gbl.var)
  
  #------------ DRAW MYDATA DOT in plot
  for(i in 1:gbl.var$mydata.samples) {
    
    gbl.var$cur.sample <- i
    split.mydata.file=gbl.var$split.mydata.file
    #-------------------- READ FILE CpG DATA
    gbl.var <-read.file.mydata(split.mydata.file,config.var, gbl.var,i)
    
    #----- FIX VALUES -------------------------
    
    fix.var <- fix.values(config.var, gbl.var)
    
    config.var <- fix.var$config.var
    gbl.var <- fix.var$gbl.var
  }
  
  #------------ DRAW LEGEND 
  if(config.var$disp.legend) {
    draw.legend(config.var, gbl.var)
  }
  
  #------------ DRAW CORMATRIXMAP
  if(config.var$disp.cormatrixmap) {
    draw.plot.cormatrix.plot(config.var, gbl.var)
  }
  
  #------------ DRAW BIOFEATURES : BIOFEAT and OTHER ANNOTATION 
  # draw annotation tracks
  draw.plot.annotation(config.var, gbl.var)
  
  
  
  config.var <- fix.var$config.var
  gbl.var <- fix.var$gbl.var
  
  #DEBUG STATEMENT
  if (config.var$verbose) cat("FINISH DRAW.PLOT.COMET.NOPVAL\n")
  
  return(gbl.var)
}

#------------------- PRINT in File ---------------
printPlot.comet.web <- function(config.var, gbl.var){
  #START IMAGE CAPTURE
  imagetype = config.var$image.type
  imagename = config.var$image.name
  imagesize = config.var$image.size
  imagetitle = config.var$image.title
  if(imagetype == "pdf") {
    
    tmp.image.name <- paste(imagename, ".pdf", sep="")
    pdf(encoding = "ISOLatin1.enc",
        file = tmp.image.name,
        onefile=FALSE,
        width=as.numeric(imagesize),
        height=as.numeric(imagesize),
        paper="special")
  } else if(imagetype == "eps") {
    
    tmp.image.name <- paste(imagename, ".eps", sep="")
    postscript(encoding = "ISOLatin1.enc",
               file = tmp.image.name,
               horizontal=FALSE,
               onefile=FALSE,
               width=as.numeric(imagesize),
               height=as.numeric(imagesize),
               paper="special",
               pagecentre=TRUE,
               fonts=c("sans"))
  } else {
    stop("Invalid image type: ", imagetype, " and image size: ", 
         imagesize, " combination.\n")
  }
  #END IMAGE CAPTURE
  
  gbl.var <- draw.plot.comet.web(config.var, gbl.var,newpage=TRUE)
  
  
  #make sure the graphics device is closed and that the remaining device is the null device
  while(dev.cur()[[1]] != 1) {
    dev.off()
  }
}

#------------------- PRINT in File draw generic ---------------
printPlot.comet <- function(config.var, gbl.var){
  #START IMAGE CAPTURE
  imagetype = config.var$image.type
  imagename = config.var$image.name
  imagesize = config.var$image.size
  imagetitle = config.var$image.title
  if(imagetype == "pdf") {
    
    tmp.image.name <- paste(imagename, ".pdf", sep="")
    pdf(encoding = "ISOLatin1.enc",
        file = tmp.image.name,
        onefile=FALSE,
        width=as.numeric(imagesize),
        height=as.numeric(imagesize),
        paper="special")
  } else if(imagetype == "eps") {
    
    tmp.image.name <- paste(imagename, ".eps", sep="")
    postscript(encoding = "ISOLatin1.enc",
               file = tmp.image.name,
               horizontal=FALSE,
               onefile=FALSE,
               width=as.numeric(imagesize),
               height=as.numeric(imagesize),
               paper="special",
               pagecentre=TRUE,
               fonts=c("sans"))
  } else {
    stop("Invalid image type: ", imagetype, " and image size: ", 
         imagesize, " combination.\n")
  }
  #END IMAGE CAPTURE
  
  gbl.var <- draw.plot.comet(config.var, gbl.var,newpage=TRUE)
  
  
  #make sure the graphics device is closed and that the remaining device is the null device
  while(dev.cur()[[1]] != 1) {
    dev.off()
  }
}


#------------------- PRINT in File draw generic and no pvalue plot ---------------
printPlot.comet.nopval <- function(config.var, gbl.var){
  #START IMAGE CAPTURE
  imagetype = config.var$image.type
  imagename = config.var$image.name
  imagesize = config.var$image.size
  imagetitle = config.var$image.title
  if(imagetype == "pdf") {
    
    tmp.image.name <- paste(imagename, ".pdf", sep="")
    pdf(encoding = "ISOLatin1.enc",
        file = tmp.image.name,
        onefile=FALSE,
        width=as.numeric(imagesize),
        height=as.numeric(imagesize),
        paper="special")
  } else if(imagetype == "eps") {
    
    tmp.image.name <- paste(imagename, ".eps", sep="")
    postscript(encoding = "ISOLatin1.enc",
               file = tmp.image.name,
               horizontal=FALSE,
               onefile=FALSE,
               width=as.numeric(imagesize),
               height=as.numeric(imagesize),
               paper="special",
               pagecentre=TRUE,
               fonts=c("sans"))
  } else {
    stop("Invalid image type: ", imagetype, " and image size: ", 
         imagesize, " combination.\n")
  }
  #END IMAGE CAPTURE
  
  gbl.var <- draw.plot.comet.nopval(config.var, gbl.var,newpage=TRUE)
  
  
  #make sure the graphics device is closed and that the remaining device is the null device
  while(dev.cur()[[1]] != 1) {
    dev.off()
  }
}

#------------------- WRITE in File the correlation between omic features ---------------
write.comet.list <- function(config.var, gbl.var){
  if (config.var$verbose)  cat("START WRITE COMET LIST\n")
  #Create the output file
  cormatrix.output = config.var$cormatrix.output
  # tmp.cormatrix.output <- paste(cormatrix.output, ".txt", sep="")
  
  #Different matrix
  cormatrix.matrix <- gbl.var$cormatrix.data.full
  cormatrix.pvalue.matrix <- gbl.var$cormatrix.pvalue.data.full
  cormatrix.CI.matrix <- gbl.var$cormatrix.CI.data
  
  matrix.rows <- dim(cormatrix.matrix)[1]
  matrix.cols <- dim(cormatrix.matrix)[2]
  
  #Create the data.frame
  listcorrelation <- NULL
  
  listcorrelation <- data.frame(omicfeature1=character(),
                                omcifeature2=character(), 
                                correlation=numeric(),
                                standartDeviation=numeric(),
                                lowerCI=numeric(),
                                upperCI=numeric(),
                                pvalue=numeric(),
                                pvalue.adjust=numeric(),
                                stringsAsFactors = FALSE) 
  
  r <- 0
  for (i in 1:(matrix.rows - 1)) {
    for (j in (i + 1):matrix.rows) {
      r <- r+1
      omic1 <- rownames(cormatrix.matrix)[i]
      omic2 <- rownames(cormatrix.matrix)[j]
      cor12 <- cormatrix.matrix[i,j]
      pval.adjust12 <- cormatrix.pvalue.matrix[i,j]
      pval12 <- cormatrix.pvalue.matrix[j,i]
      uppCI12 <- cormatrix.CI.matrix[r,"upper"]
      lowCI12 <- cormatrix.CI.matrix[r,"lower"]
      
      if(pval.adjust12 <= config.var$cormatrix.sig.level) {
        onecorrelation <-c(omic1,omic2,cor12,lowCI12,uppCI12,pval12,pval.adjust12)
        listcorrelation <-rbind(listcorrelation,onecorrelation)
      }
      listcorrelation <- data.frame(lapply(listcorrelation, as.character), stringsAsFactors=FALSE)
    }
  }
  
  colnames(listcorrelation)<- c("omicFeature1","omicFeature2","correlation",
                                "lowerCI","upperCI","pvalue","pvalue.adjusted")
  
  write.table(listcorrelation,file=cormatrix.output,row.names = FALSE,
              col.names = TRUE,sep="\t",quote = FALSE)
  
  if (config.var$verbose)  cat("END WRITE COMET LIST\n")
}