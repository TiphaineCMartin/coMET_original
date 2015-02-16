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
# File: GeneralMethodComet.R
# Author: Tiphaine Martin
# Email: tiphaine.martin@kcl.ac.uk
# Purpose: coMET allows the display of p-values from association
#           with a correlation heatmap.
# Version : 0.99.0
###########################################################################

#-------------------DRAW PLOT GRID----------------------------
draw.plot.grid.setup <- function(config.var, gbl.var) {
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("START DRAW.PLOT.GRID.SETUP\n")
  
  #set y-axis initial value
  if(config.var$DISP.PVAL.THRESHOLD != -1) {
    gbl.var$r.y[1] <- config.var$DISP.PVAL.THRESHOLD
    gbl.var$axis.y <- seq(floor(gbl.var$r.y[1]), max(gbl.var$axis.y))
  } else {
    gbl.var$r.y[1] <- 0
  }
  
  #-------------------VIEWPORT BEGINS------------------
  
  #-------------------CORMATRIX VIEWPORT---------------------
  # FOR CORMATRIX MAP
  top.vp.cormatrixmap <- viewport(layout = grid.layout(7, 3, widths = c(0.18, 0.6875, 0.1325),
                                                       heights = c(.05, .025, .18, .24, .001, .5, .05)),
                                  name = "top.vp.cormatrixmap")
  
  
  top.vp.nocormatrixmap <- viewport(layout = grid.layout(6, 3, widths = c(0.18, 0.6875, 0.1325),,
                                                         heights = c(.05, .025, .42, .001, .5, .05)),
                                    name = "top.vp.nocormatrixmap")
  
  top.vp.cormatrixmap.nopval <- viewport(layout = grid.layout(6, 3, widths = c(0.18, 0.6875, 0.1325),
                                                              heights = c(.05, .025, .42, .001, .5, .05)),
                                         name = "top.vp.cormatrixmap.nopval")
  if(config.var$DISP.PVALUEPLOT){
    if(config.var$DISP.CORMATRIXMAP) {
      gbl.var$top.vp <- top.vp.cormatrixmap
    }
    else {
      gbl.var$top.vp <- top.vp.nocormatrixmap
    }
  } else {
    gbl.var$top.vp <- top.vp.cormatrixmap.nopval
  }
  
  #-------------------TITLE VIEWPORT---------------------
  title.vp.cormatrixmap <- viewport(height = .8,
                                    width = .8,
                                    layout.pos.row = 1,
                                    layout.pos.col = 2,
                                    name = "title.vp.cormatrixmap")
  title.vp.nocormatrixmap <- viewport(height = .8,
                                      width = .8,
                                      layout.pos.row = 1,
                                      layout.pos.col = 2,
                                      name = "title.vp.nocormatrixmap")
  title.vp.cormatrixmap.nopval <- viewport(height = .8,
                                           width = .8,
                                           layout.pos.row = 1,
                                           layout.pos.col = 2,
                                           name = "title.vp.cormatrixmap.nopval")
  if(config.var$DISP.PVALUEPLOT){
    if(config.var$DISP.CORMATRIXMAP) {
      pushViewport(vpTree(top.vp.cormatrixmap, vpList(title.vp.cormatrixmap )))
    }
    else {
      pushViewport(vpTree(top.vp.nocormatrixmap, vpList(title.vp.nocormatrixmap )))
    }
  }else{
    pushViewport(vpTree(top.vp.cormatrixmap.nopval, vpList(title.vp.cormatrixmap.nopval)))
  }
  image.title.text <- textGrob(config.var$IMAGE.TITLE, x = 0.5, y = 0.5, just = c("center"),
                               default.units = "native", gp = gpar(fontsize = gbl.var$font.size + 4, fontface = "bold"))
  image.title <- gTree(children=gList(image.title.text), vp=title.vp.cormatrixmap, name="image.title")
  
  grid.draw(image.title)
  
  #--------------------- PLOT CHROMOSOME ------------------------
  popViewport()
  
  chromosome.vp.cormatrixmap <- viewport(height = 1,
                                         width = 1,
                                         layout.pos.row = 2,
                                         layout.pos.col = 2,
                                         name = "chromosome.vp.cormatrixmap")
  chromosome.vp.nocormatrixmap <- viewport(height = 1,
                                           width = 1,
                                           layout.pos.row = 2,
                                           layout.pos.col = 2,
                                           name = "chromosome.vp.nocormatrixmap")
  chromosome.vp.cormatrixmap.nopval <- viewport(height = 1,
                                                width = 1,
                                                layout.pos.row = 2,
                                                layout.pos.col = 2,
                                                name = "chromosome.vp.cormatrixmap.nopval")
  
  if(config.var$DISP.PVALUEPLOT){
    if(config.var$DISP.CORMATRIXMAP) {
      pushViewport(vpTree(top.vp.cormatrixmap, vpList(chromosome.vp.cormatrixmap )))
    }
    else {
      pushViewport(vpTree(top.vp.nocormatrixmap, vpList(chromosome.vp.nocormatrixmap )))
    }
  } else {
    pushViewport(vpTree(top.vp.cormatrixmap.nopval, vpList(chromosome.vp.cormatrixmap.nopval )))
  }
  
  itrack <- IdeogramTrack(genome=gbl.var$mydata.gen, chromosome=gbl.var$mydata.chr)
  plotTracks(itrack,from=gbl.var$min.x,to=gbl.var$max.x,panel.only=TRUE, fontsize=5,showBandId = TRUE, cex.bands = 0.5)
  
  #-------------------CONNECTOR VIEWPORT---------------------
  
  if(config.var$DISP.CORMATRIXMAP) {
    if(config.var$DISP.PVALUEPLOT){
      draw.plot.linesconnection(top.vp.cormatrixmap, config.var, gbl.var)
      
    } else {
      draw.plot.linesconnection(top.vp.cormatrixmap.nopval, config.var, gbl.var)
    }
  }
  else {
    popViewport()
  }
  
  #-------------------DATA VIEWPORT---------------------
  
  if(config.var$DISP.CORMATRIXMAP & config.var$DISP.PVALUEPLOT) {
    draw.plot.axis.data(top.vp.cormatrixmap, config.var, gbl.var)
  }
  else {
    if(config.var$DISP.CORMATRIXMAP==FALSE & config.var$DISP.PVALUEPLOT){
      draw.plot.axis.data(top.vp.nocormatrixmap, config.var, gbl.var)
    } else {
      popViewport()
    }
    
  }
  
  #--------------------- POSITION of LEGEND PHYSICAL DISTANCE ------------------------
  #display the physical distance on plots without LD maps
  if(config.var$DISP.PHYS.DIST & !config.var$DISP.CORMATRIXMAP) {
    
    if(gbl.var$total.dist > 1000) {
      grid.text(paste("Physical Distance: ", round(gbl.var$total.dist/1000, 1), " kb", sep=""),
                x = 0.5,
                y = -1.15,
                just = "center",
                gp = gpar(fontsize = gbl.var$font.size))
    } else {
      grid.text(paste("Physical Distance: ",
                      gbl.var$total.dist, " bases", sep=""),
                x = 0.5, y = -1.15,
                just = "center",
                gp = gpar(fontsize = gbl.var$font.size))
    }
  }
  
  
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("FINISH DRAW.PLOT.GRID.SETUP\n")
  
  return(gbl.var)
}

#------------------- DRAW CONNECTOR -----------------------
draw.plot.linesconnection <- function(top.vp, config.var, gbl.var) { 
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("START DRAW.PLOT.LINESCONNECTION\n")
  
  popViewport()
  # if (config.var$VERBOSE)  cat("gbl.var$r.x[1]", gbl.var$r.x[1], "\n")
  # if (config.var$VERBOSE)  cat("gbl.var$r.y[2]", gbl.var$r.y[2], "\n")
  # if (config.var$VERBOSE)  cat("gbl.var$axis.x[length(gbl.var$axis.x)]", gbl.var$axis.x[length(gbl.var$axis.x)], "\n")
  # if (config.var$VERBOSE)  cat("gbl.var$axis.y[1]", gbl.var$axis.y[1], "\n")
  #-------------------CONNECTOR VIEWPORT---------------------
  connector.vp.cormatrixmap <- viewport(xscale = c(gbl.var$r.x[1], gbl.var$axis.x[length(gbl.var$axis.x)]),
                                        yscale = c(gbl.var$axis.y[1], gbl.var$r.y[2]),
                                        height = 0.8,
                                        width = 0.8,
                                        layout.pos.row = 5,
                                        layout.pos.col = 2,
                                        name = "data.vp.cormatrixmap")
  
  connector.vp.cormatrixmap.nopval <- viewport(xscale = c(gbl.var$r.x[1], gbl.var$axis.x[length(gbl.var$axis.x)]),
                                               yscale = c(gbl.var$axis.y[1], gbl.var$r.y[2]),
                                               height = 0.8,
                                               width = 0.8,
                                               layout.pos.row = 4,
                                               layout.pos.col = 2,
                                               name = "data.vp.cormatrixmap.nopval")
  if(config.var$DISP.CORMATRIXMAP) {
    if(config.var$DISP.PVALUEPLOT) {
      pushViewport(vpTree(top.vp, vpList(connector.vp.cormatrixmap)))
    } else{
      pushViewport(vpTree(top.vp, vpList(connector.vp.cormatrixmap.nopval)))
    }
  }
  
  
  #-------------------AXIS BEGINS------------------
  # if (config.var$VERBOSE)  cat("axis\n")
  grid.lines(x = gbl.var$axis.x, y = 0, default.units = "native", gp = gpar(lty = 1, lwd = gbl.var$line.width))
  
  #----------------------- X-AXIS --------------------------
  # if (config.var$VERBOSE)  cat("x-axis\n")
  if(config.var$DISP.MULT.LAB.X) {
    #Truncate the number of X-axis labels
    if(length(gbl.var$axis.x) > 5) {
      increment <- ceiling(length(gbl.var$axis.x) / 5)
      
      axis.x.labels <- gbl.var$axis.x[seq(1, length(gbl.var$axis.x), by= increment)]
      
      if(length(gbl.var$axis.x) %% 5 != 0) {
        axis.x.labels <- c(axis.x.labels, gbl.var$axis.x[length(gbl.var$axis.x)])
      }
    }
  }
  else {
    #Truncate the number of X-axis labels
    if(length(gbl.var$axis.x) >= 2) {
      
      axis.x.labels <- c(head(gbl.var$axis.x, 1), tail(gbl.var$axis.x, 1))
    }
  }
  
  #DEBUG STATEMENT
  # if (config.var$VERBOSE)  cat("cormatrixmap\n")
  # if (config.var$VERBOSE)  cat("List",gbl.var$sorted.mydata.pos,"\n")
  # if (config.var$VERBOSE)  cat("gbl.var$sorted.mydata.pos.zoom ", gbl.var$sorted.mydata.pos.zoom, "\n")
  # if (config.var$VERBOSE)  cat("List zoom",gbl.var$sorted.mydata.pos.zoom,"\n")
  
  grid.xaxis(at = axis.x.labels, label = FALSE, gp = gpar(cex = gbl.var$cex.factor), name = "axis.x.labels")
  grid.xaxis(at = gbl.var$sorted.mydata.pos.zoom, label = FALSE, gp = gpar(cex = gbl.var$cex.factor, lwd = gbl.var$line.width))
  grid.text(axis.x.labels[1], x = -0.025, y = -0.1, just = "right", gp = gpar(fontsize = (gbl.var$font.size * 0.75)))
  grid.text(axis.x.labels[2], x = 1.015, y = -0.1, just = "left", gp = gpar(fontsize = (gbl.var$font.size * 0.75)))
  
  if(config.var$DISP.COLOR.REF) {
    #Color darkorchid1 the reference
    grid.xaxis(at = gbl.var$mydata.ref.pos, label = FALSE, gp = gpar(cex = gbl.var$cex.factor, lwd = gbl.var$line.width,col="darkorchid1")) 
  }
  
  #Remove trailing tick mark
  grid.edit("axis.x.labels", gp = gpar(col = "white"))
  
  
  
  #------------------ X-AXIS -----------------------
  #-------------------EQUIDISPOS COORDINATE BEGINS-------------
  # if (config.var$VERBOSE)  cat("equidistance\n")
  gbl.var$equidis.pos <- NULL
  
  #DEBUG STATEMENT
  # if (config.var$VERBOSE)  cat("length(gbl.var$equidis.pos) ", length(gbl.var$equidis.pos), "\n")
  # if (config.var$VERBOSE)  cat("length(gbl.var$mydata.num) ", length(gbl.var$mydata.pos), "\n")
  
  gbl.var$mydata.num <- length(gbl.var$sorted.mydata.pos)
  
  
  pos.increment <- (gbl.var$max.x - gbl.var$min.x) / gbl.var$mydata.num
  
  #Attempt to center the MYDATA labels
  start.pos <- gbl.var$min.x
  
  for(i in 0:(gbl.var$mydata.num - 1)) {
    gbl.var$equidis.pos <- c(gbl.var$equidis.pos, pos.increment * i + start.pos)
  }
  #-------------------CONNECTING LINES BEGINS-------------
  
  if(config.var$CONNECTING.LINES.ADJ == -1) {
    stop("CONNECTING.LINES.ADJ may not equal: -1")
  }
  
  tmp.correction.factor <- gbl.var$equidis.pos[1] * 10^(gbl.var$cur.exp - 9)  * (1 + config.var$CONNECTING.LINES.ADJ)
  
  #DEBUG STATEMENT
  # if (config.var$VERBOSE)  cat("tmp.correction.factor ", tmp.correction.factor, "\n")
  # if (config.var$VERBOSE)  cat("length(gbl.var$sorted.mydata.pos) ", length(gbl.var$sorted.mydata.pos), "\n")
  # if (config.var$VERBOSE)  cat("length(gbl.var$equidis.pos) ", length(gbl.var$equidis.pos), "\n")
  # if (config.var$VERBOSE)  cat("gbl.var$equidis.pos[1] ", gbl.var$equidis.pos[1], "\n")
  # if (config.var$VERBOSE)  cat("gbl.var$cur.exp ", gbl.var$cur.exp, "\n")
  # if (config.var$VERBOSE)  cat("config.var$CONNECTING.LINES.ADJ ", config.var$CONNECTING.LINES.ADJ, "\n")
  
  length.sorted.mydata.pos.zoom <- length(gbl.var$sorted.mydata.pos.zoom)
  
  if(config.var$IMAGE.SIZE == 3.5) {
    y.finish.pos <- rep(-0.9, length.sorted.mydata.pos.zoom * config.var$CONNECTING.LINES.FACTOR)
  } else if(config.var$IMAGE.SIZE == 7) {
    y.finish.pos <- rep(-1.8, length.sorted.mydata.pos.zoom * config.var$CONNECTING.LINES.FACTOR)
  }
  
  #DEBUG STATEMENT
  # if (config.var$VERBOSE)  cat("gbl.var$axis.y[1] ", gbl.var$axis.y[1], "\n")
  # if (config.var$VERBOSE)  cat("y.finish.pos ", y.finish.pos, "\n")
  
  if(config.var$CONNECTING.LINES.VERT.ADJ != -1) {
    y.start.pos <- config.var$CONNECTING.LINES.VERT.ADJ
  } else if (config.var$IMAGE.SIZE == 3.5) {
    y.start.pos <- -0.5
  } else if (config.var$IMAGE.SIZE == 7) {
    y.start.pos <- -0.7
  }else  {
    stop("Invalid image size: ", config.var$IMAGE.SIZE, "\n")
  }
  
  # if (config.var$VERBOSE)  cat("BUG HERE \n")
  
  # if (config.var$VERBOSE)  cat("BUGgg HERE \n")
  x.finish.pos <- gbl.var$equidis.pos +
    seq(1, length(gbl.var$equidis.pos) * abs(tmp.correction.factor), abs(tmp.correction.factor)) * config.var$CONNECTING.LINES.FLEX +
    gbl.var$total.dist * config.var$CONNECTING.LINES.ADJ
  
  #for reference sequence
  ref.position <- which(gbl.var$sorted.mydata.pos == gbl.var$mydata.ref.pos)
  x.finish.pos.zoom.ref <- x.finish.pos[ref.position]
  
  #for region study
  remove.position <- which(gbl.var$sorted.mydata.pos < gbl.var$min.x | gbl.var$sorted.mydata.pos > gbl.var$max.x)
  
  if(length(remove.position) > 0){
    x.finish.pos.zoom <- x.finish.pos[-remove.position]
  } else {
    x.finish.pos.zoom <- x.finish.pos
  }
  
  # if (config.var$VERBOSE)  cat("y.finish.pos ", y.finish.pos, "\n")
  # if (config.var$VERBOSE)  cat("y.start.pos ", y.start.pos, "\n")
  # if (config.var$VERBOSE)  cat("gbl.var$min.x ", gbl.var$min.x, "\n")
  # if (config.var$VERBOSE)  cat("gbl.var$max.x ", gbl.var$max.x, "\n")
  # if (config.var$VERBOSE)  cat("x.finish.pos", x.finish.pos, "\n")
  # if (config.var$VERBOSE)  cat("remove.position", remove.position, "\n")
  # if (config.var$VERBOSE)  cat("x.finish.pos.zoom ", x.finish.pos.zoom, "\n")
  # if (config.var$VERBOSE)  cat("gbl.var$sorted.mydata.pos.zoom ", gbl.var$sorted.mydata.pos.zoom, "\n")
  
  if(config.var$IMAGE.SIZE == 3.5) {
    connecting.lines <- segmentsGrob(x0 = x.finish.pos.zoom,
                                     x1 = gbl.var$sorted.mydata.pos.zoom,
                                     y0 = unit(y.finish.pos, "char"),
                                     y1 = unit(y.start.pos, "char"),
                                     default.units = "native",
                                     gp = gpar(lwd = gbl.var$line.width))
    
    ######ref
    if(config.var$DISP.COLOR.REF) {
      ######ref
      if(gbl.var$mydata.ref.pos > gbl.var$min.x & gbl.var$mydata.ref.pos < gbl.var$max.x){
        connecting.lines.ref <- segmentsGrob(x0 = x.finish.pos.zoom.ref,
                                             x1 = gbl.var$mydata.ref.pos,
                                             y0 = unit(y.finish.pos, "char"),
                                             y1 = unit(y.start.pos,"char"),
                                             default.units = "native",
                                             gp = gpar(lwd = gbl.var$line.width,col="darkorchid1"))
      }
    }
    
  } else if(config.var$IMAGE.SIZE == 7) {
    connecting.lines <- segmentsGrob(x0 = x.finish.pos.zoom,
                                     x1 = gbl.var$sorted.mydata.pos.zoom,
                                     y0 = unit(y.finish.pos, "char"),
                                     y1 = unit(y.start.pos,"char"),
                                     default.units = "native",
                                     gp = gpar(lwd = gbl.var$line.width))
    if(config.var$DISP.COLOR.REF) {
      ######ref
      if(gbl.var$mydata.ref.pos > gbl.var$min.x & gbl.var$mydata.ref.pos < gbl.var$max.x){
        connecting.lines.ref <- segmentsGrob(x0 = x.finish.pos.zoom.ref,
                                             x1 = gbl.var$mydata.ref.pos,
                                             y0 = unit(y.finish.pos, "char"),
                                             y1 = unit(y.start.pos,"char"),
                                             default.units = "native",
                                             gp = gpar(lwd = gbl.var$line.width,col="darkorchid1"))
      }
    }
  } else {
    stop("Invalid image size: ", config.var$IMAGE.SIZE, "\n")
  }
  
  
  # if (config.var$VERBOSE)  cat("BUGG HERE \n")
  #DEBUG STATEMENT
  #    if (config.var$VERBOSE)  cat("y.finish.pos ", y.finish.pos, "\n")
  #    if (config.var$VERBOSE)  cat("connecting.lines ", is.null(connecting.lines), "\n")
  #   cat(is.null(connecting.lines), "\n")
  #   cat(x.finish.pos, "\n")
  #   cat(gbl.var$sorted.mydata.pos, "\n")
  #----------------------- X-AXIS --------------------------
  grid.lines(x = gbl.var$axis.x, y = 0, default.units = "native", gp = gpar(lty = 1, lwd = gbl.var$line.width, col="gray80"))
  
  if(config.var$DISP.CONNECTING.LINES) {
    grid.draw(connecting.lines)
    if(config.var$DISP.COLOR.REF) {
      grid.draw(connecting.lines.ref)
    }
  }
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("FINISH DRAW.PLOT.LINECONNECTION\n")
}

#------------------- DRAW AXIS DATA -----------------------
draw.plot.axis.data <- function(top.vp, config.var, gbl.var) { 
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("START DRAW.PLOT.AXIS.DATA\n")
  
  popViewport()
  
  #-------------------DATA VIEWPORT---------------------
  data.vp.cormatrixmap <- viewport(xscale = c(gbl.var$r.x[1], gbl.var$axis.x[length(gbl.var$axis.x)]),
                                   yscale = c(gbl.var$axis.y[1], gbl.var$r.y[2]),
                                   height = 0.8,
                                   width = 0.8,
                                   layout.pos.row = 3,
                                   layout.pos.col = 2,
                                   name = "data.vp.cormatrixmap")
  data.vp.nocormatrixmap <- viewport(xscale = c(gbl.var$r.x[1], gbl.var$axis.x[length(gbl.var$axis.x)]),
                                     yscale = c(gbl.var$axis.y[1], gbl.var$r.y[2]),
                                     height = .8,
                                     width = .8,
                                     layout.pos.row = 3,
                                     layout.pos.col = 2,
                                     name = "data.vp.nocormatrixmap")
  if(config.var$DISP.CORMATRIXMAP) {
    pushViewport(vpTree(top.vp, vpList(data.vp.cormatrixmap)))
  }
  else {
    pushViewport(vpTree(top.vp, vpList(data.vp.nocormatrixmap)))
  }
  
  #-------------------AXIS BEGINS------------------
  
  grid.lines(x = gbl.var$axis.x, y = gbl.var$axis.y[1], default.units = "native", gp = gpar(lty = "dashed", lwd = gbl.var$line.width))
  
  
  #----------------------- X-AXIS --------------------------
  if(config.var$DISP.MULT.LAB.X) {
    #Truncate the number of X-axis labels
    if(length(gbl.var$axis.x) > 5) {
      increment <- ceiling(length(gbl.var$axis.x) / 5)
      
      axis.x.labels <- gbl.var$axis.x[seq(1, length(gbl.var$axis.x), by= increment)]
      
      if(length(gbl.var$axis.x) %% 5 != 0) {
        axis.x.labels <- c(axis.x.labels, gbl.var$axis.x[length(gbl.var$axis.x)])
      }
    }
  }
  else {
    #Truncate the number of X-axis labels
    if(length(gbl.var$axis.x) >= 2) {
      
      axis.x.labels <- c(head(gbl.var$axis.x, 1), tail(gbl.var$axis.x, 1))
    }
  }
  
  
  if(config.var$DISP.CORMATRIXMAP) {
    if(! is.null(config.var$MYDATA.LARGE.FILE)){
      grid.xaxis(at = gbl.var$sorted.mydata.large.pos.zoom, label = FALSE, gp = gpar(col="grey", cex = gbl.var$cex.factor, lwd = gbl.var$line.width))
    }
    grid.xaxis(at = axis.x.labels, label = FALSE, gp = gpar(cex = gbl.var$cex.factor), name = "axis.x.labels")
    grid.xaxis(at = gbl.var$sorted.mydata.pos.zoom, label = FALSE, gp = gpar(cex = gbl.var$cex.factor, lwd = gbl.var$line.width))
    grid.text(axis.x.labels[1], x = -0.025, y = -0.1, just = "right", gp = gpar(fontsize = (gbl.var$font.size * 0.75)))
    grid.text(axis.x.labels[2], x = 1.015, y = -0.1, just = "left", gp = gpar(fontsize = (gbl.var$font.size * 0.75)))
    
    if(config.var$DISP.COLOR.REF) {
      #Color red the reference
      grid.xaxis(at = gbl.var$mydata.ref.pos, label = FALSE, gp = gpar(cex = gbl.var$cex.factor, lwd = gbl.var$line.width, col="darkorchid1"))
    }
    
    #Remove trailing tick mark
    grid.edit("axis.x.labels", gp = gpar(col = "white"))
    
  }
  else {
    if(! is.null(config.var$MYDATA.LARGE.FILE)){
      grid.xaxis(at = gbl.var$sorted.mydata.large.pos.zoom, label = FALSE, gp = gpar(col="grey", cex = gbl.var$cex.factor, lwd = gbl.var$line.width))
    }
    
    grid.xaxis(at = axis.x.labels, label = FALSE, gp = gpar(cex = gbl.var$cex.factor), name = "axis.x.labels")
    grid.xaxis(at = gbl.var$sorted.mydata.pos.zoom, label = FALSE, gp = gpar(cex = gbl.var$cex.factor, lwd = gbl.var$line.width))
    grid.text(axis.x.labels[1], x = -0.025, y = -0.09, just = "right", gp = gpar(fontsize = (gbl.var$font.size * 0.75)))
    grid.text(axis.x.labels[2], x = 1.015, y = -0.09, just = "left", gp = gpar(fontsize = (gbl.var$font.size * 0.75)))
    
    if(config.var$DISP.COLOR.REF) {
      #Color red the reference
      grid.xaxis(at = gbl.var$mydata.ref.pos, label = FALSE, gp = gpar(cex = gbl.var$cex.factor, lwd = gbl.var$line.width, col="darkorchid1")) 
    }
    
    #Remove trailing tick mark
    grid.edit("axis.x.labels", gp = gpar(col = "white"))
    
  }
  
  #------------------ X-AXIS -----------------------
  #-------------------EQUIDISPOS COORDINATE BEGINS-------------
  
  gbl.var$equidis.pos <- NULL
  
  #DEBUG STATEMENT
  # if (config.var$VERBOSE)  cat("length(gbl.var$equidis.pos) ", length(gbl.var$equidis.pos), "\n")
  # if (config.var$VERBOSE)  cat("length(gbl.var$mydata.num) ", length(gbl.var$mydata.pos), "\n")
  
  gbl.var$mydata.num <- length(gbl.var$sorted.mydata.pos)
  
  if(config.var$DISP.CORMATRIXMAP) {
    pos.increment <- (gbl.var$max.x - gbl.var$min.x) / gbl.var$mydata.num
    
    #Attempt to center the MYDATA labels
    start.pos <- gbl.var$min.x
    
    for(i in 0:(gbl.var$mydata.num - 1)) {
      gbl.var$equidis.pos <- c(gbl.var$equidis.pos, pos.increment * i + start.pos)
    }
  } else {
    pos.increment <- (gbl.var$max.x - gbl.var$min.x) / gbl.var$mydata.num
    
    #Attempt to center the MYDATA labels
    start.pos <- gbl.var$min.x
    
    for(i in 0:(gbl.var$mydata.num - 1)) {
      gbl.var$equidis.pos <- c(gbl.var$equidis.pos, pos.increment * i + start.pos)
    }
  }
  
  #----------------- Y-AXIS --------------------
  if(config.var$LAB.Y == "ln") {
    if(config.var$DISP.MARKER.LINES) {
      
      initial.abline <- 0.1
      
      while(-log(initial.abline) < gbl.var$axis.y[1]) {
        initial.abline <- initial.abline / 10
        
        #DEBUG STATEMENT
        if (config.var$VERBOSE)  cat("initial.abline 1", initial.abline, "\n")
      }
      
      while(-log(initial.abline) < max(gbl.var$axis.y)) {
        
        #DEBUG STATEMENT
        if (config.var$VERBOSE)  cat("initial.abline 2", initial.abline, "\n")
        
        grid.lines(x = gbl.var$axis.x, y = -log(initial.abline), default.units = "native", gp = gpar(lty = "dashed", lwd = gbl.var$line.width, col="gainsboro"))
        initial.abline <- initial.abline / 10
      }
    }
  } else {
    if(config.var$DISP.MARKER.LINES) {     
      initial.abline <- 0.1
      
      while(-log10(initial.abline) < gbl.var$axis.y[1]) {
        initial.abline <- initial.abline / 10
      }
      
      while(-log10(initial.abline) < max(gbl.var$axis.y)) {
        grid.lines(x = gbl.var$axis.x, y = -log10(initial.abline), default.units = "native", gp = gpar(lty = "dashed", lwd = gbl.var$line.width, col="gainsboro"))
        initial.abline <- initial.abline / 10
      }
    }
  }
  
  #----------- X-AXIS of THRESHOLD
  if(config.var$PVAL.THRESHOLD != 0 & config.var$PVAL.THRESHOLD > config.var$DISP.PVAL.THRESHOLD){
    grid.lines(x = gbl.var$axis.x, y = config.var$PVAL.THRESHOLD, default.units = "native", gp = gpar(lty = "dashed", lwd = gbl.var$line.width , col="red"))
    
  }
  #----------- X-AXIS of base
  grid.lines(x = gbl.var$axis.x, y = gbl.var$axis.y[1], default.units = "native", gp = gpar(lwd = gbl.var$line.width))
  
  grid.yaxis(at = gbl.var$axis.y, label = gbl.var$axis.y, gp = gpar(fontsize = gbl.var$font.size, lwd = gbl.var$line.width))
  #----------- WRITE LABEL for Y-AXIS
  if(config.var$LAB.Y == "ln") {
    grid.text("-ln(p-value)",
              x = -0.1,
              y = 0.5,
              rot = 90,
              gp = gpar(fontsize = gbl.var$font.size,
                        fontface = "bold"))
  } else {
    grid.text("-log(p-value)",
              x = -0.1,
              y = 0.5,
              rot = 90,
              gp = gpar(fontsize = gbl.var$font.size,
                        fontface = "bold"))
  }
  
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("FINISH DRAW.PLOT.AXIS.DATA\n")
}


#------------------- DRAW ANNOTATION -----------------------
draw.plot.annotation <- function(config.var, gbl.var) { 
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("START DRAW.PLOT.ANNOTATION\n")
  
  popViewport()
  top.vp <- gbl.var$top.vp
  annotation.vp.cormatrixmap <- viewport(height = 0.95,
                                         width = 1,
                                         layout.pos.row = 4,
                                         layout.pos.col = 2,
                                         name = "annotation.vp.cormatrixmap")
  annotation.vp.nocormatrixmap <- viewport(height = 0.95,
                                           width = 1,
                                           layout.pos.row = 5,
                                           layout.pos.col = 2,
                                           name = "annotation.vp.nocormatrixmap")
  annotation.vp.cormatrixmap.nopval <- viewport(height = 0.95,
                                                width = 1,
                                                layout.pos.row = 3,
                                                layout.pos.col = 2,
                                                name = "annotation.vp.cormatrixmap.nopval")
  
  if(config.var$DISP.CORMATRIXMAP) {
    if(config.var$DISP.PVALUEPLOT) {
      pushViewport(vpTree(top.vp, vpList(annotation.vp.cormatrixmap)))
      pushViewport(annotation.vp.cormatrixmap)
    } else {
      pushViewport(vpTree(top.vp, vpList(annotation.vp.cormatrixmap.nopval)))
      pushViewport(annotation.vp.cormatrixmap.nopval)
    }
  }
  else {
    pushViewport(vpTree(top.vp, vpList(annotation.vp.nocormatrixmap)))
    pushViewport(annotation.vp.nocormatrixmap)
  }
  
  
  size_gviz <- 0
  if( ! is.null(gbl.var$listtracks_gviz)){
    size_gviz <- length(gbl.var$listtracks_gviz) * 2
  }
  
  size_gviz_user <- 0
  if(!is.null(config.var$BIOFEAT.USER.FILE)){
    #   gbl.var <- create.tracks.user(config.var,gbl.var)
    size_gviz_user <- length(gbl.var$listtracks_user) 
  }
  
  size_trackviewer <- 0
  if( ! is.null(gbl.var$listtracks_trackviewer)){
    size_trackviewer <- length(gbl.var$listtracks_trackviewer) 
  }
  size_ggbio <- 0
  if( ! is.null(gbl.var$listtracks_ggbio)){
    size_ggbio <- length(gbl.var$listtracks_ggbio) * 4
  }
  total_tracks <-  size_gviz + size_trackviewer + size_ggbio
  
  if (config.var$VERBOSE)  cat("Number of track", total_tracks, "GVIZ",  size_gviz,"\n")
  
  annotvp <-viewport(layout=grid.layout(4, 1, heights=unit(c(size_gviz,size_ggbio,size_trackviewer,size_gviz_user),"null")),name = "annotviewport")
  gvizviewport <- viewport(layout.pos.col = 1, layout.pos.row = 1, name = "gvizviewport")
  ggbioviewport <- viewport(layout.pos.col = 1, layout.pos.row = 2, name = "ggbioviewport")
  trackviewerviewport <- viewport(layout.pos.col = 1, layout.pos.row = 3, name = "trackviewerviewport")
  gvizuserviewport <- viewport(layout.pos.col = 1, layout.pos.row = 4, name = "gvizuserviewport")
  
  if(config.var$DISP.CORMATRIXMAP) {
    if(config.var$DISP.PVALUEPLOT) {
      pushViewport(vpTree(annotation.vp.cormatrixmap, vpList(gvizviewport, ggbioviewport, trackviewerviewport,gvizuserviewport)))
    } else {
      pushViewport(vpTree(annotation.vp.cormatrixmap.nopval, vpList(gvizviewport, ggbioviewport, trackviewerviewport,gvizuserviewport)))
    }
  }
  else {
    pushViewport(vpTree(annotation.vp.nocormatrixmap, vpList(gvizviewport, ggbioviewport, trackviewerviewport,gvizuserviewport)))
  }
  pushViewport(annotvp)
  
  if(! (is.null(gbl.var$listtracks_gviz))){
    seekViewport("annotviewport")
    pushViewport(gvizviewport)
    plotTracks(gbl.var$listtracks_gviz, from=gbl.var$min.x, to=gbl.var$max.x,panel.only=TRUE,add=TRUE, fontsize=5)
    
  }
  
  if(! (is.null(gbl.var$listtracks_ggbio))){
    seekViewport("annotviewport")
    pushViewport(ggbioviewport)
    l.g <- lapply(gbl.var$listtracks_ggbio, function(x){
      ggplotGrob(x)
    })
    # grid.rect(gp=gpar(col="blue"))
    grid.arrange(do.call(arrangeGrob, c(l.g, list(nrow = size_ggbio/4, ncol = 1))),newpage=FALSE)
    #tracks(gbl.var$listtracks_ggbio)
  }
  
  if(! (is.null(gbl.var$listtracks_trackviewer))){
    seekViewport("annotviewport")
    pushViewport(trackviewerviewport)
    viewTracks(trackList(gbl.var$listtracks_trackviewer),chromosome=gbl.var$mydata.chr,start=gbl.var$min.x, end=gbl.var$max.x,newpage=FALSE)
  }
  
  if(! (is.null(gbl.var$listtracks_user))){
    seekViewport("annotviewport")
    pushViewport(gvizuserviewport)
    plotTracks(gbl.var$listtracks_user, from=gbl.var$min.x, to=gbl.var$max.x,panel.only=TRUE,add=TRUE, fontsize=5)
  }
  
  popViewport()
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("FINISH DRAW.PLOT.ANNOTATION\n")
}


#------------------- DRAW DATA via GGBIO -----------------------
draw.plot.mydata.ggbio <- function(config.var, gbl.var,numfile) { 
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("START DRAW.PLOT.MYDATA\n")
  split.mydata.file <- gbl.var$split.mydata.file
  gbl.var <-read.file.mydata(split.mydata.file[[1]],config.var, gbl.var,numfile)
  # if (config.var$VERBOSE)  cat("READ ",i," MYDATA file \n")
  mydata.data <- gbl.var$mydata.data
  mydata.pv <- log(mydata.data[,4])
  startlist<-mydata.data[,3]
  endlist<-mydata.data[,3]
  id<-mydata.data[,1]
  cprange<-GRanges(gbl.var$mydata.chr,IRanges(startlist,endlist))
  #  startlarge<-start-50000
  #  endlarge<-end+50000
  
  cpidrange<-GRanges(gbl.var$mydata.chr,IRanges(startlist,endlist),"*",id)
  data<-DataTrack(range=cprange,start,end,data=mydata.pv,name="CpG pvalue",cex=.16)
  namedata_track <- AnnotationTrack(chromosome=gbl.var$mydata.chr,strand ="*",start=startlist,end=startlist,
                                    feature="CpG",group=id,id=id, name = "feature")
  
  return(list(data,namedata_track))
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("FINISH DRAW.PLOT.MYDATA\n")
}

#-------------------DRAW SINGLE MYDATA------------------
draw.plot.grid.mydata <- function(config.var, gbl.var) {
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("START DRAW.PLOT.GRID.MYDATA\n")
  
  # Grab from configuration is available
  if(is.null(config.var$SYMBOL.FACTOR)) {
    tmp.cex.factor.symbol = gbl.var$cex.factor.symbol
  } else {
    tmp.cex.factor.symbol = config.var$SYMBOL.FACTOR
  }
  
  
  list.colors.bar <- create.color.bar(config.var, gbl.var)
  color.cut <- list.colors.bar$color.cut
  color.cut.ref <- list.colors.bar$color.cut.ref
  #z<-color.cut[,mydata.num]
  #  gbl.var$fill.list[gbl.var$cur.sample]
  format <- as.character(gbl.var$split.format[[1]][gbl.var$cur.sample])
  if(config.var$DISP.MYDATA) {
    if(config.var$DISP.TYPE == "symbol") {
      if(gbl.var$cur.sample == 1) {
        
        #------------- DATA of INTERET
        for (i in 1:nrow(gbl.var$mydata.data)) {
          name.test <-gbl.var$mydata.data$MYDATA.NAME[i]
          position <- gbl.var$mydata.hash.names.pos[[name.test ]]
          #cat ("color",i, " ", color.cut.ref[i],"\n")
          if ((gbl.var$mydata.data$MYDATA.PVAL[i] > config.var$DISP.PVAL.THRESHOLD) & 
                (gbl.var$mydata.data$MYDATA.PVAL[i] != 0) &
                (position > (gbl.var$min.x - 1)) &
                (position < (gbl.var$max.x + 1)) ) {
            
            #---- ASSOCIATION COLOR
            mycolor <- color.cut.ref[i]
            if(! is.null(config.var$DISP.ASSOCIATION)){
              if(as.logical(gbl.var$split.association[[1]][gbl.var$cur.sample]) & (grepl("ASSO", format)[1])) {
                mycolor <- gbl.var$fill.list[gbl.var$cur.sample]
                mycolor <- gbl.var$fill.list[gbl.var$cur.sample]
                if  (exists("gbl.var$mydata.data$MYDATA.ASSO[i]") & length(gbl.var$mydata.data$MYDATA.ASSO[i]) >0 & gbl.var$mydata.data$MYDATA.ASSO[i] == "-") {
                  mycolor.tmp <- opposite(mycolor,plot=FALSE)[2]
                  mycolor <- mycolor.tmp
                } else if  (exists("gbl.var$mydata.data$MYDATA.ASSO[i]") & length(gbl.var$mydata.data$MYDATA.ASSO[i]) >0 & gbl.var$mydata.data$MYDATA.ASSO[i] == "0") {
                  mycolor.tmp <- "black"
                  mycolor <- mycolor.tmp
                }
              }
            }
            
            
            # if (config.var$VERBOSE)  cat("config.var$MYDATA.REF",config.var$MYDATA.REF,"\n")
            # if (config.var$VERBOSE)  cat("gbl.var$mydata.data$MYDATA.NAME[i]",gbl.var$mydata.data$MYDATA.NAME[i],"\n")
            # if (config.var$VERBOSE)  cat("My data",i,"vs my ref", gbl.var$mydata.best.position,"\n")
            ref.print <- "FALSE"
            
            if (exists("config.var$MYDATA.REF") & !is.null(config.var$MYDATA.REF) ) {
              if (length(config.var$MYDATA.REF) > 0 & (gbl.var$mydata.data$MYDATA.NAME[i] == config.var$MYDATA.REF)){
                ref.print <- "TRUE"
              }
            } else if ( (length(config.var$MYDATA.REF) > 0 ) & (i == gbl.var$mydata.best.position)){
              ref.print <- "TRUE"	
            }
            
            if (i == gbl.var$mydata.best.position){
              ref.print <- "TRUE"	
            }
            
            if (as.logical(ref.print)) {
              
              grid.points(position,
                          gbl.var$mydata.data$MYDATA.PVAL[i],
                          pch = gbl.var$symbol.list[gbl.var$cur.sample],
                          gp = gpar(col = mycolor,
                                    cex = tmp.cex.factor.symbol,
                                    fill = "black" ,
                                    lwd = gbl.var$line.width))
              
              if(!is.null(config.var$DISP.REGION)){
                if(as.logical(gbl.var$split.region[[1]][gbl.var$cur.sample])  & (grepl("REGION", format)[1])){
                  #  if (config.var$VERBOSE)  cat("format",format,"\n")
                  position.start <- gbl.var$mydata.hash.names.start[[name.test ]]
                  position.end <- gbl.var$mydata.hash.names.end[[name.test ]]
                  if(position.start < (gbl.var$min.x - 1)){
                    position.start <- (gbl.var$min.x - 1)
                  }
                  if(position.end  > (gbl.var$max.x + 1)){
                    position.end <- (gbl.var$max.x + 1)
                  }
                  grid.lines(x = c(position.start,position.end),
                             y = c(gbl.var$mydata.data$MYDATA.PVAL[i],gbl.var$mydata.data$MYDATA.PVAL[i]), 
                             default.units = "native", gp = gpar(lty = "solid", lwd = gbl.var$line.width, col=mycolor))
                }
              }
            }else {
              grid.points(position,
                          gbl.var$mydata.data$MYDATA.PVAL[i],
                          pch = gbl.var$symbol.list[gbl.var$cur.sample],
                          gp = gpar(col = "black",
                                    cex = tmp.cex.factor.symbol,
                                    fill = mycolor ,
                                    lwd = gbl.var$line.width))
              
              if(!is.null(config.var$DISP.REGION)){
                if(as.logical(gbl.var$split.region[[1]][gbl.var$cur.sample])  & (grepl("REGION", format)[1] )){
                  position.start <- gbl.var$mydata.hash.names.start[[name.test ]]
                  position.end <- gbl.var$mydata.hash.names.end[[name.test ]]
                  if(position.start < (gbl.var$min.x - 1)){
                    position.start <- (gbl.var$min.x - 1)
                  }
                  if(position.end > (gbl.var$max.x + 1)){
                    position.end <- (gbl.var$max.x + 1)
                  }
                  grid.lines(x = c(position.start,position.end),
                             y = c(gbl.var$mydata.data$MYDATA.PVAL[i],gbl.var$mydata.data$MYDATA.PVAL[i]), 
                             default.units = "native", gp = gpar(lty = "solid", lwd = gbl.var$line.width, col=mycolor))
                }
              }
              
            }
            
          }
        }
      } else {
        for (i in 1:nrow(gbl.var$mydata.data)) {
          name.test <-gbl.var$mydata.data$MYDATA.NAME[i]
          position <- gbl.var$mydata.hash.names.pos[[name.test ]]
          if ((gbl.var$mydata.data$MYDATA.PVAL[i] > config.var$DISP.PVAL.THRESHOLD) & 
                (gbl.var$mydata.data$MYDATA.PVAL[i] != 0) &
                (position > (gbl.var$min.x - 1)) &
                (position < (gbl.var$max.x + 1)) ) {
            
            #---- ASSOCIATION COLOR
            mycolor <- color.cut.ref[i]
            if(! is.null(config.var$DISP.ASSOCIATION)){
              if(gbl.var$split.association[[1]][gbl.var$cur.sample] & (grepl("ASSO", format)[1]) & (grepl("ASSO", format)[1])) {
                mycolor <- gbl.var$fill.list[gbl.var$cur.sample]
                if  (gbl.var$mydata.data$MYDATA.ASSO[i] == "-") {
                  mycolor.tmp <- opposite(mycolor,plot=FALSE)[2]
                  mycolor <- mycolor.tmp
                } else if  (gbl.var$mydata.data$MYDATA.ASSO[i] == "0") {
                  mycolor.tmp <- "black"
                  mycolor <- mycolor.tmp
                }
              }
            }
            
            
            if (is.na(gbl.var$symbol.list[gbl.var$cur.sample])) {
              grid.points(position, 
                          gbl.var$mydata.data$MYDATA.PVAL[i], 
                          pch = 25, 
                          gp = gpar(col = gbl.var$color.list[gbl.var$cur.sample], 
                                    cex = tmp.cex.factor.symbol,
                                    fill = mycolor, 
                                    lwd = gbl.var$line.width))
              
            }
            else { 
              grid.points(position, 
                          gbl.var$mydata.data$MYDATA.PVAL[i], 
                          pch = gbl.var$symbol.list[gbl.var$cur.sample],
                          gp = gpar(col = gbl.var$color.list[gbl.var$cur.sample], 
                                    cex = tmp.cex.factor.symbol, 
                                    fill = mycolor, 
                                    lwd = gbl.var$line.width))
            }
            
            if(! is.null(config.var$DISP.REGION)){
              if(as.logical(gbl.var$split.region[[1]][gbl.var$cur.sample])  & (grepl("REGION", format)[1] )){
                position.start <- gbl.var$mydata.hash.names.start[[name.test ]]
                position.end <- gbl.var$mydata.hash.names.end[[name.test ]]
                if(position.start < (gbl.var$min.x - 1)){
                  position.start <- (gbl.var$min.x - 1)
                }
                if(position.end  > (gbl.var$max.x + 1)){
                  position.end <- (gbl.var$max.x + 1)
                }
                grid.lines(x = c(position.start,position.end),
                           y = c(gbl.var$mydata.data$MYDATA.PVAL[i],gbl.var$mydata.data$MYDATA.PVAL[i]), 
                           default.units = "native", gp = gpar(lty = "solid", lwd = gbl.var$line.width, col=mycolor))
              }
            }
            
          }
        }
      }
    }
    else {
      stop("Invalid display type: ", config.var$DISP.TYPE, "\n")
    }
  }
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("FINISH DRAW.PLOT.GRID.MYDATA\n")
}

#-------------------DRAW LARGE MYDATA------------------
draw.plot.grid.mydata.large <- function(config.var, gbl.var) {
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("START DRAW.PLOT.GRID.MYDATA.LARGE\n")
  
  # Grab from configuration is available
  if(is.null(config.var$SYMBOL.FACTOR)) {
    tmp.cex.factor.symbol = gbl.var$cex.factor.symbol
  } else {
    tmp.cex.factor.symbol = config.var$SYMBOL.FACTOR
  }
  format <- as.character(gbl.var$large.split.format[[1]][gbl.var$cur.sample.large])
  if(config.var$DISP.MYDATA) {
    if(config.var$DISP.TYPE == "symbol") {
      #------------- DATA LARGE
      if(! is.null(config.var$MYDATA.LARGE.FILE)) {
        for (i in 1:nrow(gbl.var$mydata.large.data)) {
          
          #---- ASSOCIATION COLOR
          mycolor <- gbl.var$large.fill.list[gbl.var$cur.sample.large]
          
          if(!is.null(config.var$DISP.ASSOCIATION.LARGE)){
            if(as.logical(gbl.var$large.split.association[[1]][gbl.var$cur.sample.large]) & (grepl("ASSO", format)[1])) {
              # if (config.var$VERBOSE)  cat("Association ",gbl.var$mydata.large.data$MYDATA.ASSO[i]," \n") 
              if  (gbl.var$mydata.large.data$MYDATA.ASSO[i] == "-") {
                mycolor.fill.tmp <- opposite(mycolor,plot=FALSE)[2]
                mycolor <- mycolor.fill.tmp
              } else if  (gbl.var$mydata.large.data$MYDATA.ASSO[i] == "0") {
                mycolor.fill.tmp <- "black"
                mycolor <- mycolor.fill.tmp
              }
            }
          }
          
          name.test <-gbl.var$mydata.large.data$MYDATA.NAME[i]
          position <- gbl.var$mydata.large.hash.names.pos[[name.test]]
          # if (config.var$VERBOSE)  cat("position",position)
          if ((gbl.var$mydata.large.data$MYDATA.PVAL[i] > config.var$DISP.PVAL.THRESHOLD) & 
                (gbl.var$mydata.large.data$MYDATA.PVAL[i] != 0) &
                (position > (gbl.var$min.x - 1)) &
                (position < (gbl.var$max.x + 1)) & 
                (gbl.var$mydata.large.data$MYDATA.PVAL[i] > (gbl.var$min.y - 1)) &
                (gbl.var$mydata.large.data$MYDATA.PVAL[i] < (gbl.var$max.y + 1)) ) {
            
            # if (config.var$VERBOSE)  cat("config.var$MYDATA.REF",config.var$MYDATA.REF,"\n")
            # if (config.var$VERBOSE)  cat("gbl.var$mydata.data$MYDATA.NAME[i]",gbl.var$mydata.data$MYDATA.NAME[i],"\n")
            if (is.na(gbl.var$large.symbol.list[gbl.var$cur.sample.large])) {
              grid.points(position,
                          gbl.var$mydata.large.data$MYDATA.PVAL[i],
                          pch = 21,
                          gp = gpar(col = gbl.var$large.color.list[gbl.var$cur.sample.large],
                                    cex = tmp.cex.factor.symbol,
                                    fill = mycolor,
                                    lwd = gbl.var$line.width))
            } else {
              grid.points(position,
                          gbl.var$mydata.large.data$MYDATA.PVAL[i],
                          pch = gbl.var$large.symbol.list[gbl.var$cur.sample.large],
                          gp = gpar(col = gbl.var$large.color.list[gbl.var$cur.sample.large],
                                    cex = tmp.cex.factor.symbol,
                                    fill = mycolor,
                                    lwd = gbl.var$line.width))
            }
            #  if (config.var$VERBOSE)  cat("format ",format,"\n")
            #  if (config.var$VERBOSE)  cat("Sample",gbl.var$cur.sample.large,"\n")
            
            # if(grepl("REGION", format)[1] ){
            #    if (config.var$VERBOSE)  cat("format ",format,"\n")
            # }
            if (grepl("REGION", format)[1]){
              
              if(! is.null(config.var$DISP.REGION.LARGE))   {
                region.print <- "FALSE"
                if(! is.null(config.var$DISP.REGION.LARGE)){
                  # if(exists("gbl.var$large.split.region[[1]][gbl.var$cur.sample.large])"){
                  #if(!is.null(gbl.var$large.split.region[[1]][gbl.var$cur.sample.large])){
                  if(as.logical(gbl.var$large.split.region[[1]][gbl.var$cur.sample.large])){
                    region.print <- "TRUE"
                  }
                  
                  if(as.logical(region.print) ){
                    if (config.var$VERBOSE)  cat("Region",gbl.var$large.split.region[[1]][gbl.var$cur.sample.large],"\n")
                    position.start <- gbl.var$mydata.large.hash.names.start[[name.test]]
                    position.end <- gbl.var$mydata.large.hash.names.end[[name.test]]
                    if(position.start < (gbl.var$min.x - 1)){
                      position.start <- (gbl.var$min.x - 1)
                    }
                    if(position.end  > (gbl.var$max.x + 1)){
                      position.end <- (gbl.var$max.x + 1)
                    }
                    if (config.var$VERBOSE)  cat("position.start ",position.start,"\n")
                    if (config.var$VERBOSE)  cat("position.end ",position.end,"\n")
                    if (config.var$VERBOSE)  cat("y",gbl.var$mydata.large.data$MYDATA.PVAL[i],"\n")
                    if (config.var$VERBOSE)  cat("axis.x ", gbl.var$axis.x,"\n")
                    grid.lines(x = c(position.start,position.end),
                               y = c(gbl.var$mydata.large.data$MYDATA.PVAL[i],gbl.var$mydata.large.data$MYDATA.PVAL[i]), 
                               default.units = "native", gp = gpar(lty = "solid", lwd = gbl.var$line.width, col=mycolor)) 
                  }
                } 
              }
            }
          }
        }
      }
    }
    else {
      stop("Invalid display type: ", config.var$DISP.TYPE, "\n")
    }
  }
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("FINISH DRAW.PLOT.GRID.MYDATA.LARGE\n")
}



#-------------------DRAW CORMATRIX plot------------------

draw.plot.cormatrix.plot <- function(config.var, gbl.var) {
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("START DRAW.PLOT.CORMATRIX.PLOT\n")
  
  cormatrix.matrix <- gbl.var$cormatrix.data
  
  matrix.rows <- dim(cormatrix.matrix)[1]
  matrix.cols <- dim(cormatrix.matrix)[2]
  
  list.colors.bar <- create.color.bar(config.var, gbl.var)
  color.cut <- list.colors.bar$color.cut
  cormatrix.key <- list.colors.bar$cormatrix.key
  map.label.ldtype <- list.colors.bar$map.label.ldtype
  map.label.distance <- list.colors.bar$map.label.distance
  #cat ("color",1, " ", color.cut)
  
  x.x <- (1:matrix.cols)/matrix.cols
  y.y <- (1:matrix.rows)/matrix.rows
  right <- rep(x.x, nrow(gbl.var$general.data))
  top <- rep(y.y, each=matrix.cols)
  
  image.rect <- rectGrob(x=right,
                         y=top,
                         width=1/matrix.cols,
                         height=1/matrix.rows,
                         just=c("right", "top"),
                         gp=gpar(col=NULL,
                                 fill=color.cut,
                                 lty="blank"),
                         name="image.rect")
  
  #seq.x <- c(1*(1/gbl.var$mydata.num), 1)
  #seq.y <- c(0, (1 - (1/gbl.var$mydata.num)))
  
  tmp.mydata.num <- gbl.var$mydata.num
  
  seq.x.names <- seq((1/tmp.mydata.num), 1, by=(1/tmp.mydata.num))
  seq.y.names <- seq(0, (1 - (1/tmp.mydata.num)), by=(1/tmp.mydata.num))
  
  #The condition for the tmp.mydata.num is met in the original construction of the sequence
  if(tmp.mydata.num > 50) {
    tmp.x.factor <- seq.x.names[1]
    
    seq.x.names <- seq.x.names + tmp.x.factor/2
    seq.y.names <- seq.y.names - tmp.x.factor/2
  } else if(tmp.mydata.num < 50) {
    tmp.x.factor <- (seq.x.names[1] - 0.02)/2 + 0.02
    tmp.y.factor <- (seq.x.names[1] - 0.02)/2
    
    seq.x.names <- seq(tmp.x.factor, (1 - tmp.y.factor), by=(1/tmp.mydata.num))
    seq.y.names <- seq(tmp.y.factor, (1 - tmp.y.factor), by=(1/tmp.mydata.num))
  }
  
  mydata.names.pos.x <- (1:gbl.var$mydata.num)/gbl.var$mydata.num
  mydata.names.pos.y <- ((1:gbl.var$mydata.num)/gbl.var$mydata.num - 1/gbl.var$mydata.num)
  
  if(is.null(config.var$FONT.FACTOR)) {
    tmp.font.factor = (gbl.var$cex.factor - 0.5)
  } else {
    tmp.font.factor = config.var$FONT.FACTOR
  }
  
  if (config.var$DISP.MYDATA.NAMES) {
    if( length(gbl.var$sorted.mydata.names) < 30) {
      list.name <- paste("   ",gbl.var$sorted.mydata.names)
      name.ref <- paste("   ",gbl.var$sorted.mydata.names[gbl.var$mydata.ref.pos])
    } else {
      list.name <- paste("",gbl.var$sorted.mydata.names)
      name.ref <- paste("",gbl.var$sorted.mydata.names[gbl.var$mydata.ref.pos])
    }
    
    rot.name <- -45
    #if (config.var$CORMATRIX.FORMAT == "RAW") {
    # rot.name <- 135
    #}
    mydata.names.cormatrix.sec <- textGrob(rev(list.name),
                                           x = seq.x.names[1:gbl.var$mydata.num],
                                           y = seq.y.names[1:gbl.var$mydata.num],
                                           rot=rot.name,
                                           gp=gpar(cex = tmp.font.factor),
                                           just =  c("left"),
                                           name="mydata.names.cormatrix.sec")
    if(config.var$DISP.COLOR.REF) {
      #REFERENCE
      mydata.names.cormatrix.sec.ref <- textGrob(rev(name.ref),
                                                 x = seq.x.names[gbl.var$mydata.ref.pos],
                                                 y = seq.y.names[gbl.var$mydata.ref.pos],
                                                 rot=rot.name,
                                                 gp=gpar(cex = tmp.font.factor,col="darkorchid1"),
                                                 just =  c("left"),
                                                 name="mydata.names.cormatrix.sec")
    }
    
  }
  
  popViewport()
  
  #67.5
  ##### VIEW PORT FOR LD MAP
  cormatrix.map.vp <- viewport(y = unit(0.415, "npc"),
                               x = unit(0.471875, "npc"),
                               width = unit(0.7, "snpc"),
                               height = unit(0.7, "snpc"),
                               angle = 67.5,
                               just = c("center", "center"),
                               name = "cormatrix.map.vp")
  if (config.var$DISP.MYDATA.NAMES) {
    if(config.var$DISP.COLOR.REF) {
      cormatrix.map <- gTree(children=gList(image.rect, mydata.names.cormatrix.sec, mydata.names.cormatrix.sec.ref), just=c("center", "bottom"), vp=cormatrix.map.vp, name="cormatrix.map")
    } else {
      cormatrix.map <- gTree(children=gList(image.rect, mydata.names.cormatrix.sec), just=c("center", "bottom"), vp=cormatrix.map.vp, name="cormatrix.map")
    }
  } else {
    cormatrix.map <- gTree(children=gList(image.rect), just=c("center"), vp=cormatrix.map.vp, name="cormatrix.map")
  }
  
  gene.map.vp <- viewport(name = "gene.map.vp")
  gene.map <- gTree(children=gList(map.label.ldtype, map.label.distance), vp = gene.map.vp, name="gene.map")
  
  pushViewport(cormatrix.map.vp)
  
  grid.draw(cormatrix.map)
  
  popViewport()
  
  pushViewport(gene.map.vp)
  
  grid.draw(gene.map)
  
  if(config.var$DISP.COLOR.BAR) {
    #DEBUG STATEMENT
    # if (config.var$VERBOSE)  cat("is.null(cormatrix.key) ", is.null(cormatrix.key), "\n")
    
    grid.draw(cormatrix.key)
  }
  
  popViewport()
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("FINISH DRAW.PLOT.CORMATRIX.PLOT\n")
}

#-------------------DRAW LEGEND----------------------------
draw.legend <- function(config.var, gbl.var) {
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("START DRAW.LEGEND\n")
  
  #DEBUG STATEMENT
  # if (config.var$VERBOSE)  cat("gbl.var$font.size ", gbl.var$font.size, "\n")
  if(config.var$DISP.CORMATRIXMAP) {
    legend.vp <- viewport(x = 0.5, y = -2, height = 1, width = 1, name = "legend.vp")
  }
  else {
    legend.vp<- viewport(x = 0.5, y = 0.05, height = 1, width = 1, name = "legend.vp")
  }
  
  pushViewport(legend.vp)
  
  
  #-------------------LEGEND SAMPLES BEGINS-----------------------
  
  legend.samples.list <- NULL
  legend.samples <- gTree(children=NULL,
                          just=c("center", "bottom"), vp=legend.vp,
                          name="legend.samples")
  
  samples <- gbl.var$samples
  
  if(samples <= 8) {
    y.pos.sec.row.correction <- 0
    x.pos.sec.row.correction <- 0
    
    
    
    #--------------- LEGEND FOR REFERENCE
    #  if (config.var$VERBOSE)  cat("Reference \n")
    
    #Ensure that labels do not overlap by moving them
    if(config.var$DISP.CORMATRIXMAP) {
      y.pos <- .8 
      x.pos <- 0.75 
      
    } else {
      y.pos <- -0.3125 + y.pos.sec.row.correction
      x.pos <- -0.19 + x.pos.sec.row.correction
    }
    
    label.text.ref <- gbl.var$mydata.data$MYDATA.NAME[which(gbl.var$mydata.data$LOC == gbl.var$mydata.ref.pos)]
    label.sample.text <- textGrob(label.text.ref,
                                  x = (x.pos + 0.035),
                                  y = y.pos,
                                  just=("left"),
                                  gp = gpar(fontsize = gbl.var$font.size *0.75))
    
    
    #FOR LD adjust is .3
    #FOR NO-LD adjust is
    label.sample.symbol <- pointsGrob(x = x.pos,
                                      y = y.pos,
                                      pch = 21,
                                      gp = gpar(col = "red",
                                                cex = (gbl.var$cex.factor - 0.5),
                                                lwd = gbl.var$line.width,
                                                fill = "black"))
    
    legend.samples.list[[(length(legend.samples.list) + 1)]] <- label.sample.text
    legend.samples.list[[(length(legend.samples.list) + 1)]] <- label.sample.symbol
    
    #All samples
    # if (config.var$VERBOSE)  cat("All samples \n")
    for(i in 1:samples) {
      if(i > 3) {
        x.pos.sec.row.correction <- -1.5
        y.pos.sec.row.correction <- -.04
      }
      
      #label.text <- substr(gbl.var$split.sample.labels[[1]][i],1, 30)
      label.text <- gbl.var$split.sample.labels[[1]][i]
      
      #Ensure that labels do not overlap by moving them
      if(config.var$DISP.CORMATRIXMAP) {
        y.pos <- .8 + 0.1 * i  
        x.pos <- 0.75
        
      } else {
        y.pos <- -0.3125 + y.pos.sec.row.correction
        x.pos <- -0.19 + 0.35* i  + x.pos.sec.row.correction
      }
      
      label.sample.text <- textGrob(label.text,
                                    x = (x.pos + 0.035),
                                    y = y.pos,
                                    just=("left"),
                                    gp = gpar(fontsize = (gbl.var$font.size)*0.75 ))
      
      if(is.na(gbl.var$symbol.list[i])) {
        gbl.var$symbol.list[i] <- "circle-fill"
      }
      
      #FOR LD adjust is .3
      #FOR NO-LD adjust is
      label.sample.symbol <- pointsGrob(x = x.pos,
                                        y = y.pos,
                                        pch = gbl.var$symbol.list[i],
                                        gp = gpar(col = "black",
                                                  cex = (gbl.var$cex.factor - 0.5),
                                                  lwd = gbl.var$line.width,
                                                  fill = gbl.var$fill.list[i]))
      
      legend.samples.list[[(length(legend.samples.list) + 1)]] <- label.sample.text
      legend.samples.list[[(length(legend.samples.list) + 1)]] <- label.sample.symbol
    }
    
    #--------------- LEGEND LARGE DATA
    if(! is.null(config.var$MYDATA.LARGE.FILE)){
      large.samples <- gbl.var$large.samples
      for(i in 1:large.samples) {
        if((samples + i) > 3) {
          x.pos.sec.row.correction <- -1.5
          y.pos.sec.row.correction <- -.04
        }
        
        #label.text <- substr(gbl.var$large.split.sample.labels[[1]][i], 1, 30)
        label.text <- gbl.var$large.split.sample.labels[[1]][i]
        
        #Ensure that labels do not overlap by moving them
        if(config.var$DISP.CORMATRIXMAP) {
          y.pos <- .8 + 0.1 * (samples + i ) 
          x.pos <- 0.75
          
        } else {
          y.pos <- -0.3125 + y.pos.sec.row.correction
          x.pos <- -0.19 + 0.35* ( samples + i ) + x.pos.sec.row.correction
        }
        
        label.sample.text <- textGrob(label.text,
                                      x = (x.pos + 0.035),
                                      y = y.pos,
                                      just=("left"),
                                      gp = gpar(fontsize = (gbl.var$font.size) *0.75 ))
        
        if(is.na(gbl.var$large.symbol.list[i])) {
          gbl.var$large.symbol.list[i] <- 24
        }
        
        #FOR LD adjust is .3
        #FOR NO-LD adjust is
        label.sample.symbol <- pointsGrob(x = x.pos,
                                          y = y.pos,
                                          pch = gbl.var$large.symbol.list[i],
                                          gp = gpar(col = gbl.var$large.color.list[i],
                                                    cex = (gbl.var$cex.factor - 0.5),
                                                    lwd = gbl.var$line.width,
                                                    fill = gbl.var$large.fill.list[i]))
        
        legend.samples.list[[(length(legend.samples.list) + 1)]] <- label.sample.text
        legend.samples.list[[(length(legend.samples.list) + 1)]] <- label.sample.symbol
      }
    }
    
    class(legend.samples.list) <- c("gList")
    
    legend.glist.samples <- legend.samples.list
    
    legend.samples <- gTree(children=legend.glist.samples,
                            just=c("left", "bottom"),
                            vp=legend.vp,
                            name="legend.samples")
    
    grid.draw(legend.samples)
  }
  
  popViewport()
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("FINISH DRAW.LEGEND\n")
}

#-------------------DRAW MYDATA NAME----------------------------

#draw MYDATA labels
draw.plot.grid.mydata.names <- function(config.var, gbl.var) {
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("START DRAW.PLOT.GRID.MYDATA.NAMES\n")
  
  tmp.correction.factor <- gbl.var$equidis.pos[1] * 10^(gbl.var$cur.exp - 9)	* (1 + config.var$CONNECTING.LINES.ADJ)
  
  #DEBUG STATEMENT
  # if (config.var$VERBOSE)  cat("config.var$DISP.CORMATRIXMAP ", config.var$DISP.CORMATRIXMAP, "\n")
  if (config.var$VERBOSE)  cat("length(gbl.var$sorted.mydata.names) ", length(gbl.var$sorted.mydata.names), "\n")
  # if (config.var$VERBOSE)  cat("tmp.correction.factor ", tmp.correction.factor, "\n")
  # if (config.var$VERBOSE)  cat("gbl.var$axis.y[1] ", gbl.var$axis.y[1], "\n")
  
  if(!config.var$DISP.CORMATRIXMAP) {
    if( length(gbl.var$sorted.mydata.names) < 30) {
      list.name <- paste("   ",gbl.var$sorted.mydata.names)
      name.ref <- paste("   ",gbl.var$sorted.mydata.names[gbl.var$mydata.ref.pos])
    } else {
      list.name <- paste("",gbl.var$sorted.mydata.names)
      name.ref <- paste("",gbl.var$sorted.mydata.names[gbl.var$mydata.ref.pos])
    }
    
    rot.name <- 90
    if (config.var$CORMATRIX.FORMAT == "RAW") {
      rot.name <- -90
    }
    if(config.var$IMAGE.SIZE == 3.5) {
      if( length(list.name) < 30) {
        mydata.names.grob <- textGrob(list.name, x = gbl.var$equidis.pos, y = unit(-8, "char"),
                                      rot=rot.name, default.units = "native",
                                      gp=gpar(fontsize = (gbl.var$font.size - 1)), just=c("left"), name="mydata.names")
        
        if(config.var$DISP.COLOR.REF) {
          #Ref
          mydata.names.grob.ref <- textGrob(name.ref, x = gbl.var$equidis.pos[gbl.var$mydata.ref.pos], y = unit(-8, "char"),
                                            rot=rot.name, default.units = "native",
                                            gp=gpar(fontsize = (gbl.var$font.size - 1),col="darkorchid1"), just=c("left"), name="mydata.names.ref")
        }
      } else {
        mydata.names.grob <- textGrob(list.name, x = gbl.var$equidis.pos, y = unit(-15, "char"),
                                      rot=rot.name, default.units = "native",
                                      gp=gpar(fontsize = (gbl.var$font.size - 1)), just=c("left"), name="mydata.names")
        if(config.var$DISP.COLOR.REF) {
          #Ref
          mydata.names.grob.ref <- textGrob(name.ref, x = gbl.var$equidis.pos[gbl.var$mydata.ref.pos], y = unit(-15, "char"),
                                            rot=rot.name, default.units = "native",
                                            gp=gpar(fontsize = (gbl.var$font.size - 1),col="darkorchid1"), just=c("left"), name="mydata.names.ref")
        }
      }
      
    } else if (config.var$IMAGE.SIZE == 7) {
      if( length(gbl.var$sorted.mydata.names) < 30) {
        mydata.names.grob <- textGrob(list.name, x = gbl.var$equidis.pos, y = unit(-8, "char"),
                                      rot=rot.name, default.units = "native",
                                      gp=gpar(fontsize = (gbl.var$font.size - 2)), just=c("left"), name="mydata.names")
        
        if(config.var$DISP.COLOR.REF) {        
          #Ref
          mydata.names.grob.ref <- textGrob(name.ref, x = gbl.var$equidis.pos[gbl.var$mydata.ref.pos], y = unit(-8, "char"),
                                            rot=rot.name, default.units = "native",
                                            gp=gpar(fontsize = (gbl.var$font.size - 2),col="darkorchid1"), just=c("left"), name="mydata.names.ref")
        }
      }else {
        mydata.names.grob <- textGrob(list.name, x = gbl.var$equidis.pos, y = unit(-15, "char"),
                                      rot=rot.name, default.units = "native",
                                      gp=gpar(fontsize = (gbl.var$font.size - 2)), just=c("left"), name="mydata.names")
        
        if(config.var$DISP.COLOR.REF) {        
          #Ref
          mydata.names.grob.ref <- textGrob(name.ref, x = gbl.var$equidis.pos[gbl.var$mydata.ref.pos], y = unit(-15, "char"),
                                            rot=rot.name, default.units = "native",
                                            gp=gpar(fontsize = (gbl.var$font.size - 2),col="darkorchid1"), just=c("left"), name="mydata.names.ref")
        }
      }
    }
    
    
    if(config.var$DISP.MYDATA.NAMES) {
      grid.draw(mydata.names.grob)
      if(config.var$DISP.COLOR.REF) {
        grid.draw(mydata.names.grob.ref)
      }
    }
  }
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("FINISH DRAW.PLOT.GRID.MYDATA.NAMES\n")
}

#------------------------------------------------------
#------------- SET PARAMETER DRAW (color, symbole)-----

#------------- CREATE COLOR BAR ------------------------
create.color.bar <- function(config.var, gbl.var) {
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("START CREATE COLOR BAR\n")
  cormatrix.matrix.ref <-  gbl.var$ref
  cormatrix.matrix <- gbl.var$cormatrix.data
  if(config.var$USE.COLORS) {
    if(config.var$CORMATRIX.COLOR.SCHEME == "heat") {
      cormatrix.colors <- heat.colors(gbl.var$palette.size)
    } else if(config.var$CORMATRIX.COLOR.SCHEME == "cm") {
      cormatrix.colors <- cm.colors(gbl.var$palette.size)
    } else if(config.var$CORMATRIX.COLOR.SCHEME == "custom") {
      if(!is.null(config.var$PALETTE.FILE)) {
        tmp.colors <- read.table(config.var$PALETTE.FILE, as.is=TRUE, header=FALSE, blank.lines.skip = TRUE, fill=TRUE)
        gbl.var$palette.size <- length(tmp.colors$V1)
        custom.colors <- paste("#", tmp.colors$V1, sep="")
        cormatrix.colors <- rev(custom.colors)
      }
      else {
        stop("PALETTE.FILE must be specified\n")
      }
    }
    else if(config.var$CORMATRIX.COLOR.SCHEME == "topo") {
      cormatrix.colors <- topo.colors(gbl.var$palette.size)
    }
    else if(config.var$CORMATRIX.COLOR.SCHEME == "gray") {
      cormatrix.colors <- gray.colors(gbl.var$palette.size)
    } else if(config.var$CORMATRIX.COLOR.SCHEME == "bluetored"){
      cormatrix.colors <- colorRampPalette(c("blue","red"))(gbl.var$mydata.num)
    } else if(config.var$CORMATRIX.COLOR.SCHEME == "bluewhitered"){
      cormatrix.colors <- colorRampPalette(c("blue","white","red"))(gbl.var$mydata.num)
    }
    else {
      stop("Invalid color scheme: ", config.var$CORMATRIX.COLOR.SCHEME, "\n")
    }
  }
  else {
    cormatrix.colors <- gray.colors(gbl.var$palette.size)
  }
  
  color.bar.colors <- c(rep(NA, gbl.var$palette.size), cormatrix.colors[gbl.var$palette.size:1])
  
  #DEBUG STATEMENT
  # if (config.var$VERBOSE)  cat("color.bar.colors ", color.bar.colors, "\n")
  #height=0.125
  color.bar <- rectGrob(x=rep(seq(1:gbl.var$palette.size)/gbl.var$palette.size, 2),
                        y=rep(seq(1:2)/2, each=gbl.var$palette.size),
                        width=(1/gbl.var$palette.size),
                        height=0.075,
                        just=c("right", "top"),
                        gp=gpar(col=NULL,
                                fill=color.bar.colors,
                                cex = (gbl.var$cex.factor + 0.2),
                                lty="blank"),
                        name = "color.bar")
  
  #  color.bar.labels <- textGrob(paste(c(-1,-0.6,-0.2,0.2,0.6,1)), x=0.2*0:5, y=0.8, gp=gpar(fontsize = gbl.var$font.size), name="color.bar.labels")
  color.bar.labels <- textGrob(paste(c(1,0.6,0.2,-0.2,-0.6,-1)), x=0.2*0:5, y=0.8, gp=gpar(fontsize = gbl.var$font.size), name="color.bar.labels")
  
  cormatrix.key.vp <- viewport(x=0.2, y=-0.03, width = 0.25, height = 0.25)
  cormatrix.key <- gTree(children=gList(color.bar, color.bar.labels), name = "cormatrix.key", vp=cormatrix.key.vp)
  
  if(config.var$DISP.PHYS.DIST) {
    if(gbl.var$total.dist > 1000) {
      map.label.text.distance <- paste("Physical Distance: ", round(gbl.var$total.dist/1000, 1), " kb", sep="")
    } else {
      map.label.text.distance <- paste("Physical Distance: ", gbl.var$total.dist, " bases", sep="")
    }
  } else {
    map.label.text.distance <- " "
  }
  
  if (config.var$CORMATRIX.METHOD == "spearman") {
    map.label.text.ldtype <- "Correlation Matrix Map Type: spearman"
  } else if (config.var$CORMATRIX.METHOD == "pearson") {
    map.label.text.ldtype <- "Correlation Matrix Map Type: pearson"
  } else if (config.var$CORMATRIX.METHOD == "kendall") {
    map.label.text.ldtype <- "Correlation Matrix Map Type: kendall"
  } else {
    stop("Invalid CORMATRIX metric : ", config.var$CORMATRIX.METHOD, "\n")
  }
  
  map.label.distance <- textGrob(map.label.text.distance,
                                 x = 0.2,
                                 y = 0.1425,
                                 just=("center"),
                                 gp = gpar(fontsize = gbl.var$font.size),
                                 name="map.label.distance")
  map.label.ldtype <- textGrob(map.label.text.ldtype,
                               x = 0.2,
                               y = 0.11150,
                               just=("center"),
                               gp = gpar(fontsize = gbl.var$font.size),
                               name="map.label.ldtype")
  
  matrix.rows <- dim(cormatrix.matrix)[1]
  matrix.cols <- dim(cormatrix.matrix)[2]
  
  intervals<-2/(length(cormatrix.colors))
  
  color.intervals <- seq(-1,1,by=intervals)
  # cat(length(color.intervals),"\n")
  # cat(length(cormatrix.colors),"\n")
  tmp.color.cut <- as.character(cut(cormatrix.matrix, color.intervals, labels=as.character(cormatrix.colors)))
  color.cut.ref <- as.character(cut(cormatrix.matrix.ref, color.intervals, labels=as.character(cormatrix.colors)))
  
  
  tmp.vector <- NULL
  
  for(i in seq(0, matrix.rows - 1)) {
    for(j in seq(matrix.rows, 1, by =-1)) {
      tmp.vector <- c(tmp.vector, (matrix.rows*j-i))
    }
  }
  
  color.cut <- rep(NULL, matrix.rows^2)
  
  for(i in 1:matrix.rows^2) {
    color.cut[i] <- tmp.color.cut[tmp.vector[i]]
  }
  
  if (config.var$VERBOSE)  cat("FINISH CREATE.COLOR.BAR\n")
  
  return(list(color.cut = color.cut, color.cut.ref =color.cut.ref, cormatrix.key = cormatrix.key, map.label.ldtype = map.label.ldtype, map.label.distance = map.label.distance))
}


#-------------------CREATE LIST of COLOR----------------------------

create.color.list <- function(config.var,gbl.var) {
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("START CREATE.COLOR.LIST\n")
  
  if(config.var$USE.COLORS) {
    tmp.color.list <-  c("red", "blue", "green", "black", "orange")
  }
  else {
    tmp.color.list <-  c("grey0", "grey20", "grey40", "grey60", "grey80")
  }
  
  if(!is.null(config.var$COLOR.LIST)) {
    split.tmp.color.list <- strsplit(config.var$COLOR.LIST, ",")
    
    if(gbl.var$mydata.samples != length(split.tmp.color.list[[1]])) {
      stop("The color list must have ", gbl.var$mydata.samples, " colors separated by commas without spaces.\n")
    }
  } else {
    if(gbl.var$mydata.samples <= 5) {
      split.tmp.color.list <- list(head(tmp.color.list, gbl.var$mydata.samples))
    } else {
      stop("coMET includes a default set of colors for 5 samples. For more samples, please specify COLOR.LIST\n")
    }
  }
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("FINISH CREATE.COLOR.LIST\n")
  
  return(split.tmp.color.list)
}


#-------------------CREATE LIST of COLOR FOR LARGE DATA----------------------------

create.color.list.large <- function(config.var,gbl.var) {
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("START CREATE.COLOR.LIST.LARGE\n")
  
  if(config.var$USE.COLORS) {
    tmp.color.list <-  c("red", "blue", "green", "black", "orange")
  }
  else {
    tmp.color.list <-  c("grey0", "grey20", "grey40", "grey60", "grey80")
  }
  
  if(!is.null(config.var$COLOR.LIST.LARGE)) {
    large.split.tmp.color.list <- strsplit(config.var$COLOR.LIST.LARGE, ",")
    
    if(gbl.var$large.mydata.samples != length(large.split.tmp.color.list[[1]])) {
      stop("The color list must have ", gbl.var$large.mydata.samples, " colors separated by commas without spaces.\n")
    }
  } else {
    if(gbl.var$large.mydata.samples <= 5) {
      large.split.tmp.color.list <- list(head(tmp.color.list, gbl.var$large.mydata.samples))
    } else {
      stop("coMET includes a default set of colors for 5 samples. For more samples, please specify COLOR.LIST.LARGE\n")
    }
  }
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("FINISH CREATE.COLOR.LIST.LARGE\n")
  
  return(large.split.tmp.color.list)
}


#-------------------CREATE LIST of SYMBOLE----------------------------
create.symbol.list <- function(config.var, split.color.list, gbl.var) {
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("START CREATE.SYMBOL.LIST\n")
  
  tmp.symbol.list <- c("circle-fill", "square-fill", "diamond-fill", "triangle-fill", "circle")
  
  if(!is.na(config.var$SYMBOLS)) {
    split.tmp.symbol.list <- strsplit(config.var$SYMBOLS, ",")
    
    if(gbl.var$mydata.samples != length(split.tmp.symbol.list[[1]])) {
      stop("The symbol list must have ", gbl.var$mydata.samples, " symbol(s) separated by commas without spaces.\n")
    }
  } else {
    if(gbl.var$mydata.samples <= 5) {
      split.tmp.symbol.list <- list(head(tmp.symbol.list, gbl.var$mydata.samples))
    } else {
      stop("coMET includes a default set of colors for 5 samples. For more samples, please specify SYMBOLS\n")
    }
  }
  
  split.symbol.list <- NULL
  
  for(i in 1:gbl.var$mydata.samples) {
    if(any(grep("circle", split.tmp.symbol.list[[1]][i], ignore.case = TRUE))) {
      split.symbol.list <- c(split.symbol.list, 21)
    } else if(any(grep("square", split.tmp.symbol.list[[1]][i], ignore.case = TRUE))) {
      split.symbol.list <- c(split.symbol.list, 22)
    } else if(any(grep("diamond", split.tmp.symbol.list[[1]][i], ignore.case = TRUE))) {
      split.symbol.list <- c(split.symbol.list, 23)
    } else if(any(grep("triangle", split.tmp.symbol.list[[1]][i], ignore.case = TRUE))) {
      split.symbol.list <- c(split.symbol.list, 24)
    } else if(any(grep("NA", split.tmp.symbol.list[[1]][i], ignore.case = TRUE))) {
      split.symbol.list <- c(split.symbol.list, NA)
    } else {
      stop("Unknown symbol: ", split.tmp.symbol.list[[1]][i], "\n")
    }
  }
  
  split.fill.list <- list(rep("white", gbl.var$mydata.samples))
  
  fill.pos <- grep("fill", split.tmp.symbol.list[[1]], ignore.case = TRUE)
  
  #DEBUG STATEMENT
  # if (config.var$VERBOSE)  cat("fill.pos ", fill.pos, "\n")
  
  split.fill.list[[1]][fill.pos] <- split.color.list[[1]][fill.pos]
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("FINISH CREATE.SYMBOL.LIST\n")
  
  return(list(split.symbol.list = split.symbol.list, split.fill.list = split.fill.list))
}


#-------------------CREATE LIST of SYMBOLE FOR LARGE DATA----------------------------
create.symbol.list.large <- function(config.var, large.split.color.list, gbl.var) {
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("START CREATE.SYMBOL.LIST.LARGE\n")
  
  large.tmp.symbol.list <- c("circle-fill", "square-fill", "diamond-fill", "triangle-fill", "circle")
  
  if(!is.na(config.var$SYMBOLS.LARGE)) {
    large.split.tmp.symbol.list <- strsplit(config.var$SYMBOLS.LARGE, ",")
    
    if(gbl.var$large.mydata.samples != length(large.split.tmp.symbol.list[[1]])) {
      stop("The symbol list must have ", gbl.var$large.mydata.samples, " symbol(s) separated by commas without spaces.\n")
    }
  } else {
    if(gbl.var$large.mydata.samples <= 5) {
      large.split.tmp.symbol.list <- list(head(large.tmp.symbol.list, gbl.var$large.mydata.samples))
    } else {
      stop("coMET includes a default set of colors for 5 samples. For more samples, please specify SYMBOLS\n")
    }
  }
  
  large.split.symbol.list <- NULL
  
  for(i in 1:gbl.var$large.mydata.samples) {
    if(any(grep("circle", large.split.tmp.symbol.list[[1]][i], ignore.case = TRUE))) {
      large.split.symbol.list <- c(large.split.symbol.list, 21)
    } else if(any(grep("square", large.split.tmp.symbol.list[[1]][i], ignore.case = TRUE))) {
      large.split.symbol.list <- c(large.split.symbol.list, 22)
    } else if(any(grep("diamond", large.split.tmp.symbol.list[[1]][i], ignore.case = TRUE))) {
      large.split.symbol.list <- c(large.split.symbol.list, 23)
    } else if(any(grep("triangle", large.split.tmp.symbol.list[[1]][i], ignore.case = TRUE))) {
      large.split.symbol.list <- c(large.split.symbol.list, 24)
    } else if(any(grep("NA", large.split.tmp.symbol.list[[1]][i], ignore.case = TRUE))) {
      large.split.symbol.list <- c(large.split.symbol.list, NA)
    } else {
      stop("Unknown symbol: ", large.split.tmp.symbol.list[[1]][i], "\n")
    }
  }
  
  large.split.fill.list <- list(rep("white", gbl.var$large.mydata.samples))
  
  large.fill.pos <- grep("fill", large.split.tmp.symbol.list[[1]], ignore.case = TRUE)
  
  #DEBUG STATEMENT
  # if (config.var$VERBOSE)  cat("large.fill.pos ", large.fill.pos, "\n")
  
  large.split.fill.list[[1]][large.fill.pos] <- large.split.color.list[[1]][large.fill.pos]
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("FINISH CREATE.SYMBOL.LIST.LARGE\n")
  
  return(list(large.split.symbol.list = large.split.symbol.list, large.split.fill.list = large.split.fill.list))
}


#-------------------DRAW GENES NAMES for WEB PAGE----------------------------
draw.name.genes.web <- function(config.var,gbl.var) {
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("START DRAW.NAME.GENES.WEB\n")
  
  #DEBUG STATEMENT
  # if (config.var$VERBOSE)  cat("gbl.var$font.size ", gbl.var$font.size, "\n")
  
  popViewport()
  top.vp <- gbl.var$top.vp
  gene.name.vp.cormatrixmap <- viewport(height = 1,
                                        width = 1,
                                        layout.pos.row = 4,
                                        layout.pos.col = 3,
                                        just=c("left", "top"),
                                        name = "gene.name.vp.cormatrixmap")
  gene.name.vp.nocormatrixmap <- viewport(height = 1,
                                          width = 1,
                                          layout.pos.row = 5,
                                          layout.pos.col = 3,
                                          just=c("left", "top"),
                                          name = "gene.namen.vp.nocormatrixmap")
  
  if(config.var$DISP.CORMATRIXMAP) {
    pushViewport(vpTree(top.vp, vpList(gene.name.vp.cormatrixmap)))
    gene.name.vp <- gene.name.vp.cormatrixmap
  }
  else {
    pushViewport(vpTree(top.vp, vpList(gene.name.vp.nocormatrixmap)))
    gene.name.vp <- gene.name.vp.nocormatrixmap
  }
  
  #-------------------NAME of GENES BEGINS-----------------------
  
  name.genes <- genesNameENSEMBL(config.var$GENOME,gbl.var$mydata.chr,gbl.var$min.x,gbl.var$max.x,config.var$DATASET.GENE)
  
  legend.name.genes.list <- NULL
  if( !is.null(name.genes)) {
    
    if(nrow(name.genes) <= 10) {
      
      #All name.genes
      # if (config.var$VERBOSE)  cat("All name.genes \n")
      for(i in 1:nrow(name.genes)) {
        
        label.text <- name.genes[i,2]
        #  if (config.var$VERBOSE)  cat("Gene",i,":", label.text,"\n")
        #Ensure that labels do not overlap by moving them
        
        y.pos <- 1.8 - 0.3 * ( i - 1) 
        x.pos <- 0.75
        
        label.name.genes.text <- textGrob(label.text,
                                          x = x.pos,
                                          y = y.pos,
                                          just=("left"),
                                          gp = gpar(fontsize = (gbl.var$font.size * 0.7 ),fontface = "bold"))
        
        
        legend.name.genes.list[[(length(legend.name.genes.list) + 1)]] <- label.name.genes.text
      }
      
      class(legend.name.genes.list) <- c("gList")
      
      legend.glist.name.genes <- legend.name.genes.list
      #                               just=c("left", "top"),
      legend.name.genes <- gTree(children=legend.glist.name.genes,
                                 vp=gene.name.vp,
                                 name="legend.gene.name")
      
      grid.draw(legend.name.genes)
    }
    
    popViewport()
  }
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("FINISH DRAW.NAME.GENES.WEB\n")
  
  return(gbl.var)
}

#-------------------DRAW NAMES OF TRACKS for WEB PAGE----------------------------
draw.name.tracks.web <- function(config.var,gbl.var) {
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("START CREATE.NAME.TRACKS.WEB\n")
  
  #DEBUG STATEMENT
  # if (config.var$VERBOSE)  cat("gbl.var$font.size ", gbl.var$font.size, "\n")
  
  popViewport()
  top.vp <- gbl.var$top.vp
  name.tracks.vp.cormatrixmap <- viewport(height = 1,
                                          width = 1,
                                          layout.pos.row = 4,
                                          layout.pos.col = 1,
                                          name = "name.tracks.vp.cormatrixmap")
  name.tracks.vp.nocormatrixmap <- viewport(height = 1,
                                            width = 1,
                                            layout.pos.row = 5,
                                            layout.pos.col = 1,
                                            name = "name.tracks.vp.nocormatrixmap")
  name.tracks.vp.cormatrixmap.nopval <- viewport(height = 1,
                                                 width = 1,
                                                 layout.pos.row = 3,
                                                 layout.pos.col = 1,
                                                 name = "name.tracks.vp.nocormatrixmap.nopval")
  
  if(config.var$DISP.CORMATRIXMAP) {
    if(config.var$DISP.PVALUEPLOT) {
      pushViewport(vpTree(top.vp, vpList(name.tracks.vp.cormatrixmap)))
      name.tracks.vp <- name.tracks.vp.cormatrixmap
    } else {
      pushViewport(vpTree(top.vp, vpList(name.tracks.vp.cormatrixmap.nopval)))
      name.tracks.vp <- name.tracks.vp.cormatrixmap.nopval
    }
  }
  else {
    pushViewport(vpTree(top.vp, vpList(name.tracks.vp.nocormatrixmap)))
    name.tracks.vp <- name.tracks.vp.nocormatrixmap
  }
  
  #-------------------NAME of TRACKS BEGINS-----------------------
  
  legend.tracks.list <- NULL
  
  y.label.pos <- c(1.25,0, -0.35,-0.85,-1.15,-1.45,-1.75)
  num.tracks <- 1
  if(config.var$IMAGE.SIZE == 3.5) {
    #All tracks cf order define via create.tracks.web
    # if (config.var$VERBOSE)  cat("All tracks \n")
    y.label.pos <- c(1.25,0.15, -0.15,-0.75,-1.25,-1.65,-1.95)
  } else if(config.var$IMAGE.SIZE == 7) {
    y.label.pos <- c(1.25,0,30, -0.1, -0.7, -1.3, -1.7, -2.1)
  }
  #--- GENES ENSEMBL 
  if(has.key("geneENSEMBL", gbl.var$split.list.tracks)) {
    label.tracks.text.ensembl <- textGrob("ENSEMBL Genes",
                                          x = -0.25,
                                          y = y.label.pos[num.tracks],
                                          just=c("right"),
                                          gp = gpar(fontsize = (gbl.var$font.size*0.75),
                                                    fontface = "bold"))
    legend.tracks.list[[num.tracks]] <- label.tracks.text.ensembl
    num.tracks <- num.tracks + 1
  }
  if(has.key("transcriptENSEMBL", gbl.var$split.list.tracks)) {
    label.tracks.text.ensembl <- textGrob("ENSEMBL Transcripts",
                                          x = -0.25,
                                          y = y.label.pos[num.tracks],
                                          just=c("right"),
                                          gp = gpar(fontsize = (gbl.var$font.size*0.75),
                                                    fontface = "bold"))
    legend.tracks.list[[num.tracks]] <- label.tracks.text.ensembl
    num.tracks <- num.tracks + 1
  }
  
  #--- genes UCSC 
  if(has.key("genesUCSC",gbl.var$split.list.tracks) ){
    label.tracks.text.snp <- textGrob("genes UCSC",
                                      x = -0.25,
                                      y = y.label.pos[num.tracks],
                                      just=c("right"),
                                      gp = gpar(fontsize = (gbl.var$font.size *0.75 ),fontface = "bold"))
    legend.tracks.list[[num.tracks]] <- label.tracks.text.snp
    num.tracks <- num.tracks + 1
  }
  
  #--- xenogenes UCSC 
  if(has.key("xenogenesUCSC",gbl.var$split.list.tracks) ){
    label.tracks.text.snp <- textGrob("xeno genes UCSC",
                                      x = -0.25,
                                      y = y.label.pos[num.tracks],
                                      just=c("right"),
                                      gp = gpar(fontsize = (gbl.var$font.size *0.75 ),fontface = "bold"))
    legend.tracks.list[[num.tracks]] <- label.tracks.text.snp
    num.tracks <- num.tracks + 1
  }
  
  #--- CG Island
  if(has.key("CGI",gbl.var$split.list.tracks)) {
    label.tracks.text.cgisland <- textGrob("CG Island",
                                           x = -0.25,
                                           y = y.label.pos[num.tracks],
                                           just=c("right"),
                                           gp = gpar(fontsize = (gbl.var$font.size *0.75  ),
                                                     fontface = "bold"))
    legend.tracks.list[[num.tracks]] <- label.tracks.text.cgisland
    num.tracks <- num.tracks + 1
  }
  
  
  #--- ChromatinHMM
  if(has.key("ChromHMM",gbl.var$split.list.tracks)) {
    label.tracks.text.chromatinhmm <- textGrob("Broad ChromHMM",
                                               x = -0.25,
                                               y = y.label.pos[num.tracks],
                                               just=c("right"),
                                               gp = gpar(fontsize = (gbl.var$font.size *0.75 ),
                                                         fontface = "bold"))
    legend.tracks.list[[num.tracks]] <- label.tracks.text.chromatinhmm
    num.tracks <- num.tracks + 1
  }
  
  
  #--- DNAse
  if(has.key("DNAse",gbl.var$split.list.tracks)) {
    label.tracks.text.dnase <- textGrob("DNase Clusters",
                                        x = -0.25,
                                        y = y.label.pos[num.tracks],
                                        just=c("right"),
                                        gp = gpar(fontsize = (gbl.var$font.size *0.75 ),fontface = "bold"))
    
    legend.tracks.list[[num.tracks]] <- label.tracks.text.dnase
    num.tracks <- num.tracks + 1
  }
  
  #--- Regulation ENSEMBL
  if(has.key("RegENSEMBL",gbl.var$split.list.tracks)) {
    label.tracks.text.reg <- textGrob("Regulation ENSEMBL",
                                      x = -0.25,
                                      y = y.label.pos[num.tracks],
                                      just=c("right"),
                                      gp = gpar(fontsize = (gbl.var$font.size *0.75 ),fontface = "bold"))
    legend.tracks.list[[num.tracks]] <- label.tracks.text.reg
    num.tracks <- num.tracks + 1
  }
  
  #--- SNP
  if(has.key("SNP",gbl.var$split.list.tracks) ){
    label.tracks.text.snp <- textGrob("SNP UCSC",
                                      x = -0.25,
                                      y = y.label.pos[num.tracks],
                                      just=c("right"),
                                      gp = gpar(fontsize = (gbl.var$font.size *0.75 ),fontface = "bold"))
    legend.tracks.list[[num.tracks]] <- label.tracks.text.snp
    num.tracks <- num.tracks + 1
  }
  
  
  #--- SNPstoma
  if(has.key("SNPstoma",gbl.var$split.list.tracks)) {
    label.tracks.text.chromatinhmm <- textGrob("SNP stomatic cells",
                                               x = -0.25,
                                               y = y.label.pos[num.tracks],
                                               just=c("right"),
                                               gp = gpar(fontsize = (gbl.var$font.size *0.75 ),fontface = "bold"))
    
    legend.tracks.list[[num.tracks]] <- label.tracks.text.chromatinhmm
    num.tracks <- num.tracks + 1
  }
  
  
  #--- SNPstru
  if(has.key("SNPstru",gbl.var$split.list.tracks)) {
    label.tracks.text.dnase <- textGrob("structural SNP",
                                        x = -0.25,
                                        y = y.label.pos[num.tracks],
                                        just=c("right"),
                                        gp = gpar(fontsize = (gbl.var$font.size *0.75 ),fontface = "bold"))
    
    legend.tracks.list[[num.tracks]] <- label.tracks.text.dnase
    num.tracks <- num.tracks + 1
  }
  #--- SNPstrustoma
  if(has.key("SNPstrustoma",gbl.var$split.list.tracks)) {
    label.tracks.text.dnase <- textGrob(" stomatic stuctural SNP",
                                        x = -0.25,
                                        y = y.label.pos[num.tracks],
                                        just=c("right"),
                                        gp = gpar(fontsize = (gbl.var$font.size *0.75 ),fontface = "bold"))
    
    legend.tracks.list[[num.tracks]] <- label.tracks.text.dnase
    num.tracks <- num.tracks + 1
  }
  
  
  
  #--- ISCA
  if(has.key("ISCA",gbl.var$split.list.tracks)) {
    label.tracks.text.reg <- textGrob("ISCA",
                                      x = -0.25,
                                      y = y.label.pos[num.tracks],
                                      just=c("right"),
                                      gp = gpar(fontsize = (gbl.var$font.size *0.75 ),fontface = "bold"))
    
    legend.tracks.list[[num.tracks]] <- label.tracks.text.reg
    num.tracks <- num.tracks + 1
  }
  
  #--- COSMIC
  if(has.key("COSMIC",gbl.var$split.list.tracks)) {
    label.tracks.text.reg <- textGrob("COSMIC",
                                      x = -0.25,
                                      y = y.label.pos[num.tracks],
                                      just=c("right"),
                                      gp = gpar(fontsize = (gbl.var$font.size *0.75 ),fontface = "bold"))
    
    legend.tracks.list[[num.tracks]] <- label.tracks.text.reg
    num.tracks <- num.tracks + 1
  }
  
  
  #--- GAD
  if(has.key("GAD",gbl.var$split.list.tracks)) {
    label.tracks.text.reg <- textGrob("GAD",
                                      x = -0.25,
                                      y = y.label.pos[num.tracks],
                                      just=c("right"),
                                      gp = gpar(fontsize = (gbl.var$font.size *0.75 ),fontface = "bold"))
    
    legend.tracks.list[[num.tracks]] <- label.tracks.text.reg
    num.tracks <- num.tracks + 1
  }
  
  
  #--- ClinVar
  if(has.key("ClinVar",gbl.var$split.list.tracks) ){
    label.tracks.text.snp <- textGrob("ClinVar",
                                      x = -0.25,
                                      y = y.label.pos[num.tracks],
                                      just=c("right"),
                                      gp = gpar(fontsize = (gbl.var$font.size *0.75 ),fontface = "bold"))
    legend.tracks.list[[num.tracks]] <- label.tracks.text.snp
    num.tracks <- num.tracks + 1
  }
  
  #--- GeneReviews
  if(has.key("GeneReviews",gbl.var$split.list.tracks) ){
    label.tracks.text.snp <- textGrob("GeneReviews",
                                      x = -0.25,
                                      y = y.label.pos[num.tracks],
                                      just=c("right"),
                                      gp = gpar(fontsize = (gbl.var$font.size *0.75 ),fontface = "bold"))
    legend.tracks.list[[num.tracks]] <- label.tracks.text.snp
    num.tracks <- num.tracks + 1
  }
  
  #--- GWAS
  if(has.key("GWAS",gbl.var$split.list.tracks) ){
    label.tracks.text.snp <- textGrob("GWAS",
                                      x = -0.25,
                                      y = y.label.pos[num.tracks],
                                      just=c("right"),
                                      gp = gpar(fontsize = (gbl.var$font.size *0.75 ),fontface = "bold"))
    legend.tracks.list[[num.tracks]] <- label.tracks.text.snp
    num.tracks <- num.tracks + 1
  }
  
  #--- ClinVarCNV
  if(has.key("ClinVarCNV",gbl.var$split.list.tracks) ){
    label.tracks.text.snp <- textGrob("ClinVarCNV",
                                      x = -0.25,
                                      y = y.label.pos[num.tracks],
                                      just=c("right"),
                                      gp = gpar(fontsize = (gbl.var$font.size *0.75 ),fontface = "bold"))
    legend.tracks.list[[num.tracks]] <- label.tracks.text.snp
    num.tracks <- num.tracks + 1
  }
  
  
  #--- GC content
  if(has.key("GCcontent",gbl.var$split.list.tracks) ){
    label.tracks.text.snp <- textGrob("GC content",
                                      x = -0.25,
                                      y = y.label.pos[num.tracks],
                                      just=c("right"),
                                      gp = gpar(fontsize = (gbl.var$font.size *0.75 ),fontface = "bold"))
    legend.tracks.list[[num.tracks]] <- label.tracks.text.snp
    num.tracks <- num.tracks + 1
  }
  
  #--- Repeat Element
  if(has.key("RepeatElt",gbl.var$split.list.tracks) ){
    label.tracks.text.snp <- textGrob("Repeat elements",
                                      x = -0.25,
                                      y = y.label.pos[num.tracks],
                                      just=c("right"),
                                      gp = gpar(fontsize = (gbl.var$font.size *0.75 ),fontface = "bold"))
    legend.tracks.list[[num.tracks]] <- label.tracks.text.snp
    num.tracks <- num.tracks + 1
  }
  

  
  class(legend.tracks.list) <- c("gList")
  
  legend.glist.tracks <- legend.tracks.list
  #                         just=c("left", "top"),
  legend.tracks <- gTree(children=legend.glist.tracks,
                         vp=name.tracks.vp,
                         name="legend.tracks")
  
  grid.draw(legend.tracks)
  
  popViewport()
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("FINISH DRAW.NAME.TRACKS.WEB\n")
  
  return(gbl.var)
}

#------------------ CREATE TRACKS for USER ------------------------------
create.tracks.user <- function(config.var,gbl.var){
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("START CREATE.TRACKS.USER\n")
  
  if(!is.null(config.var$BIOFEAT.USER.FILE)){
    split.biofeature.data.user.file <- gbl.var$split.biofeature.data.user.file
    split.biofeature.data.user.type <- gbl.var$split.biofeature.data.user.type
    split.biofeature.data.user.type.plot <- gbl.var$split.biofeature.data.user.type.plot
    
    num.data=0
    listtracks_user <- list()
    track.biofeat <- NULL
    cur.biofeat <- 0
    for(i in 1:length(split.biofeature.data.user.file[[1]])) {
      cur.biofeat <- i
      cur.biofeat.data <- split.biofeature.data.user.file[[1]][i]
      cur.biofeat.type <- split.biofeature.data.user.type[[1]][i]
      if(cur.biofeat.type == "GeneRegion"){
        track.biofeat <-  GeneRegionTrack(range=cur.biofeat.data, genome=config.var$GENOME, chromosome=gbl.var$mydata.chr)
      }else if (cur.biofeat.type == "Annotation") {
        track.biofeat <-  AnnotationTrack(range=cur.biofeat.data, genome=config.var$GENOME, chromosome=gbl.var$mydata.chr)
      }else if (cur.biofeat.type == "DATA") {
        num.data <- num.data + 1
        cur.biofeat.type.plot <- split.biofeature.data.user.type.plot[[1]][num.data]
        track.biofeat <-  DataTrack(range=cur.biofeat.data, genome=config.var$GENOME, type=cur.biofeat.type.plot, chromosome=gbl.var$mydata.chr)
      }else {
        stop("TYPE of TRACK UNKNOWN \n")
      }
      if (cur.biofeat == 0){
        listtracks_user=list(track.biofeat)
      }else{
        listtracks_user = c(listtracks_user, track.biofeat)
      }
    }
    
    gbl.var$listtracks_user <- c(listtracks_user)
  }
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("FINISH CREATE.TRACKS.USER\n")
  
  return(gbl.var)
  
}


#-------------------CREATE TRACKS for WEB PAGE----------------------------
create.tracks.web <- function(config.var,gbl.var) {
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("START CREATE.TRACKS.WEB\n")
  
  
  listtracks_gviz <-list()
  # if (config.var$VERBOSE)  cat("test",is.hash(gbl.var$split.list.tracks))
  #---- Genome Axis 
  if(has.key("genomeAxis", gbl.var$split.list.tracks)) {
    gtrack <- GenomeAxisTrack()
    
    if(length(listtracks_gviz) == 0) {
      listtracks_gviz <- list(gtrack)
    } else {
      listtracks_gviz <- c(listtracks_gviz,gtrack)
    }
  }
  
  #--- GENES ENSEMBL
  if(has.key("geneENSEMBL", gbl.var$split.list.tracks)) {
    ENSEMBLtrack <- genesENSEMBL(config.var$GENOME,gbl.var$mydata.chr,
                                 gbl.var$min.x,gbl.var$max.x,showId=TRUE)
    
    if(length(listtracks_gviz) == 0) {
      listtracks_gviz <- list(ENSEMBLtrack)
    } else {
      listtracks_gviz <- c(listtracks_gviz,ENSEMBLtrack)
    }
  }
  
  #--- transcript ENSEMBL
  if(has.key("transcriptENSEMBL",gbl.var$split.list.tracks)) {
    tENSEMBLtrack <- transcriptENSEMBL(config.var$GENOME,gbl.var$mydata.chr
                                       ,gbl.var$min.x,gbl.var$max.x,showId=TRUE)
    if(length(listtracks_gviz) == 0) {
      listtracks_gviz <- list(tENSEMBLtrack)
    } else {
      listtracks_gviz <- c(listtracks_gviz,tENSEMBLtrack)
    }
  }
  
  #--- gene UCSC
  if(has.key("genesUCSC",gbl.var$split.list.tracks)) {
    dnasetrack<-knownGenesUCSC(config.var$GENOME,gbl.var$mydata.chr,gbl.var$min.x,
                               gbl.var$max.x,showId=TRUE)
    if(length(listtracks_gviz) == 0) {
      listtracks_gviz <- list(dnasetrack)
    } else {
      listtracks_gviz <- c(listtracks_gviz,dnasetrack)
    }
  }
  
  #------ xeno ref
  if(has.key("xenogenesUCSC",gbl.var$split.list.tracks) ){
    gctrack <-xenorefGenesUCSC(config.var$GENOME,gbl.var$mydata.chr,gbl.var$min.x,
                               gbl.var$max.x,showId=TRUE)
    if(length(listtracks_gviz) == 0) {
      listtracks_gviz <- list(gctrack)
    } else {
      listtracks_gviz <- c(listtracks_gviz,gctrack)
    }
  }
  
  #--- CG Island
  if(has.key("CGI",gbl.var$split.list.tracks)) {
    cgitrack<-cpgIslandsUCSC(config.var$GENOME,gbl.var$mydata.chr,gbl.var$min.x,gbl.var$max.x)
    if(length(listtracks_gviz) == 0) {
      listtracks_gviz <- list(cgitrack)
    } else {
      listtracks_gviz <- c(listtracks_gviz,cgitrack)
    }
  }
  
  #--- ChromatinHMM
  if(has.key("ChromHMM",gbl.var$split.list.tracks)) {
    chromatintrack <- chromatinHMMAll(config.var$GENOME,gbl.var$mydata.chr,
                                      gbl.var$min.x,gbl.var$max.x,
                                      gbl.var$mySession,track.name="Broad ChromHMM",
                                      pattern=config.var$PATTERN.REGULATION)
    if(length(listtracks_gviz) == 0) {
      listtracks_gviz <- list(chromatintrack)
    } else {
      listtracks_gviz <- c(listtracks_gviz,chromatintrack)
    }
  }
  
  #--- Broad Histone
  if(has.key("BroadHistone",gbl.var$split.list.tracks)) {
    chromatintrack <- HistoneAll(config.var$GENOME,gbl.var$mydata.chr,
                                 gbl.var$min.x,gbl.var$max.x,
                                 gbl.var$mySession,track.name="Broad histone",
                                 pattern=config.var$PATTERN.REGULATION)
    if(length(listtracks_gviz) == 0) {
      listtracks_gviz <- list(chromatintrack)
    } else {
      listtracks_gviz <- c(listtracks_gviz,chromatintrack)
    }
  }
  
  #--- DNAse
  if(has.key("DNAse",gbl.var$split.list.tracks)) {
    dnasetrack<-DNAseUCSC(config.var$GENOME,gbl.var$mydata.chr,gbl.var$min.x,
                          gbl.var$max.x,gbl.var$mySession)
    if(length(listtracks_gviz) == 0) {
      listtracks_gviz <- list(dnasetrack)
    } else {
      listtracks_gviz <- c(listtracks_gviz,dnasetrack)
    }
  }
  
  
  #---- Regulation
  if(has.key("RegENSEMBL",gbl.var$split.list.tracks)) {
    regtrack<-regulationBiomart(config.var$GENOME,gbl.var$mydata.chr,gbl.var$min.x,gbl.var$max.x)
    if(length(listtracks_gviz) == 0) {
      listtracks_gviz <- list(regtrack)
    } else {
      listtracks_gviz <- c(listtracks_gviz,regtrack)
    }
  }
  
  #---- structural variation
  if(has.key("SNPstru",gbl.var$split.list.tracks)) {
    structrack<-structureBiomart(gbl.var$mydata.chr,gbl.var$min.x,gbl.var$max.x,"*",
                                 config.var$DATASET.STRU)
    if(length(listtracks_gviz) == 0) {
      listtracks_gviz <- list(structrack)
    } else {
      listtracks_gviz <- c(listtracks_gviz,structrack)
    }
  }
  
  #---- stomatic structural variation
  if(has.key("SNPstrustoma",gbl.var$split.list.tracks)) {
    structrack<-structureBiomart(gbl.var$mydata.chr,gbl.var$min.x,gbl.var$max.x,"*",
                                 config.var$DATASET.STRU.STOMA)
    if(length(listtracks_gviz) == 0) {
      listtracks_gviz <- list(structrack)
    } else {
      listtracks_gviz <- c(listtracks_gviz,structrack)
    }
  }
  
  #---- SNP stomatic cell
  if(has.key("SNPstoma",gbl.var$split.list.tracks)) {
    snpstomaENSEMBLtrack<-snpBiomart(gbl.var$mydata.chr,gbl.var$min.x,gbl.var$max.x,
                                     config.var$DATASET.SNP.STOMA,
                                     title="Stomatic Short Variation")
    if(length(listtracks_gviz) == 0) {
      listtracks_gviz <- list(snpstomaENSEMBLtrack)
    } else {
      listtracks_gviz <- c(listtracks_gviz,snpstomaENSEMBLtrack)
    }
  }
  
  #---- SNP
  if(has.key("SNP",gbl.var$split.list.tracks) ){
    snptrack <- snpLocationsUCSC(config.var$GENOME,gbl.var$mydata.chr,gbl.var$min.x,
                                 gbl.var$max.x,config.var$VERSION.DBSNP)
    if(length(listtracks_gviz) == 0) {
      listtracks_gviz <- list(snptrack)
    } else {
      listtracks_gviz <- c(listtracks_gviz,snptrack)
    }
  }
  
  #------ ISCA
  if(has.key("ISCA",gbl.var$split.list.tracks) ){
    iscatrack <-ISCATrack(config.var$GENOME,gbl.var$mydata.chr,gbl.var$min.x,gbl.var$max.x,
                          gbl.var$mySession, table.name="iscaPathogenic",showId=FALSE)
    if(length(listtracks_gviz) == 0) {
      listtracks_gviz <- list(iscatrack)
    } else {
      listtracks_gviz <- c(listtracks_gviz,iscatrack)
    }
  }
  
  
  #------ COSMIC
  if(has.key("COSMIC",gbl.var$split.list.tracks) ){
    cosmictrack <-COSMICTrack(config.var$GENOME,gbl.var$mydata.chr,
                              gbl.var$min.x,gbl.var$max.x,showId=FALSE)
    if(length(listtracks_gviz) == 0) {
      listtracks_gviz <- list(cosmictrack)
    } else {
      listtracks_gviz <- c(listtracks_gviz,cosmictrack)
    }
  }
  
  #------ GAD
  if(has.key("GAD",gbl.var$split.list.tracks) ){
    iscatrack <-GADTrack(config.var$GENOME,gbl.var$mydata.chr,gbl.var$min.x,gbl.var$max.x,
                         showId=FALSE)
    if(length(listtracks_gviz) == 0) {
      listtracks_gviz <- list(iscatrack)
    } else {
      listtracks_gviz <- c(listtracks_gviz,iscatrack)
    }
  }
  
  #------ clinic Variant
  if(has.key("ClinVar",gbl.var$split.list.tracks) ){
    clinVariant<-ClinVarMainTrack(config.var$GENOME,gbl.var$mydata.chr,gbl.var$min.x,
                                  gbl.var$max.x)
    if(length(listtracks_gviz) == 0) {
      listtracks_gviz <- list(clinVariant)
    } else {
      listtracks_gviz <- c(listtracks_gviz,clinVariant)
    }
  }
  
  if(has.key("ClinVarCNV",gbl.var$split.list.tracks) ){
    clinCNV<-ClinVarCnvTrack(config.var$GENOME,gbl.var$mydata.chr,gbl.var$min.x,gbl.var$max.x)
    if(length(listtracks_gviz) == 0) {
      listtracks_gviz <- list(clinCNV)
    } else {
      listtracks_gviz <- c(listtracks_gviz,clinCNV)
    }
  }
  
  #------ GWAS Variant
  if(has.key("GWAS",gbl.var$split.list.tracks) ){
    gwastrack <-GWASTrack(config.var$GENOME,gbl.var$mydata.chr,gbl.var$min.x,gbl.var$max.x)
    if(length(listtracks_gviz) == 0) {
      listtracks_gviz <- list(gwastrack)
    } else {
      listtracks_gviz <- c(listtracks_gviz,gwastrack)
    }
  }
  
  #------ GeneReviews
  if(has.key("GeneReviews",gbl.var$split.list.tracks) ){
    geneRtrack <-GeneReviewsTrack(config.var$GENOME,gbl.var$mydata.chr,gbl.var$min.x,
                                  gbl.var$max.x)
    if(length(listtracks_gviz) == 0) {
      listtracks_gviz <- list(geneRtrack)
    } else {
      listtracks_gviz <- c(listtracks_gviz,geneRtrack)
    }
  }
  
  #------ GC content
  if(has.key("GCcontent",gbl.var$split.list.tracks) ){
    gctrack <-gcContent(config.var$GENOME,gbl.var$mydata.chr,gbl.var$min.x,
                        gbl.var$max.x)
    if(length(listtracks_gviz) == 0) {
      listtracks_gviz <- list(gctrack)
    } else {
      listtracks_gviz <- c(listtracks_gviz,gctrack)
    }
  }
  
  
  
  #------ BioUSER
  if(!is.null(config.var$BIOFEAT.USER.FILE)){
    if(length(gbl.var$listtracks_user) == 0) {
      listtracks_gviz <- gbl.var$listtracks_user
    } else {
      listtracks_gviz <- c(listtracks_gviz,gbl.var$listtracks_user)
    }
  }
  
  #------ Repeat element 
  if(has.key("RepeatElt",gbl.var$split.list.tracks) ){
    reptrack <-RepeatMaskerTrack(config.var$GENOME,gbl.var$mydata.chr,gbl.var$min.x,
                        gbl.var$max.x,showId=TRUE)
    if(length(listtracks_gviz) == 0) {
      listtracks_gviz <- list(reptrack)
    } else {
      listtracks_gviz <- c(listtracks_gviz,reptrack)
    }
  }
  
  
  #   if (config.var$ZOOM) {
  #     listtracks_gviz=c(list(tENSEMBLtrack,cgitrack),chromatintrack,dnasetrack,regtrack,snptrack)
  #   }else {
  #     listtracks_gviz=c(list(ENSEMBLtrack,cgitrack),chromatintrack,dnasetrack,regtrack,snptrack)
  #   }
  
  
  gbl.var$listtracks_gviz <- listtracks_gviz
  
  #DEBUG STATEMENT
  if (config.var$VERBOSE)  cat("FINISH CREATE.TRACKS.WEB\n")
  
  return(gbl.var)
  
}
