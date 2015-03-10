#Load library
library("shiny")
library("grid")
library("grDevices")
library("rtracklayer")
library("biomaRt")
library("Gviz")
library("ggbio")
library("colortools")
library("hash")
library("EBImage")
library("psych")
library("coMET")

#Normally Need to load coMET package not like that !!!!
## Need to wait that it accept in Bioconductor
#source("/home/ubuntu/git_iop/comet/Rpackage/comet/R/cometWeb.R")
#source("/home/ubuntu/git_iop/comet/Rpackage/comet/R/AnalyseFile.R")
#source("/home/ubuntu/git_iop/comet/Rpackage/comet/R/BiofeatureGraphics.R")
#source("/home/ubuntu/git_iop/comet/Rpackage/comet/R/DrawPlot.R")
#source("/home/ubuntu/git_iop/comet/Rpackage/comet/R/GeneralMethodComet.R")


# By default, the file size limit is 5MB. It can be changed by
# setting this option. Here we'll raise limit to 9MB.
options(shiny.maxRequestSize = 9*1024^2)


shinyServer(function(input, output,session) {
  filenameImage <- NULL
  
  #### CREATE LIST OF ELEMENTS
  output$listCpG <- renderUI ({
    inFile1 <- input$datafile
    listCpG <- c("")
    if (is.null(inFile1)) {
      listCpG <- c("")
      sizeListCpG <- length(listCpG) - 1
      # print(sizeListCpG,digit= 10)
      
    } else {
      dataInput <- read.csv(inFile1$datapath, header = TRUE,
                            sep = "\t",quote="")
      listDI<-dataInput[,1]
      listCpG <- sapply(listDI,as.character)
      sizeListCpG <- length(listCpG) 
      #print(sizeListCpG,digit=10);
      
    }
    selectInput("reflistCpG", paste ("among",sizeListCpG) , 
                choices = listCpG,selected = NULL,
                multiple = FALSE)
    
  })  
  
  ## Start position
  output$startCpG <- renderUI ({
    inFile1 <- input$datafile
    startCpG <- 1
    if (is.null(inFile1)) {
      p('Start position to visualise (nt):')
    } else {
      dataInput <- read.csv(inFile1$datapath, header = TRUE,
                            sep = "\t",quote="")
      if (input$dataformat == 'region' | input$dataformat == 'region_asso') {
        startCpG <- min(dataInput[,3])
      } else  {
        startCpG <- min(dataInput[,3])
      }
      numericInput("start", "Start position to visualise (nt):", startCpG)
    } 
    
  })  
  
  ## Stop position
  output$stopCpG <- renderUI ({
    inFile1 <- input$datafile
    stopCpG <- 1
    if (is.null(inFile1)) {
      p('Stop position to visualise (nt):')
    } else {
      dataInput <- read.csv(inFile1$datapath, header = TRUE,
                            sep = "\t",quote="")
      if (input$dataformat == 'region' | input$dataformat == 'region_asso') {
        stopCpG <- max(dataInput[,4])
      } else  {
        stopCpG <- max(dataInput[,3])
      }
      numericInput("stop", "Stop position to visualise (nt):", stopCpG)
    }
    
  }) 
  
  #### READ DIFFERENT FILES 
  ## DATASET
  output$table <- renderTable({
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    
    inFile <- input$datafile
    
    if (is.null(inFile))
      return(NULL)
    
    read.csv(inFile$datapath, header = TRUE,
             sep = "\t",quote="")
  })
  
  ## INFO DATA for help
  output$infoHelp <- renderTable({
    infohelp <- "/var/shiny-server/www/coMET/cyp1b1info.txt" 
    
    if (is.null(infohelp))
      return(NULL)
    
    data_infohelp <-read.csv(infohelp, header = TRUE,
                             sep = "\t", quote = "")
    data_infohelp[1:6,]
  })
  
  ## INFO DATA for help,region format
  output$inforegionHelp <- renderTable({
    infohelp <- "/var/shiny-server/www/coMET/cyp1b1info_expr_reduce_region.txt" 
    
    if (is.null(infohelp))
      return(NULL)
    
    data_infohelp <-read.csv(infohelp, header = TRUE,
                             sep = "\t", quote = "")
    data_infohelp
  })
  
  ## CONFING FILE for help,region format
  output$configFileHelp <- renderTable({
    infohelp <- "/var/shiny-server/www/coMET/config_cyp1b1_zoom.txt" 
    
    if (is.null(infohelp))
      return(NULL)
    
    data_infohelp <-read.csv(infohelp, header = TRUE,
                             sep = "\t", quote = "")
    data_infohelp
  })

 ## INFO DATA to download
 # output$downloadINFO <- downloadHandler(
#	filename =  "cyp1b1_infofile.txt",

 #  content = function(file){
#	dataINFO = "/var/shiny-server/www/coMET/cyp1b1_infofile.txt",
#	write.csv(dataINFO, file)
#	}
# )


  
  ## CORRELATION MATRIX
  output$cortable <- renderTable({
    corinFile <- input$corfile
    
    if (is.null(corinFile))
      return(NULL)
    
    read.csv(corinFile$datapath, header = TRUE,
             sep = "\t", quote = "")
  })
  
  ## CORRELATION MATRIX for help
  output$corHelp <- renderTable({
    corhelp <- "/var/shiny-server/www/coMET/cyp1b1res_37.txt" 
    
    if (is.null(corhelp))
      return(NULL)
    
    data_corhelp <-read.csv(corhelp, header = TRUE,
                            sep = "\t", quote = "")
    data_corhelp[1:6,1:6]
  })
  
  
  ## EXTRA DATASET
  output$odatatable <- renderTable({
    oinFile <- input$datalargefile
    
    if (is.null(oinFile))
      return(NULL)
    
    read.csv(oinFile$datapath, header = TRUE,
             sep = "\t", quote = "")
  })
  
  ## ANNOTATION DATA
  output$annottable <- renderTable({
    annotinFile <- input$annotfile
    
    if (is.null(annotinFile))
      return(NULL)
    
    read.csv(annotinFile$datapath, header = TRUE,
             sep = "\t", quote = "")
  })
  
  ## CONFIGURATION DATA
  output$configtable <- renderTable({
    configdatainFile <- input$configfile
    
    if (is.null(configdatainFile))
      return(NULL)
    
    configdata<-read.csv(configdatainFile$datapath, header = FALSE,
                         sep = "=", quote = "")
    colnames(configdata)<- c("Parameter","Value")
    configdata
  })
  
  
  ##### DRAW THE PLOT
  #Define the name of file
  filenameplot <- function() {
      filename = paste(input$plotfilename, "png", sep='.')
      filenameImagetmp <- tempfile(fileext=filename)
      filenameImagetmp
  }

  #Create the png
  output$cometplotImage <- renderImage({
    isize <- as.numeric(input$imagesize) * 100
    filenameImage <- filenameplot()
        png(file = filenameImage,
                  width=isize,
                  height=isize,
                  fonts=c("sans"))
      plotPrintInput()
      dev.off()
    
    # Return a list containing the filename
      listFile <- list(src = filenameImage,
           contentType = 'image/png',
           width = isize,
           height = isize,
           alt = "This is alternate text")
    listFile 
  }, deleteFile = FALSE)
  
  cometplotPrint <- renderPrint({ plotPrintInput() })
  
  plotPrintInput <- function() {
    
    datainFile <- input$datafile
    if (is.null(datainFile))
      return(NULL)
    
    cordatainFile <- input$corfile
    if (is.null(cordatainFile))
      return(NULL)
    
    largedatainFile <- input$datalargefile
    largedatainFilepath <- NULL
    if (!is.null(largedatainFile)) 
      largedatainFilepath <-  largedatainFile$datapath
    
    annotdatainFile <- input$annotfile
    annotdatainFilepath <- NULL
    if (!is.null(annotdatainFile))
      annotdatainFilepath <- annotdatainFile$datapath
    
    configdatainFile <- input$configfile
    
    listTrack <- "geneENSEMBL,CGI,ChromHMM,DNAse,RegENSEMBL,SNP"
    if (is.null(input$trackAnnot)){
      listTrack <-"geneENSEMBL,CGI,ChromHMM,DNAse,RegENSEMBL,SNP"
    } else if(!is.null(input$trackAnnot) & input$defineTrack){
      listTrack<-paste(input$trackAnnot,collapse=',')
    }
    
    dataformat <- NULL
    if(!is.null(input$dataformat) & (input$defineParm))
      dataformat <-input$dataformat
    
    dispAsso <- "FALSE"
    if(!is.null(input$dispAsso) & (input$defineParm) & input$dispAsso)
      dispAsso <- "TRUE"
    
    dispReg <- "FALSE"
    if( !is.null(input$dispReg) & (input$defineParm) & input$dispReg)
      dispReg <- "TRUE"
    
    datalab <- NULL
    if(!is.null(input$datalabel) & (input$defineParm))
      datalab <- input$datalabel
    
    datasymb <- NULL
    if(!is.null(input$datasymb) & (input$defineParm))
      datasymb <- input$datasymb
    
    datacolor <- "red"
    if(!is.null(input$datacolor) & (input$defineParm))
      datacolor <- input$datacolor
    
    
    datalargeformat <- NULL
    if(!is.null(input$datalargeformat) & (input$defineextraParm))
      datalargeformat <-input$datalargeformat
    
    displargeAsso <- "FALSE"
    if(!is.null(input$dispAssolarge) & (input$defineextraParm) & input$dispAssolarge)
      displargeAsso <- "TRUE"
    
    dispReglarge <- "FALSE"
    if( !is.null(input$dispReglarge) & (input$defineextraParm) & input$dispReglarge)
      dispReglarge <- "TRUE"
    disReglarge <- as.character(dispReglarge)
    
    datalablarge <- NULL
    if(!is.null(input$datalargelabel) & (input$defineextraParm))
      datalablarge <- input$datalargelabel
    
    datalargesymb <- NULL
    if(!is.null(input$datalargesymb) & (input$defineextraParm))
      datalargesymb <- input$datalargesymb
    
    datalargecolor <- "green"
    if(!is.null(input$datalargecolor) & (input$defineextraParm))
      datalargecolor <- input$datalargecolor
    
    cormethod <- "spearman"
    if(!is.null(input$cormethod) & (input$definecorParm))
      cormethod <- input$cormethod
    
    corformat <- "raw"
    if(!is.null(input$corformat) & (input$definecorParm))
      corformat <- input$corformat
    
    coralphaCI <- 0.05
    if(!is.null(input$coralphaCI) & (input$definecorParm))
      coralphaCI <- input$coralphaCI
    
    copvalThres <- 0.05
    if(!is.null(input$corpvalThres) & (input$definecorParm))
      corpvalThres <- input$corpvalThres
    
    coradjmethod <- "none"
    if(!is.null(input$coradjmethod) & (input$definecorParm))
      coradjmethod <- input$coradjmethod
    
    corcolor <- "bluewhitered"
    if(!is.null(input$corcolor) & (input$definecorParm))
      corcolor <- input$corcolor
    
    genomeCpG <- "hg19"
    if(!is.null(input$genome) & (input$defineplotParm))
      genomeCpG <- input$genome
    
    startCpG <- NULL
    if(!is.null(input$startCpG) & (input$defineplotParm))
      startCpG <- input$startCpG
    
    stopCpG <- NULL
    if(!is.null(input$stopCpG) & (input$defineplotParm))
      stopCpG <- input$stopCpG
    
    myrefCpG <- NULL
    if(!is.null(input$reflistCpG) & (!(input$reflistCpG == "" )) & (input$defineplotParm) & input$refCpG)
      myrefCpG <- as.character(input$reflistCpG)
    
    dispCpG <- "FALSE"
    if(!is.null(input$refCpGcolor) & (input$defineplotParm) & input$refCpGcolor) 
      dispCpG <- "TRUE"
    
    if(!is.null(input$refCpGcolor) & (!input$defineplotParm)) 
      dispCpG <- "TRUE"
    
    
    annotformat <- NULL
    if(!is.null(input$annotformat) & input$annotformat != "NULL"
       &(input$defineTrack))
      annotformat <- as.character(input$annotformat)
    
    annotplot <- NULL
    if(!is.null(input$annotplot) & input$annotplot != "NULL"
       & (input$defineTrack))
      annotplot <- as.character(input$annotplot)
    
    imagetitle <- "coMET plot"
    if(!is.null(input$imagetitle))
      imagetitle <- as.character(input$imagetitle)
    
    if (is.null(configdatainFile)) {
      #plot(1,2)
      #comet.web(mydata.file=datainFile$datapath,mydata.format=as.character(input$dataformat))
      comet.web(mydata.file=datainFile$datapath,mydata.format=dataformat,disp.association=dispAsso,
                disp.region=dispReg, sample.labels=datalab, symbols=datasymb, color.list=datacolor,
                cormatrix.file=cordatainFile$datapath,cormatrix.method=cormethod, 
                cormatrix.format=corformat, cormatrix.adjust=coradjmethod,
                cormatrix.conf.level=coralphaCI, cormatrix.sig.level=corpvalThres,
                cormatrix.color.scheme=corcolor,mydata.large.file=largedatainFilepath,
                mydata.large.format=datalargeformat,disp.association.large=displargeAsso,
                disp.region.large=dispReglarge, sample.labels.large=datalablarge, 
                symbols.large=datalargesymb, color.list.large=datalargecolor,genome=genomeCpG,
                start=startCpG,end=stopCpG,mydata.ref=myrefCpG,disp.color.ref=dispCpG,
                list.tracks=listTrack, biofeat.user.file=annotdatainFilepath, 
                biofeat.user.type=annotformat, biofeat.user.type.plot=annotplot, image.title=imagetitle, 
                print.image=FALSE, verbose=TRUE)
    } else {
      #plot(1,2)
      #mydata.large.file=largedatainFilepath,
      #configFile="/home/tmartin/git_iop/comet/Rpackage/comet/data/smoking/config_cyp1b1_zoom_local.txt"
      #comet.web(config.file=configFile,mydata.file=datainFile$datapath,cormatrix.file=cordatainFile$datapath, print.image=FALSE)
      
      comet.web(config.file=configdatainFile$datapath, mydata.file=datainFile$datapath,
                cormatrix.file=cordatainFile$datapath, mydata.large.file=largedatainFilepath, 
                biofeat.user.file=annotdatainFilepath, print.image=FALSE, verbose=TRUE)
    }
    
  }
  
  output$downloadPlot <- downloadHandler(
    filename = function() { paste(input$plotfilename, input$imageformat, sep='.') },
    content = function(filename) {
      if(input$imageformat == "pdf"){
        pdf(encoding = "ISOLatin1.enc",
            file = filenameImage,
            onefile=FALSE,
            width=as.numeric(input$imagesize),
            height=as.numeric(input$imagesize),
            paper="special")
      } 
      if(input$imageformat == "eps"){
        postscript(encoding = "ISOLatin1.enc",
                   file = filenameImage,
                   horizontal=FALSE,
                   onefile=FALSE,
                   width=as.numeric(input$imagesize),
                   height=as.numeric(input$imagesize),
                   paper="special",
                   pagecentre=TRUE,
                   fonts=c("sans"))
      }
      if(input$imageformat == "png"){
        png(file = filenameImage,
            width=as.numeric(input$imagesize),
            height=as.numeric(input$imagesize),
            fonts=c("sans"))
      }
      #lena = readImage(filenameImage)
      #display(lena)
      cometplotPrint()
      dev.off()
    }
  )
  
  
  observe({
    if(input$goPlot) {
      if(is.null(input$datafile)){
        output$cometplotUI <- renderUI({
          h5("Need to upload data file",style = "color:red")
        })
      } else if (is.null(input$corfile)) {
        output$cometplotUI <- renderUI({
          h5("Need to upload correlation file",style = "color:red")
        })
      } else if (is.null(input$configfile) &&  ((input$defineParm == FALSE) || (input$definecorParm == FALSE) || (input$defineplotParm == FALSE))) {
        output$cometplotUI <- renderUI({
          h5("Need to upload configuration file to define parameters of data file, of correlation file and of plot or define from the website",style = "color:red")
        })
      } else {
        output$cometplotUI <- renderUI({
          
          tagList(
            p('your plot is running, please wait....'),
            p('Connexion to UCSC and Ensembl'),
            hr(),    
            imageOutput("cometplotImage"),
            hr(),    
            hr(),  
            hr(),  
            hr(),  
            hr(),  
            hr(),
            hr(),
            hr(),    
            hr(),  
            hr(),  
            hr(),  
            hr(),  
            hr(),
            hr(),
            h5("Download your image",style = "color:red"),
            p('Click on the image and save the image in png or redo in pdf or eps'),
            p('It is going to take time because it need to rerun coMET and need to connect to UCSC and ENSEMBL'),
            selectInput("imageformat", "Define the format of plot:" , 
                        choices = c("pdf","eps","png")),
            downloadButton('downloadPlot', 'Download')
          )
        })
      }
    }else {          output$cometplotUI <- renderUI({
      h5("Need to upload different files, parameters, and to click on button that launch the plot",style = "color:red")
    })
    }
  })
  
  
  #### PARAMETERS of coMET PLOT
  output$configtableUI <- renderUI({
    tagList(
      h1('Parameters used in coMET plot'),
      h3('Parameters from web interface :'),
      
      #### List of parameters
      h4('Info file :'),
      p("Info File :",input$datafile),
      
      if(!is.null(input$dataformat) & (input$defineParm))
        p('Format of info file:',input$dataformat),
      
      if(!is.null(input$dispAsso) & (input$defineParm) & input$dispAsso)
        p('Show association of info file:',input$dispAsso),
      
      if(!is.null(input$dispAsso) & (input$defineParm) &  !input$dispAsso)
        p('Show association of info file: FALSE'),
      
      if( !is.null(input$dispReg) & (input$defineParm) & input$dispReg)
        p('Show region of info file:',TRUE),
      
      if( !is.null(input$dispReg) & (input$defineParm) & !input$dispReg)
        p('Show region of info file:',FALSE),
      
      if(!is.null(input$datalabel) & (input$defineParm))
        p('Label of info file:',input$datalabel),
      
      if(!is.null(input$datasymb) & (input$defineParm))
        p('Symbol of info file:', input$datasymb),
      
      if(!is.null(input$datacolor) & (input$defineParm))
        p('Color of info file:', input$datacolor),
      
      ##Correlation
      h4('Info file :'),
      p('Correlation File :', input$corfile),
      
      if(!is.null(input$cormethod) & (input$definecorParm))
        p('Method of Correlation File:',input$cormethod),
      
      if(!is.null(input$corformat) & (input$definecorParm))
        p('Format of Correlation File:',input$corformat),
      
      if(!is.null(input$corcolor) & (input$definecorParm))
        p('Color Correlation File:',input$corcolor),
      
      #Supplementary file
      h4('Supplementary file :'),
      if (!is.null(input$datalargefile)) 
        p('Supplementary File :',input$datalargefile),
      
      if(!is.null(input$datalargeformat) & (input$defineextraParm))
        p('Format of supplementary file:',input$datalargeformat),
      
      if(!is.null(input$dispAssolarge) & (input$defineextraParm) & input$dispAssolarge)
        p('Show association of supplementary file:',input$dispAssolarge),
      
      if( !is.null(input$dispReglarge) & (input$defineextraParm) & input$dispReglarge)
        p('Show region of supplementary file:',input$dispReglarge),
      
      if(!is.null(input$datalargelabel) & (input$defineextraParm))
        p('Label of supplementary file:',input$datalargelabel),
      
      if(!is.null(input$datalargesymb) & (input$defineextraParm))
        p('Symbole of supplementary file:',input$datalargesymb),
      
      if(!is.null(input$datalargecolor) & (input$defineextraParm))
        p('Color of supplementary file:',input$datalargecolor),
      
      #Configuration
      h4('Configuration file :'),
      p('Configuration File :',input$configfile),
      
      if (is.null(input$trackAnnot))
        p('list of tracks:geneENSEMBL,CGI,ChromHMM,DNAse,RegENSEMBL,SNP'),
      if(!is.null(input$trackAnnot) & input$defineTrack)
        p('list of tracks:',paste(input$trackAnnot,collapse=',')),
      
      if(!is.null(input$startCpG) & (input$defineplotParm))
        p('Start of genomic region :',input$startCpG),
      
      if(!is.null(input$stopCpG) & (input$defineplotParm))
        p('Stop of genomic region :',input$stopCpG),
      
      if(!is.null(input$reflistCpG) & (!(input$reflistCpG == "" )) & (input$defineplotParm) & input$refCpG)
        p('Name of reference element:',as.character(input$reflistCpG)),
      
      if(!is.null(input$refCpGcolor) & (input$defineplotParm) & input$refCpGcolor) 
        p('Show reference element:',input$refCpGcolor),
      if(!is.null(input$refCpGcolor) & (!input$defineplotParm))
        p('Show reference element:',input$refCpGcolor),
      
      #Annotation
      h4('user file for annotation :'),
      if (!is.null(input$annotfile))
        p('Annotation File :',input$annotfile),
      
      if(!is.null(input$annotformat) & (input$defineTrack))
        p('Format of annotation file :',as.character(input$annotformat)),
      
      if(!is.null(input$annotplot) & (input$defineTrack))
        p('Type of plot :',as.character(input$annotplot)),
      
      if(!is.null(input$imagetitle))
        p('Title of plot :',as.character(input$imagetitle)),
      ####
      
      h3('Parametres from configuration file'),
      tableOutput("configtable")
    )
  })
  

  ##### EXPLANATIONS ABOUT COMET
  output$home <- renderUI({
    tagList(
      h1('Welcome to coMET'),
      h3('Overview'),
      p('The coMET package is a web-based plotting tool and R-based package to visualize different genome-wide association scans such as EWAS (epigenome-wide association scan) results in a genomic region of interest. coMET provides a plot of the EWAS association signal and visualisation of the methylation correlation between CpG sites (co-methylation). The coMET package also provides the option to annotate the region using functional genomic information, including both user-defined features and pre-selected features based on the ',
        a(href = 'http://genome.ucsc.edu/ENCODE/', 'Encode'),  
        'project. The plot can be customized with different parameters, such as plot labels, colours, symbols, heatmap colour scheme, significance thresholds, and including reference CpG sites. Finally, the tool can also be applied to display the correlation patterns of other genomic data, e.g. gene expression array data.'
      ),
      h3('coMET webservice'),
      p('The webservice is the pre-formated web service of coMET with a reduction of parameters availlable.'),
      p('Only 100 omic features can be visualised in the correlation matrix.'),
      p('If the region is large or ENSEMBL and UCSC is busy, the creation of plot can take time.'),
      h3('Developpers'),
      p('coMET is developed by Tiphaine C. Martin in collaboration with Idil Yet, Pei-Chien Tsai, Jordana T.Bell, Department of Twin Research, Kings College London.'),
      h3('Cite'),
      p('Martin, T.C, Erte, I, Tsai, P-C, Bell, J.T.,coMET: an R plotting package to visualize regional plots of epigenome-wide association scan results,',
	a(href='http://quantgen.soc.srcf.net/qg14/', 'QG14'),
	', 2014.'),
       h3('Example of coMET plot'),
         img(src='http://comet.epigen.kcl.ac.uk:3838/minimal-cometwebPlot.jpeg'),
      h3('Contacts'),
      p('For any question, you can send an email to',
        a(href='mailto:tiphaine.martin@kcl.ac.uk;jordana.bell@kcl.ac.uk?Subject=CoMET', 'Tiphaine Martin and Jordana Bell')),
      h3('More information'),
      p('Go to the website',
        a(href='http://epigen.kcl.ac.uk/comet', ' Department of Twin Research')),
      h3('Download'),
      p('Want to download the R package :',
        a(href='https://github.com/TiphaineCMartin/coMET', ' Department of Twin Research'))
	
    )
    
  })
  
  output$help <- renderUI({
    tagList(
      h1('Welcome to coMET help'),
      h3('Format of info file (mandatory)'),
      p('It is mandatory and has to be a file in tabular format with an header.'),
      p('The info file can be a list of CpG sites with/without Beta (or direction sign). If it is a CpG-site file then it is mandatory to have 4 columns as shown below with headers in the same order. The beta or direction can be included in the 5th column (optional).'),  
      tableOutput("infoHelp"),
      
      p('Alternatively, the info file can be region-based and if so, the region-based info file must have the 5 columns (see below) with headers in this order. The beta or direction can be included in the 6th column (optional).'),
      
      tableOutput("inforegionHelp"),
      
      p('There are 4 different options for mydata.format: site (4 columns), region (5 columns), site_asso (5 columns), region_asso (6 columns).'),
      
      h3('Format of correlation matrix (mandatory)'),
      p('It is mandatory and has to be a file in tabular format with an header.'),
      p('The correlation matrix dataset can either be a pre-calculated correlation matrix or the raw data. There are two format for raw data. The format called RAW is put if the CpG sites/regions are by column and the samples are by row whereas the format called raw_rev is put if the CpG sites/regions are by row and the samples are by column. If it is a raw data then you can select the type of correlation method (spearman, kendall or pearson).'),
      p('Example of data at RAW format:'),
      tableOutput("corHelp"),
      
      h3('Format of extra info file'),
      p('It is optional file or list of files separatated by comma. Files should be in tabular format with an header'),
      p('This can be another type of info file (e.g Expression data or replication data) and it follows the same rules than the info file.'),
      
      h3('Format of annotation file'),
      p('format accepted by GViz such as BED, GTF, and GFF3 format. It is optional file.'),
      
      h3('Option of config.file'),
      p('If you would like to make your own changes to the plot you can download the configuration file from this site. After you make the relative changes you can upload it to the server again and plot.'),
      p('It is a file where each line is one option. The name of option is in capital and is separated to its value by "=". If there are multiple values such as for the option list.tracks or the options for exta data.'),
      tableOutput("configFileHelp"),
      
      
      hr(),
      h3('Option of comet.web'),
      ?comet.web
    )
  })
  
  
  
  
})
