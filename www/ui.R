library(shiny)
library(coMET)

shinyUI(fluidPage(
  
  titlePanel("coMET"),
  
  sidebarLayout(
    sidebarPanel(
      h5("Info file",style = "color:red"),
      p(span("You need to load a file", style = "color:red")),
      p("Example file:",a(href="http://comet.epigen.kcl.ac.uk:3838/cyp1b1_infofile.txt", "cyp1b1_infofile.txt", target="_blank")),
   
      fileInput('datafile', 'Choose info file to upload (mandatory, max 100 omic features) :',
                accept=c('text/csv', 
                         'text/comma-separated-values,
		text/plain', 
                         '.csv')
      ),
      
      checkboxInput(inputId="defineParm",label="Define parameter for data file", FALSE),
      conditionalPanel(condition = "input.defineParm", 
                       p(span("You need to define your format", style = "color:red")),
                       selectInput("dataformat", "Format of data file:",
                                   c(" "="NULL",
                                     "site" = "site",
                                     "region" = "region",
                                     "Site with direction values" = "site_asso",
                                     "Region with directions values" = "region_asso")),
                       
                       conditionalPanel(condition = "input.dataformat == 'site_asso' | input.dataformat == 'region_asso' ",
                                        checkboxInput(inputId="dispAsso",label="Display the direction of association", TRUE),
                                        conditionalPanel(condition = "input.dispAsso",
                                                         textInput("datacolor", "Give a color:", "red")
                                        )
                       ),
                       conditionalPanel(condition = "input.dataformat == 'region' | input.dataformat == 'region_asso' ",
                                        checkboxInput(inputId="dispReg",label="Display the region ", TRUE)
                       ),
                       
                       textInput("datalabel", "Give a lable:", "Discovery"),
                       
                       selectInput("datasymb", "Define the symbol of data:" , 
                                   choices = c("circle-fill","square-fill","diamond-fill"))
      ),
      
      hr(),
      h5("Correlation matrix or Raw data",style = "color:red"),
      p(span("You need to load a file", style = "color:red")),
      p("Example file:",a(href="http://comet.epigen.kcl.ac.uk:3838/cyp1b1_res37_rawMatrix.txt", "cyp1b1_res37_rawMatrix.txt", target="_blank")),
      fileInput('corfile', 'Choose raw or correlation file to upload (mandatory):',
                accept=c('text/csv', 
                         'text/comma-separated-values,
		text/plain', 
                         '.csv')
      ),
      
      checkboxInput(inputId="definecorParm",label="Define parameter for raw or correlation file", FALSE),
      conditionalPanel(condition = "input.definecorParm",
                       p(span("You need to define your format", style = "color:red")),
                       selectInput("corformat", "Format of the correlation matrice:",
                                   c(" " = "NULL",
                                     "Correlation" = "cormatrix",
                                     "raw" = "raw",
                                     "inverse Raw" = "raw_rev")),
                       conditionalPanel(condition = "input.corformat == 'raw' | input.corformat == 'raw_rev'",
                                        selectInput("cormethod", "Method to compute the correlation:",
                                                    c(
                                                      "Spearman" = "spearman",
                                                      "Pearson" = "pearson",
                                                      "Kendall" = "kendall")),
                                        numericInput("coralphaCI", "alpha value of confidence level (for example 0.05):",0.05),
                                        numericInput("corpvalThres", " higher P-value threshold displayed (for example 0.05):", 0.05),
                                        selectInput("coradjmethod", "Method to adjust for multiple test:",
                                                    c(
                                                      "none" = "none",
                                                      "Holm" = "holm",
                                                      "Hochberg" = "hochberg",
                                                      "Hommel" = "hommel",
                                                      "Bonferroni" = "bonferroni",
                                                      "Benjamini-Hochberg" = "BH",
                                                      "Benjamini–Hochberg–Yekutieli" = "BY",
                                                      "FDR" = "fdr"
                                                      )
                                                    )
                       ),
			selectInput("corcolor", "Define the color scheme of correlation matrix:",
                                                    c(
                                                      "Blue to Red via white" = "bluewhitered",
                                                      "Heat" = "heat",
                                                      "cm" = "cm",
                                                      "topo" = "topo",
                                                      "Gray" = "gray",
                                                      "Blue to Red"= "bluetored"))

      ),
      
      hr(),
      h5("Configuration file"),
      p('If you do not want to define the parameters of coMET via the interface, you can download the example configuration file and modify according to your data. But you have to modify at least 3 parameters (', span("mydata.format (format of info file), CORMATRIX.FORMAT ( format of raw or correlation matrice), cormatrix.method (method to analyse the raw data if it is raw format)", style = "color:red"),')'),
      p("Example file:",a(href="http://comet.epigen.kcl.ac.uk:3838/config_cyp1b1_zoom_4webserver.txt", target="_blank", 
                          "config_cyp1b1_zoom_4webserver.txt")),
      fileInput('configfile', 'Choose configuration file to upload:',
                accept=c('text/plain', '.txt')
      ),
      
      hr(),
      h5("Parameters for P-value plot"),
      selectInput("genome", "Define the genome:",
			            c(
			              "Hg19" = "hg19",
			              "GRCh37" = "grch37",
			              "GRCh38" = "grch38")),
      checkboxInput(inputId="defineplotParm",label="Define parameter for plot", FALSE),
      conditionalPanel(condition = "input.defineplotParm", 
                       uiOutput("startCpG"),
                       uiOutput("stopCpG"),
                       numericInput("pvalThres", "Significant P-value threshold (for example 1e-8):",0.0000001),
                       numericInput("disppvalThres", " higher P-value threshold displayed (for example 0):", 0),
                       selectInput("scale", "Select your scale of y-values" , 
                                   choices = c("log10","ln")),
                       checkboxInput(inputId="refCpG",label="Select a CpG reference", FALSE),
                       
                       
                       # Display this only if the data file exist
                       conditionalPanel(condition = "input.refCpG == true",
                                        checkboxInput(inputId="refCpGcolor",label="Visualize the connective line of the reference CpG in purple", TRUE),
                                        uiOutput("listCpG")
                       ),
                       
                       hr(),
                       h5("Parameters of image"),
                       textInput("imagetitle", "Define the title of your plot:", "Comet plot")
      ),
      
      
      hr(),
      h5("Info file for the second set of data (e.g.= gene expression, validation data) (optional, no limitation)"),
      p("Example file:",a(href="http://comet.epigen.kcl.ac.uk:3838/cyp1b1_infofile_exprGene_region.txt", target="_blank", 
                          "cyp1b1_infofile_exprGene_region.txt")),
      fileInput('datalargefile', 'Choose other data file to upload',
                accept=c('text/csv', 
                         'text/comma-separated-values,
		text/plain', 
                         '.csv')
      ),
      
      checkboxInput(inputId="defineextraParm",label="Define parameter for extra data file", FALSE),
      conditionalPanel(condition = "input.defineextraParm", 
                       selectInput("datalargeformat", "Format of data file:",
                                   c(" " = "NULL",
                                     "site" = "site",
                                     "region" = "region",
                                     "site with direction values" = "site_asso",
                                     "region with directions values" = "region_asso")),
                       
                       conditionalPanel(condition = "input.datalargeformat == 'site_asso' | input.datalargeformat == 'region_asso' ",
                                        checkboxInput(inputId="dispAssolarge",label="Display the direction of association", TRUE),
                                        conditionalPanel(condition = "input.dispAssolarge",
                                                         textInput("datalargecolor", "Give a color:", "green")
                                        )
                       ),
                       conditionalPanel(condition = "input.datalargeformat == 'region' | input.datalargeformat == 'region_asso' ",
                                        checkboxInput(inputId="dispReglarge",label="Display the region ", TRUE)
                       ),
                       
                       textInput("datalargelabel", "Give a lable:", "Validation"),
                       selectInput("datalargesymb", "Define the symbol of data:" , 
                                   choices = c("circle-fill","square-fill","diamond-fill"))
      ),
      
      
      hr(),
      h5("Annotation track (optional)"),
      checkboxInput(inputId="defineTrack",label="Define annotation track (optional)", FALSE),
      conditionalPanel(condition = "input.defineTrack", 
                       selectInput("trackAnnot", "Select multiple annotation tracks (max 6):",
                                   c("Genes (ENSEMBL)" = "geneENSEMBL",
                                     "Transcript (ENSEMBL)" = "transcriptENSEMBL",
                                     "Genes (UCSC)" = "genesUCSC",
                                     "CpG Island (UCSC)" = "CGI",
                                     "ChromHMM Broad (UCSC)" = "ChromHMM",
                                     "DNAse (UCSC)" = "DNAse",
                                     "Regulation (ENSEMBL)" = "RegENSEMBL",
                                     "SNP (version dbSNP 142)" = "SNP",
                                     "ISCA" = "ISCA",
                                     "ClinVar Main" = "ClinVar",
                                     "ClinVar CNV" = "ClinVarCNV",
                                     "GWAS catalog" = "GWAS",
                                     "GAD variants" = "GAD",
                                     "GeneReviews" = "GeneReviews",
                                     "Genome axis" = "genomeAxis",
                                     "structural SNP" = "SNPstru",
                                     "stomatic SNP " = "SNPstoma",
                                     "stomatic structural SNP" = "SNPstrustoma",
                                     "xeno genes (UCSC)" = "xenogenesUCSC",
                                     "GC content" = "GCcontent",
                                     "COSMIC" = "COSMIC"),
                                   multiple=TRUE),
                       checkboxInput(inputId="addTrack",label="Add user-customised annotation track", FALSE),
                       conditionalPanel(condition = "input.addTrack",
                                        fileInput('annotfile', 'Choose your annotation file to upload',
                                                  accept=c('text/csv', 
                                                           'text/comma-separated-values,
		text/plain', 
                                                           '.csv',
                                                           '.bed')
                                        ),
                                        
                                        selectInput("annotformat", "Format of your annotation track:",
                                                    c(" " = "NULL",
                                                      "Gene track" = "Generegion",
                                                      "Annotation track" = "Annotation",
                                                      "Data track" = "Data")),
                                        
                                        conditionalPanel(condition = "input.annotformat == 'Data' ",
                                                         selectInput("annotplot", "Format of your Data plot:",
                                                                     c(" "= "NULL",
                                                                       "dot plot" = "p",
                                                                       "line plot" = "l",
                                                                       "dot and lines plots" = "b",
                                                                       "lines plot of average values" = "a",
                                                                       "stair steps (horizontal first)" = "s",
                                                                       "stair steps (vertical first)" = "S",
                                                                       "add grid lines" = "g",
                                                                       "add linear regression line" = "r",
                                                                       "histogram lines" = "h",
                                                                       "add loess curve" = "smooth",
                                                                       "histogram" = "histogram",
                                                                       "mountain-type plot" = "mountain",
                                                                       "polygon-type plot" = "polygon",
                                                                       "box and whisker plot" = "boxplot",
                                                                       "false color image of the summarized values" = "gradient",
                                                                       "false color image of the individual values" = "heatmap",
                                                                       "Horizon plot" = "horizon"))
                                        )
                       )
      ),
      
      hr(),
			h5("Save your image",style = "color:red"),
			textInput('plotfilename', "Filename of your plot","coMET"),
			selectInput("imagesize", "Define the size of plot:" , 
			            choices = c("7","3.5")),
      hr(),
      h5("Submit the data:",style = "color:red"),
      p('You need to click on the plot button and go the coMET plot tab.'),
      p('The creation of plot takes time relative to the time to connect UCSC and ENSEMBL'),
      actionButton("goPlot", "Plot")
    ),
    
    mainPanel(
      tabsetPanel(type = "tabs", 
                  tabPanel("coMET Home ", htmlOutput("home")),
                  tabPanel("coMET Plot", htmlOutput("cometplotUI")),
                  tabPanel("Info data", tableOutput("table")),
                  tabPanel("Correlation Matrix", tableOutput("cortable")),
                  tabPanel("Second info data", tableOutput("odatatable")),
                  tabPanel("Annotation Data", tableOutput("annottable")),
                  tabPanel("Configuration File", htmlOutput("configtableUI")),
                  tabPanel("coMET Help ", htmlOutput("help"))
      )
    )
    
  )
  
))
