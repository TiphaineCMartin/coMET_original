\name{comet}
\alias{comet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualize EWAS results in a genomic region of interest
}
\description{
coMET is an R-based package to visualize EWAS (epigenome-wide association scans) results in a genomic region of interest. The main feature of coMET is to plot the the significance level of EWAS results in the selected region, along with correlation in DNA methylation values between CpG sites in the region. The coMET package generates plots of phenotype-association, co-methylation patterns, and a series of annotation tracks.
}
\usage{
comet(mydata.file = NULL, mydata.format = "site", mydata.type = "file",
    mydata.large.file = NULL, mydata.large.format = "site",
    mydata.large.type = "listfile", cormatrix.file = NULL,
    cormatrix.method = "spearman", cormatrix.format = "raw",
    cormatrix.color.scheme = "bluewhitered",cormatrix.conf.level=0.05,
    cormatrix.sig.level= 1, cormatrix.adjust="none",
    cormatrix.type = "listfile", mydata.ref = NULL,
    start = NULL, end = NULL, zoom = FALSE, lab.Y = "log", 
    pval.threshold = 1e-05,pval.threshold.2 = 0,disp.pval.threshold = 1, 
    disp.association = FALSE, disp.association.large = FALSE,
    disp.region = FALSE, disp.region.large = FALSE,
    disp.beta.association = FALSE, disp.beta.association.large = FALSE, factor.beta = 0.3,
    symbols = "circle-fill", symbols.large = NA, sample.labels = NULL, sample.labels.large = NULL,
    use.colors = TRUE , disp.color.ref = TRUE, color.list = NULL, color.list.large = NULL,
    disp.mydata = TRUE, biofeat.user.file = NULL, biofeat.user.type = NULL,
    biofeat.user.type.plot = NULL, genome = "hg19", dataset.gene = "hsapiens_gene_ensembl",
    tracks.gviz = NULL, tracks.ggbio = NULL, tracks.trackviewer = NULL,
    disp.mydata.names = TRUE, disp.color.bar = TRUE, disp.phys.dist = TRUE,
    disp.legend = TRUE, disp.marker.lines = TRUE, disp.cormatrixmap = TRUE,
    disp.pvalueplot =TRUE, disp.type = "symbol", disp.mult.lab.X = FALSE,
    disp.connecting.lines = TRUE, palette.file = NULL, image.title = NULL,
    image.name = "coMET", image.type = NULL, image.size = 3.5, fontsize.gviz=5, font.factor = 1,
    symbol.factor = NULL, print.image = TRUE, connecting.lines.factor = 1.5,
    connecting.lines.adj = 0.01, connecting.lines.vert.adj = -1,
    connecting.lines.flex = 0, config.file = NULL, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{mydata.file}{
 Name of the info file describing the coMET parameters
}
    \item{mydata.format}{
    Format of the input data in mydata.file. There are 4 different options: site, region, site_asso, region_asso.
}
    \item{mydata.type}{
    Format of mydata.file. There are 2 different options: FILE or MATRIX.
}
    \item{mydata.large.file}{
    Name of additional info files describing the coMET parameters. File names should be comma-separated. It is optional, but if you add some, they need to be file(s) in tabular format with a header. Additional info file can be a list of CpG sites with/without Beta value (DNA methylation level) or direction sign. If it is a site file then it is mandatory to have the 4 columns as shown below with headers in the same order. Beta can be the 5th column(optional) and it can be either a numeric value (positive or negative values) or only direction sign ("+", "-"). The number of columns and their types are defined but the option mydata.large.format.
}
    \item{mydata.large.format}{
    Format of additional data to be visualised in the p-value plot. Format should be comma-separated. There are 4 different options for each file: site, region, site_asso, region_asso.
}
    \item{mydata.large.type}{
    Format of mydata.large.file. There are 2 different options: listfile or listdataframe.
}
    \item{cormatrix.file}{
    Name of the raw data file or the pre-computed correlation matrix file. It is mandatory and has to be a file in tabular format with an header.
}
    \item{cormatrix.method}{
    Options for calculating the correlation matrix: spearman, pearson and kendall
}
    \item{cormatrix.format}{
    Format of the input cormatrix.file. TThere are two options: raw file (raw if CpG sites are by column and samples by row or raw_rev if CpG site are by row and samples by column) and pre-computed correlation matrix (cormatrix)
}
    \item{cormatrix.color.scheme}{
    Color scheme options: heat, bluewhitered, cm, topo, gray, bluetored
}
    \item{cormatrix.conf.level}{
    Alpha level for the confidence interval. Default value= 0.05. CI will be the alpha/2 lower and upper values.
}
    \item{cormatrix.sig.level}{
    Significant level to visualise the correlation. If the correlation has a pvalue under the significant level, the correlation will be colored in "goshwhite", else the color is related to the correlation level and the color scheme choosen.Default value =1.
}
    \item{cormatrix.adjust}{
    indicates which adjustment for multiple tests should be used. "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".Default value="none"
}
    \item{cormatrix.type}{
    Format of cormatrix.file. There are 2 different options: listfile or listdataframe.
}
    \item{mydata.ref}{
 The name of the referenceomic feature (e.g. CpG-site) listed in mydata.file
}
    \item{start}{
    The first nucleotide position to be visualised. It could be bigger or smaller than the first position of our list of omic features.
}
    \item{end}{
    the last nucleotide position to be visualised. It has to be bigger than the value in the option start, but it could be smaller or bigger than the last position of our list of omic features.
}
    \item{zoom}{
 Default=False
}
    \item{lab.Y}{
 Scale of the y-axis. Options: log or ln
}
    \item{pval.threshold}{
 Significance threshold to be displayed as a red dashed line
}
  \item{pval.threshold.2}{
 the second significance threshold to be displayed as a orange dashed line
}
    \item{disp.pval.threshold}{
 Display only the findings that pass the value put in disp.pval.threshold
}
    \item{disp.association}{
    This logical option works only if mydata.file contains the effect direction (mydata.format=site_asso or region_asso). The value can be TRUE or FALSE: if FALSE (default), for each point of data in the p-value plot, the color of symbol is the color of co-methylation pattern between the point and the reference site; if TRUE, the effect direction is shown. If the association is positive, the color is the one defined with the option color.list. On the other hand, if the association is negative, the color is the opposed color.
}
    \item{disp.association.large}{
    This logical option works only if mydata.large.file contains the effect direction (mydata.large.format=site_asso or region_asso). The value can be TRUE or FALSE: if FALSE (default), for each point of data in the p-value plot, the color of symbol is the color of co-methylation pattern between the point and the reference site; if TRUE, the effect direction is shown. If the association is positive, the color is the one defined with the option color.list.large. On the other hand, if the association is negative, the color is the opposed color.
}
   \item{disp.region}{
    This logical option works only if mydata.file contains regions (mydata.format=region or region_asso). The value can be TRUE or FALSE (default). If TRUE, the genomic element will be shown by a continuous line with the color of the element, in addition to the symbol at the center of the region. If FALSE, only the symbol is shown.
}
    \item{disp.region.large}{
    This logical option works only if mydata.large.file contains regions (mydata.large.format=region or region_asso). The value can be TRUE or FALSE (default). If TRUE, the genomic element will be shown by a continuous line with the color of the element, in addition to the symbol at the center of the region. If FALSE, only the symbol is shown.
}
    \item{disp.beta.association}{
    This logical option works only if mydata.file contains the effect direction (mydata.format=site_asso or region_asso). The value can be TRUE or FALSE: if FALSE (default), for each point of data in the p-value plot, the size of symbol is the default size of symbole; if TRUE, the effect direction is shown.
}
    \item{disp.beta.association.large}{
    This logical option works only if mydata.large.file contains the effect direction (mydata.large.format=site_asso or region_asso). The value can be TRUE or FALSE: if FALSE (default), for each point of data in the p-value plot, the size of symbol is ththe default size of symbole; if TRUE, the effect direction is shown.
}
    \item{factor.beta}{
    Factor to visualise the size of beta. Default value = 0.3.
    }
    \item{symbols}{
    The symbol shown in the p-value plot. Options: circle, square, diamond, triangle. symbols can be filled by appending -fill, e.g. square-fill. Example: circle,diamond-fill,triangle
}
    \item{symbols.large}{
    The symbol to visualise the data defined in mydata.large.file. Options: circle, square, diamond, triangle; symbols can either be filled or not filled by appending -fill e.s., square-fill. Example: circle,diamond-fill,triangle
}
    \item{sample.labels}{
    Labels for the sample described in mydata.file to include in the legend
}
    \item{sample.labels.large}{
    Labels for the sample described in mydata.large.file to include in the legend
}
    \item{use.colors}{
    Use the colors defined or use the grey color scheme
}
    \item{disp.color.ref}{
    Logical option TRUE or FALSE (TRUE default). if TRUE, the connection line related to the reference probe is in purple, if FALSE if the connection line related to the reference probe stay black.
}
    \item{color.list}{
    List of colors for displaying the P-value symbols related to the data in mydata.file
}
    \item{color.list.large}{
    List of colors for displaying the P-value symbols related to the data in mydata.large.file
}
    \item{disp.mydata}{
    logical option TRUE or FALSE. TRUE (default). If TRUE, the P-value plot is shown; if FALSE the plot will be defined by GViz
}
    \item{biofeat.user.file}{
    Name of data file to visualise in the tracks. File names should be comma-separated.
}
    \item{biofeat.user.type}{
    Track type, where multiple tracks can be shown (comma-separated): DataTrack, AnnotationTrack, GeneregionTrack.
}
    \item{biofeat.user.type.plot}{
    Format of the plot if the data are shown with the Gviz's function called DataTrack (comma-separated)
}
    \item{genome}{
    The human genome reference file. e.g. "hg19" for Human genome 19 (NCBI 37), "grch37" (GRCh37),"grch38" (GRCh38)
}
    \item{dataset.gene}{
    The gene names from ENSEMBL. e.g. hsapiens_gene
}
    \item{tracks.gviz}{
    list of tracks created by Gviz.
}
    \item{tracks.ggbio}{
    list of tracks created by ggbio.
}
    \item{tracks.trackviewer}{
    list of tracks created by track viewer.
}
    \item{disp.mydata.names}{
    logical option TRUE or FALSE. If True (default), the names of the CpG sites are displayed.
}
    \item{disp.color.bar}{
    Color legend for the correlation matrix (range -1 to 1). Default: blue-white-red
}
    \item{disp.phys.dist}{
    logical option (TRUE or FALSE). TRUE (default).Display the bp distance on the plots
}
\item{disp.legend}{
    logical option TRUE or FALSE. TRUE (default) Display the sample labels and corresponding symbols on the lower right side
}
    \item{disp.marker.lines}{
    logical option TRUE or FALSE. TRUE (default), if FALSE the red line for pval.threshold is not shown
}
    \item{disp.cormatrixmap}{
    logical option TRUE or FALSE. TRUE (default), if FALSE correlation matrix is not shown
}
    \item{disp.pvalueplot}{
    logical option (TRUE or FALSE). TRUE (default), if FALSE the pvalue plot is not shown
}
    \item{disp.type}{
    Default: symbol
}
    \item{disp.mult.lab.X}{
    logical option TRUE or FALSE. FALSE (default).Display evenly spaced X-axis labels; up to 5 labels are shown.
}
    \item{disp.connecting.lines}{
    logical option TRUE or FALSE. TRUE (default) displays connecting lines between p-value plot and correlation matrix
}
    \item{palette.file}{
    File that contains color scheme for the heatmap. Colors are hexidecimal HTML color codes; one color per line; if you do not want to use this option, use the color defined by the option cormatrix.color.scheme
}
    \item{image.title}{
    Title of the plot
}
    \item{image.name}{
    The path and the name of the plot file without extension. The extension will be added by coMET depending on the option image.type.
}
    \item{image.type}{
    Options: pdf or eps
}
    \item{image.size}{
    Default: 3.5 inches. Possible sizes : 3.5 or 7
}
    \item{fontsize.gviz}{
    Font size of writing in annotation track. Default value =5
   }
    \item{font.factor}{
    Font size of the sample labels. Range: 0-1
}
    \item{symbol.factor}{
    Size of the symbols. Range: 0-1
}
    \item{print.image}{
    Print image in file or not.
}
    \item{connecting.lines.factor}{
    Length of the connecting lines. Range: 0-2
}
    \item{connecting.lines.adj}{
    Position of the connecting lines horizontally. Negative values shift the connecting lines to the left and positive values shift the lines to the right. Range: (-1;1) option -1 means no connecting lines.
}
    \item{connecting.lines.vert.adj}{
    Position of the connecting lines vertically. Can be used to vertically adjust the position of the connecting lines in relation to the CpG-site names. Negative value shift the connecting lines down. Range: (-0.5 - 0), option -1 mean the default value related to the plot size (-0.5 for 3.5 plot size; -0.7 for 7.5 plot size)
}
    \item{connecting.lines.flex}{
    Adjusts the spread of the connecting lines. Range: 0-2
}
    \item{config.file}{
    Configuration file contains the values of these options instead of defining these by command line. It is a file where each line is one option. The name of option and its value are separated by "=". If there are multiple values such as for the option list.tracks or the options for additional data, you need to separated them by a "comma" and not extra space. (i.e. list.tracks=geneENSEMBL,CGI,ChromHMM,DNAse,RegENSEMBL,SNP)
}
    \item{verbose}{
       logical option TRUE or FALSE. TRUE (default). If TRUE, shows comments.
}
}
\details{
The function is limited to visualize 120 omic features.
}
\value{
Create a plot in pdf or eps format depending to some options
}
\references{
http://epigen.kcl.ac.uk/comet/
}
\author{
Tiphaine Martin
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \code{\link{comet.web}},\code{\link{comet.list}}
}
\examples{
extdata <- system.file("extdata", package="coMET",mustWork=TRUE)
configfile <- file.path(extdata, "config_cyp1b1_zoom_4comet.txt")
myinfofile <- file.path(extdata, "cyp1b1_infofile.txt")
myexpressfile <- file.path(extdata, "cyp1b1_infofile_exprGene_region.txt")
mycorrelation <- file.path(extdata, "cyp1b1_res37_rawMatrix.txt")

chrom <- "chr2"
start <- 38290160
end <- 38303219
gen <- "hg38"

if(interactive()){
    cat("interactive")
    genetrack <-genesENSEMBL(gen,chrom,start,end,showId=TRUE)
    snptrack <- snpBiomart(chrom, start, end,
                dataset="hsapiens_snp_som",showId=FALSE)
    strutrack <- structureBiomart(chrom, start, end,
                strand, dataset="hsapiens_structvar_som")
    clinVariant<-ClinVarMainTrack(gen,chrom,start,end)
    clinCNV<-ClinVarCnvTrack(gen,chrom,start,end)
    gwastrack <-GWASTrack(gen,chrom,start,end)
    geneRtrack <-GeneReviewsTrack(gen,chrom,start,end)
    listgviz <- list(genetrack,snptrack,strutrack,clinVariant,
                 clinCNV,gwastrack,geneRtrack)
    comet(config.file=configfile, mydata.file=myinfofile, mydata.type="file",
      cormatrix.file=mycorrelation, cormatrix.type="listfile",
      mydata.large.file=myexpressfile, mydata.large.type="listfile",
      tracks.gviz=listgviz, verbose=FALSE, print.image=FALSE,disp.pvalueplot=FALSE)
} else {
    cat("Non interactive")
    data(geneENSEMBLtrack)
    data(snpBiomarttrack)
    data(ISCAtrack)
    data(strucBiomarttrack)
    data(ClinVarCnvTrack)
    data(clinVarMaintrack)
    data(GWASTrack)
    data(GeneReviewTrack)
    listgviz <- list(genetrack,snptrack,strutrack,clinVariant,
                clinCNV,gwastrack,geneRtrack)
    comet(config.file=configfile, mydata.file=myinfofile, mydata.type="file",
       cormatrix.file=mycorrelation, cormatrix.type="listfile",
        mydata.large.file=myexpressfile,  mydata.large.type="listfile",
        tracks.gviz=listgviz, verbose=FALSE, print.image=FALSE,disp.pvalueplot=FALSE)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}% __ONLY ONE__ keyword per line
