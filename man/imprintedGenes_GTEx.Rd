\name{imprintedGenes_GTEx}
\alias{imprintedGenes_GTEx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a imprinted genes track from GTEx
}
\description{
Creates a track of imprinted genes from GTEx using the Gviz bioconductor package. A complete list of features and their associated colours can be found in the user guide.
}
\usage{
imprintedGenes_GTEx(gen="hg19", chr,start, end, tissues="all", classification="all",showId=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{gen}{
the name of the genome. Default value=hg19
}
  \item{chr}{
The chromosome of interest
}
  \item{start}{
The starting position in the region of interest (the smallest value)
}
  \item{end}{
The end position in the region of interest (the largest value)
}
 \item{tissues}{
  list of tissues among 33 tissues in GTEx
}
  \item{classification}{
 list of classification from 5 types (biallelic, consistent with biallelic, consistent with imprinting, imprinted, NC)
}
\item{showId}{
logical. say if we write the name of group
}
}
%% \details{
%%  ~~ If necessary, more details than the description above ~~
%% }
\value{
An AnnotationTrack object of Gviz
}
\references{
http://bioconductor.org/packages/release/bioc/html/Gviz.html

Got to BindingMotifsBiomart binding motif biomart
}
\author{
Tiphaine Martin
}
%% \note{
%%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%% }
\examples{
library("Gviz")
gen<-"hg19"
chr<- "chr6"
start <- 144251437
end <- 144330541

%all
if(interactive()){
   allIGtrack <- imprintedGenes_GTEx(gen,chr,start, end, tissues="all", classification="imprinted",showId=TRUE)
   allimprintedIGtrack <- imprintedGenes_GTEx(chr,start, end, tissues="all", classification="imprinted",showId=TRUE)
   StomachIGtrack <-imprintedGenes_GTEx(chr,start, end, tissues="Stomach", classification="all",showId=TRUE)
   PancreasIGtrack <- imprintedGenes_GTEx(chr,start, end, tissues="Pancreas", classification="all",showId=TRUE)
   PancreasimprintedIGtrack <- imprintedGenes_GTEx(chr,start, end, tissues="Pancreas", classification="biallelic",showId=TRUE)
   
   imprintinglist <- list(allIGtrack,allimprintedIGtrack,StomachIGtrack,PancreasIGtrack,PancreasimprintedIGtrack)
   
  plotTracks(,imprintinglist, from = start, to = end)
  
} else {

  data(allIGtrack)
  data(allimprintedIGtrack)
  data(StomachIGtrack)
  data(PancreasIGtrack)
  data(PancreasimprintedIGtrack)
  
   imprintinglist <- list(allIGtrack,allimprintedIGtrack,StomachIGtrack,PancreasIGtrack,PancreasimprintedIGtrack)
   
  plotTracks(imprintinglist, from = start, to = end)
}

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
