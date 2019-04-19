manual()

# assign("a",get("x"))

# phrapl ####

# installer gmp https://gmplib.org/ avec dll, tar puis in repertory
#	./configure
#make
#make check		<= VERY IMPORTANT!!
#	make install

#sudo apt-get install libfreetype6-dev
devtools::install_github('bomeara/phrapl')


if(!require("ape")){install.packages("ape",repos="http://cran.rstudio.com")}
if(!require("partitions")){install.packages("partitions",repos="http://cran.rstudio.com")}
if(!require("lattice")){install.packages("lattice",repos="http://cran.rstudio.com")}
if(!require("polynom")){install.packages("polynom",repos="http://cran.rstudio.com")}
if(!require("gmp")){install.packages("gmp",repos="http://cran.rstudio.com")}
if(!require("rgenoud")){install.packages("rgenoud",repos="http://cran.rstudio.com")}
if(!require("parallel")){install.packages("parallel",repos="http://cran.rstudio.com")}
if(!require("optimx")){install.packages("optimx",repos="http://cran.rstudio.com")}
if(!require("igraph")){install.packages("igraph",repos="http://cran.rstudio.com")}
if(!require("numDeriv")){install.packages("numDeriv",repos="http://cran.rstudio.com")}
if(!require("nloptr")){install.packages("nloptr",repos="http://cran.rstudio.com")}
if(!require("Matrix")){install.packages("Matrix",repos="http://cran.rstudio.com")}
if(!require("rgl")){install.packages("rgl",repos="http://cran.rstudio.com")}
if(!require("RColorBrewer")){install.packages("RColorBrewer",repos="http://cran.rstudio.com")}
if(!require("igraph")){install.packages("igraph",repos="http://cran.rstudio.com")}
if(!require("diagram")){install.packages("diagram",repos="http://cran.rstudio.com")}
if(!require("binom")){install.packages("binom",repos="http://cran.rstudio.com")}
if(!require("P2C2M")){install.packages("P2C2M",repos="http://cran.rstudio.com")}

library(phrapl)
library(rgl)

data(TestData)
migrationArray[[4]]

#rajout des growth map car le jeu de base ne les contient pas!!!
for (i in 1:6) {
migrationArray[[i]]$growthMap = matrix(c(0,0,0,0,NA,0),ncol =2)
}

PlotModel( migrationArray[[4]])

# full matrix ####

library(ape)
library("adegenet")
library("hierfstat")

#library(seqinr)
#library(phangorn)
#fullmat=read.phyDat("data/MAP13_COV10_fullmatrix.phy",format="phylip", type="DNA")
fullmat=read.dna(file='data/MAP13_COV10_fullmatrix.fasta',format='fasta')
fullmat

row.names(fullmat) = substr(as.character(n),1,3)

dim(fullmat) #88 ind x 1200171 SNP
int = c('AP1',
        'AMB','AML','AOL',
        'CS1','CP1','CP4',
        'DGB','DRL',
        'DMB','HC1','HGL','HS2','HP1','HPB',
        'PT1','PV1','GA2','GA4',
        'VR3','VR1','VL2','VB1')
Erythro.fullmat=fullmat[which(row.names(fullmat)%in%int),] # with one outgroup: AP1 (P. lutea)
Erythro.fullmat
dim(Erythro.fullmat)
write.dna(Erythro.fullmat,file='data/Erythrodrosum_23samples_MAP13_COV10_fullmatrix.fasta',format='fasta')


ind <- as.character(row.names(fullmat)) # use later with adegenet (individual labels)
#population <- as.character(Mydata$state) # use later with adegenet (population labels)
#county <- Mydata$county



