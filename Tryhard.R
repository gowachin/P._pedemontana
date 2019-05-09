manual()

library(pbapply) # uberclass but maybe not usefull

# abba baba ####

table = data.frame(AML = c("A","A","B","A"),
                   PT1 = c("A","B","A","B"),
                   VL2 = c("A","B","B","A"),
                   HS3 = c("B","A","A","B"))
table

ab_ba = function(table) {
  r = "no.inf"
  if ( table[1] == table[4] & table[2] == table[3] & table[1] != table[2] ) {r = "abba"
  } else if ( table[1] == table[3] & table[2] == table[4] & table[1] != table[2] ) {r = "baba"
  } #else {r = NA}
  return(r)
}

x = as.factor(unlist(apply(table, 1, function(table) ab_ba(table))))

abba = sum(x == "abba") ; abba
baba = sum(x == "baba") ; baba
x

d <- (abba - baba) / (abba + baba) ; d


file$.csv
y = file$.pop ; y


table = substr(as.matrix(readr::read_delim(file$.csv,"\t", escape_double = FALSE, trim_ws = TRUE)[,-c(1:9)]),1,3)
colnames(table)
file$.ind

table = rbind(as.character(file$.pop),table)


c1 = sample(colnames(table[,which(table[1,]=="apennina")]),1) ;c1
c2 = sample(colnames(table[,which(table[1,]=="pedemontana")]),1) ;c2
c3 = sample(colnames(table[,which(table[1,]=="valgau")]),1) ;c3
root = sample(colnames(table[,which(table[1,]=="cottia")]),1) ; root

manus = cbind(table[,which(colnames(table) == c1 )],
              table[,which(colnames(table) == c2 )],
              table[,which(colnames(table) == c3 )],
              table[,which(colnames(table) == root )] )
colnames(manus) = c(c1,c2,c3,root)
head(manus)

hap1 = substr(as.matrix(manus),1,1)
hap2 = substr(as.matrix(manus),3,3)
homo = hap1 == hap2
homo
manus = manus[]
manus[which(homo==F)] = "."
manus = as.data.frame(manus)
manus$inf = as.factor(unlist(apply(manus, 1, function(manus) ab_ba(manus))))
summary(manus$inf)




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



