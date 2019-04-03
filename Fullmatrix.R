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

