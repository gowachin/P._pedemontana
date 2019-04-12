#install.packages("vcfR")
library(vcfR)
library(readr)

#vcf.cov10.fullmat.SNP = read.vcfR("data_vcf/freebayes_-F0_3-n10-m13_-q20_mincov10_90samples_SNPs_only.vcf", checkFile = F)
#vcf.cov20.fullmat.SNP = read.vcfR("data_vcf/freebayes_-F0_3-n10-m30_-q20_mincov20_90samples_SNPs_only.vcf", checkFile = F)

mincov10_90samples_CSV <- read_delim("data_vcf/freebayes_-F0_3-n10-m13_-q20_mincov10_90samples_SNPs_onlyCSV.csv",
                                     "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(mincov10_90samples_CSV)
mincov20_90samples_CSV <- read_delim("data_vcf/freebayes_-F0_3-n10-m30_-q20_mincov20_90samples_SNPs_onlyCSV.csv",
                                     "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(mincov20_90samples_CSV)

#colnames(mincov20_90samples_CSV) == colnames(mincov10_90samples_CSV) # mêmes colnames pour deux tableaux
#vcf = colnames(mincov20_90samples_CSV)[1:9]

vcf = c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")
Eryth = c('AMB','AML','AOL'
          ,'CS1','CP1','CP4'
          ,'DGB','DRL'
          ,'DMB','HC1','HGL','HS2','HP1','HPB'
          ,'PT1','PV1','GA2','GA4'
          ,'VR3','VR1','VL2','VB1'
          )
pop = c("apennina", "apennina"
        ,"apennina"
        ,"cottia","cottia","cottia"
        ,"hirsuta"
        ,"daonensis","daonensis"
        ,"hirsuta","hirsuta","hirsuta","hirsuta","hirsuta"
        ,"pedemontana","pedemontana"
        ,"villosa","villosa","villosa","villosa"
        ,"valgau","valgau"
)
Eryth = c(vcf,Eryth)
# without the outgroup: AP1 (P. lutea)
Erythro_mincov10 = mincov10_90samples_CSV [,which(colnames(mincov10_90samples_CSV) %in% Eryth)]

Erythro_mincov20 = mincov20_90samples_CSV [,which(colnames(mincov20_90samples_CSV) %in% Eryth)]

rm(mincov10_90samples_CSV,mincov20_90samples_CSV )

#permet de virer lignes vides ####
clean = function(data) {
  n.col= dim(data)[2]
  compte = apply(data[,10:n.col],1,function(data) length( levels(as.factor(substr(as.character(data),1,3)))[!levels(as.factor(substr(as.character(data),1,3))) %in% "."] ) >1)
  #exterminate = apply(data[,10:n.col],1,paste,collapse ="")
  #compte = exterminate == paste(rep(".",23),collapse = "")
  print(summary(compte))
  resum = c(pourc.clean = (1-sum(compte)/length(compte) )*100)
  print(resum)

  end = data[which(compte == T),]
  #end = data[which(compte == F),]
  return(end)
}

inform_mincov10 = clean(Erythro_mincov10) ; dim(inform_mincov10)
#   Mode   FALSE    TRUE
#logical  608327   16481
inform_mincov20 = clean(Erythro_mincov20) ; dim(inform_mincov20)
#Mode   FALSE    TRUE
#logical   84836     410

# travail ####
#inform_mincov20[1,-c(1:9)][c(1,9)]
#data = Erythro_mincov20
#n.col= dim(data)[2]
#summary(apply(data[,10:n.col],1,function(data) length(levels(as.factor(substr(as.character(data),1,3)))[!levels(as.factor(substr(as.character(data),1,3))) %in% "."]) >1) )
#levels(as.factor(substr(as.character(inform_mincov20[1,-c(1:9)]),1,3)))[!levels(as.factor(substr(as.character(inform_mincov20[1,-c(1:9)]),1,3))) %in% "."]

rm(Erythro_mincov10,Erythro_mincov20)

# analyse des positions le long du genome ####

x = as.numeric(paste(substr(as.character(inform_mincov20$CHROM),7,20) ,  as.character(inform_mincov20$POS) ,sep="."))
plot(x = log(x), y= log(inform_mincov20$QUAL) ,cex =0.1)

#permet de virer les lignes avec seulement un certain nombre de variants ####
tri = function(data,n.r = 0,n.c = 0,quiet = F, r= F,p = F) { # data est un data frame type inform_mincov10 , n = pourcentage de présence du snp , quiet est le rendu
  #tri sur les lignes
  n.col = dim(data)[2] ; n.row = dim(data)[1]
  matrix = data[,10:n.col]
  matrix = matrix == "."

  row  = apply(matrix,1,sum)
  row = 1-row/(n.col-9) >= n.r

  data = data[which(row == T),]

  #tri sur les individus
  n.col = dim(data)[2] ; n.row = dim(data)[1]
  matrix = data[,10:n.col]
  matrix = matrix == "."

  col  = apply(matrix,2,sum)
  col = 1-col/(n.row) >= n.c

  matrix = data[,10:n.col]
  matrix = matrix[,which(col == T)]
  data = cbind(data[,1:9],matrix)

  resum = c(r.pourc.end = sum(row)/length(row)*100 ,
            r.n.site = sum(row),
            c.pourc.end = sum(col)/length(col)*100 ,
            c.n.ind = sum(col))

  if(quiet == F) print(resum)
  if (r== T & p == F) return(data)
  if (r== T & p == T) return(resum)
}

test = tri(data = inform_mincov10,n.r=0.95,n.c = 0.8, r = T) # avec ces seuils on vire l'individu AML car trop d'info manquantes pour cet individu
colnames(test[,-c(1:9)])
plot(log(test$QUAL), cex = 0.1, ylim = c(log(0.001),20)) # virer des points avec pas assez bon Phred??? qu'est-ce qu'il représente??

write.table(test,"data_vcf/tryhard.csv",sep = "\t", quote = F, row.names=F)
system(" ./CSV_to_VCF.sh data_vcf/tryhard.csv ; mv data_vcf/tryhard.csv data_vcf/tryhard.vcf" )

# DAPC ####
par(mfrow = c(1,1))
library(adegenet)
library(vcfR)
library(poppr, parallel)
tryhard = read.vcfR("data_vcf/tryhard.vcf", checkFile = T)
tryhard2 = vcfR2genlight(tryhard, n.cores = 7)
tryhard2@ind.names
pop(tryhard2) <- as.factor(pop)

myFreq <- glMean(tryhard2)
hist(myFreq, proba=TRUE, col="gold", xlab="Allele frequencies",
     main="Distribution of (second) allele frequencies")
temp <- density(myFreq)
lines(temp$x, temp$y*1.8,lwd=3)

myFreq <- c(myFreq, 1-myFreq)
hist(myFreq, proba=TRUE, col="darkseagreen3", xlab="Allele frequencies",
     main="Distribution of allele frequencies", nclass=20)
temp <- density(myFreq, bw=.05)
lines(temp$x, temp$y*2,lwd=3)

pca1 <- glPca(tryhard2)
scatter(pca1, posi="bottomleft")
library(ape)
tre <- nj(dist(as.matrix(tryhard2)))
tre
plot(tre, typ="fan", cex=0.7)

myCol <- colorplot(pca1$scores,pca1$scores, transp=TRUE, cex=4)
abline(h=0,v=0, col="grey")
add.scatter.eig(pca1$eig[1:40],2,1,2, posi="bottomleft", inset=.05, ratio=.3)
plot(tre, typ="fan", show.tip=FALSE)
tiplabels(pch=20, col=myCol, cex=4)
title("NJ tree of the US influenza data")

dapc1 <- dapc(tryhard2, n.pca=10, n.da=2)
scatter(dapc1,scree.da=FALSE, bg="white", posi.pca="topright", legend=TRUE,
        txt.leg=paste("group", 1:7))#, col=c("red","blue"))

myCol <- c("darkblue","purple","green","orange","red","blue","darkgreen")

compoplot(dapc1, col=myCol,lab="", txt.leg=paste("group", 1:7), ncol=2)


library( hierfstat)
tryhard3 = vcfR2genind(tryhard)
pop(tryhard3) <- as.factor(pop)
fstat(tryhard3, fstonly=TRUE)

pairwise.fst(tryhard3)

#While a large number of loci are nearly fixed (frequencies close to 0 or 1), there is an
#appreciable number of alleles with intermediate frequencies and therefore susceptible to
#contain interesting biological signal.

#Warning message: In vcfR2genlight(tryhard, n.cores = 7) :
#Found 895 loci with more than two alleles.Objects of class genlight only support loci with two alleles.
#895 loci will be omitted from the genlight object.

#Missing data should ideally be i) not too numerous and ii) randomly distributed in the dataset. In a situation like yours,
#individuals are more precious than markers, so I would discard loci with a majority of NAs, and briefly check the structure
#of the remaining missing entries.

#NAs are basically replaced to the mean allele frequency. This means individuals with NAs will tend to
#be placed closer to the origin. Also, individuals with similar patterns of NAs will be seen as more
#similar than they probably are in reality.

#If you really have a big missing value problem, and lot of NAs you cannot discard, one possibility
#would be to get a matrix of 1 and 0 where '1' indicate NAs, and do the PCA of this.
#If you obtain a structure, then this is a sign of problem - your NAs are not randomly distributed.


popNames(tryhard2)
ploidy(tryhard2) = 2
glPlot(tryhard2, posi="topleft")

#calcul de distance
tryhard2.pop.dist <- poppr::bitwise.dist(tryhard2)
tryhard2.pop.dist

grp = find.clusters(tryhard2, max.n.clust = round(nInd(tryhard2))-1)#, n.pca = length(pop(tryhard2)))
names(grp) ; grp$Kstat ; grp$stat ; grp$grp ; grp$size
table(pop(tryhard2), grp$grp) ; table.value(table(pop(tryhard2), grp$grp), col.lab=paste("inf", 1:6),row.lab=paste("ori", 1:6))
dapc1 <- dapc(tryhard2, grp$grp)
dapc1 ; scatter(dapc1)

dapc1 = dapc(tryhard2, var.contrib = TRUE, scale = FALSE, n.pca = length(pop(tryhard2)), n.da = nPop(tryhard2) - 1)

dapc1 <- dapc(tryhard2)
dapc1 ;

assignplot(dapc1, only.grp=NULL, subset=NULL, new.pred=NULL, cex.lab=.75,pch=3)

scatter(dapc1, posi.da="topleft", bg="white", pch=17:22)

myCol <- c("darkblue","purple","green","orange","red","blue","darkgreen","yellow","cyan")
scatter(dapc1, posi.da="topright", bg="white",
        pch=17:22, cstar=0, col=myCol, scree.pca=TRUE,
        posi.pca="topleft")

# PCAdapt analysis ####
#install.packages("pcadapt")
library(pcadapt)
filename <- read.pcadapt("data_vcf/tryhard.vcf", type = "vcf")
x <- pcadapt(input = filename, K = 15)
plot(x, option = "screeplot")
plot(x, option = "scores", pop = pop)
plot(x, option = "scores", i = 2, j = 3, pop = pop)
plot(x, option = "scores", i = 3, j = 4, pop = pop)

#Another option to choose the number of PCs is based on the ‘score plot’ that displays
#population structure. The score plot displays the projections of the individuals onto
#the specified principal components. Using the score plot, the choice of K can be limited
#to the values of K that correspond to a relevant level of population structure.

#When population labels are known, individuals of the same populations can be displayed
#with the same color using the pop argument, which should contain the list of indices of
#the populations of origin. In the geno3pops example, the first population is composed of
#the first 50 individuals, the second population of the next 50 individuals, and so on.
#Thus, a vector of indices or characters (population names) that can be provided to the argument pop should look like this:

# LEA analysis ####
# PGDSpider
system("cd ; cd Téléchargements/PGDSpider_2.1.1.5/ ; ./PGDSpider2.sh")

spider = function(input,inFORM,output,outFORM) {
  #realise cette commande ci :
  #system("cd ; java -Xmx1024m -Xms512M -jar Téléchargements/PGDSpider_2.1.1.5/PGDSpider2-cli.jar -inputfile Bureau/BEE/Stage/Pedemontana/data_vcf/tryhard.vcf -inputformat VCF -outputfile Bureau/BEE/Stage/Pedemontana/data_vcf/tryhard.str -outputformat STRUCTURE -spid Bureau/BEE/Stage/Pedemontana/data_vcf/Spid_VCF_STRUCTURE.spid")

  if (inFORM == "VCF" & outFORM == "STRUCTURE") {spid = "Spid_VCF_STRUCTURE.spid"}
  if (inFORM == "VCF" & outFORM == "PED") {spid = "Spid_VCF_PED.spid"}


  command = paste("cd ; java -Xmx1024m -Xms512M -jar Téléchargements/PGDSpider_2.1.1.5/PGDSpider2-cli.jar -inputfile Bureau/BEE/Stage/Pedemontana/",
                  input, " -inputformat ", inFORM ," -outputfile Bureau/BEE/Stage/Pedemontana/", output, " -outputformat ", outFORM ,
                  " -spid Bureau/BEE/Stage/Pedemontana/data_vcf/",spid,sep = "")
  print(command)
  system(command)
}

spider("data_vcf/tryhard.vcf","VCF","data_vcf/tryhard.str","STRUCTURE")

system("sed -i '1d' data_vcf/tryhard.str ")

library(LEA)
library(ade4)
library(vegan)

source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")

struct2geno(file = "data_vcf/tryhard.str", TESS = FALSE, diploid = T, FORMAT = 2,extra.row = 0, extra.col = 2, output = "data_vcf/tryhard.geno")
#ne marche pas car premiere ligne en trop
#permet de créer le fichier format geno 23 individuals and 25086 markers. (SNP)

obj  <- snmf("data_vcf/tryhard.geno", K = 1:14, entropy = T, ploidy = 2,
             CPU = 7,repetitions = 10, project= "new", alpha=100)
# Choix du K optimal (20 runs)
ID =Eryth

par(mfrow = c(1,1))
plot(obj, col = "blue", pch=1,cex=0.5)
beep(3)
color = c("orange","violet","lightgreen","red","blue","green","cyan","grey","black","yellow","darkgreen")

obj.snmf = snmf("data_vcf/tryhard.geno", K = 11, alpha = 100, project = "new")
qmatrix = Q(obj.snmf, K = 11)
barplot(t(qmatrix), col = color, border = NA, space = 0, xlab = "Individuals", ylab = "Admixture coefficients",
        names.arg =ID, las = 2)

Pop = function(K) {
obj.snmf = snmf("data_vcf/tryhard.geno", K = K, alpha = 100, project = "new",iterations = 200,
                CPU = 7)
qmatrix = Q(obj.snmf, K = K)
barplot(t(qmatrix), col = color, border = NA, space = 0,xlab = "Individuals", ylab = "Admixture coefficients",
        names.arg =ID, las = 2)}

par(mfrow = c(3,4))
for (i in 1:12) Pop(i) ;beep(3)
# adegenet ####
library(parallel)


# admixture ####

#install.packages("bedr")
library(bedr)
tryhard.vcf = read.vcf("data_vcf/tryhard.vcf")
tryhard.bed = vcf2bed(tryhard.vcf)

# je n'arrive pas à faire de fichier .bed

 admixture = function(input,K) {
  #realise cette commande ci :
  #system(" cd ; Téléchargements/admixture_linux-1.3.0/admixture Bureau/BEE/Stage/Pedemontana/data_vcf/tryhard.ped 7 ")
  command = paste("cd ; Téléchargements/admixture_linux-1.3.0/admixture Bureau/BEE/Stage/Pedemontana/",
                  input, " ", K,sep = "")
  print(command)
  system(command)
}
 spider("data_vcf/tryhard.vcf","VCF","data_vcf/tryhard.ped","PED")
admixture("data_vcf/tryhard.ped",7)

# plot de l'évolution du pourcentage de données en fonction des seuils ####
par(mfrow = c(2,2))
#seuil de lignes
perc = seq(0,1,by = 0.05)
res = matrix(perc,nrow = length(perc), ncol = 5) ; colnames(res) = c("stringence","info.garde","n.sites","info.garde","n.ind")
for( i in 1:length(perc)) {
  res[i,2:5] = tri(data = inform_mincov10,n.r=perc[i], r = T,p=T)
}
plot(x = res[,1], y = res[,2], xlab = "seuil de stringence", ylab = "info gardée", main = "seuil par lignes")
abline(h = 5, col = "red") ; text(x = res[,1], y = res[,2], res[,3], pos = 3)

perc = seq(0,1,by = 0.05)
res = matrix(perc,nrow = length(perc), ncol = 5) ; colnames(res) = c("stringence","info.garde","n.sites","info.garde","n.ind")
for( i in 1:length(perc)) {
  res[i,2:5] = tri(data = inform_mincov20,n.r=perc[i], r = T,p=T)
}
plot(x = res[,1], y = res[,2], xlab = "seuil de stringence", ylab = "info gardée", main = "seuil par lignes")
abline(h = 5, col = "red") ; text(x = res[,1], y = res[,2], res[,3], pos = 3)

#seuil d'ind
perc = seq(0,1,by = 0.05)
res = matrix(perc,nrow = length(perc), ncol = 5) ; colnames(res) = c("stringence","info.garde","n.sites","info.garde","n.ind")
for( i in 1:length(perc)) {
  res[i,2:5] = tri(data = inform_mincov10,n.c=perc[i], r = T,p=T)
}
plot(x = res[,1], y = res[,4], xlab = "seuil de stringence", ylab = "info gardée", main = "seuil par ind")
abline(h = 5, col = "red") ; text(x = res[,1], y = res[,4], res[,5], pos = 3)

perc = seq(0,1,by = 0.05)
res = matrix(perc,nrow = length(perc), ncol = 5) ; colnames(res) = c("stringence","info.garde","n.sites","info.garde","n.ind")
for( i in 1:length(perc)) {
  res[i,2:5] = tri(data = inform_mincov20,n.c=perc[i], r = T,p=T)
}
plot(x = res[,1], y = res[,4], xlab = "seuil de stringence", ylab = "info gardée", main = "seuil par ind")
abline(h = 5, col = "red") ; text(x = res[,1], y = res[,4], res[,5], pos = 3)
par(mfrow = c(1,1))





