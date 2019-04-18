#install.packages("vcfR")
library(vcfR,readr)
library(ade4,adegenet)
library(poppr)
library(parallel)
library(ape)
library(hierfstat)
library(pcadapt)
library(LEA)
library(vegan)
library(pegas)
library(Pedemontana)

# creation sous jeu de donnees ####

# #vcf.cov10.fullmat.SNP = read.vcfR("data_vcf/freebayes_-F0_3-n10-m13_-q20_mincov10_90samples_SNPs_only.vcf", checkFile = F)
# #vcf.cov20.fullmat.SNP = read.vcfR("data_vcf/freebayes_-F0_3-n10-m30_-q20_mincov20_90samples_SNPs_only.vcf", checkFile = F)

#mincov10_90samples_CSV <- readr::read_delim("data_vcf/freebayes_-F0_3-n10-m13_-q20_mincov10_90samples_SNPs_onlyCSV.csv","\t", escape_double = FALSE, trim_ws = TRUE)
#mincov20_90samples_CSV <- readr::read_delim("data_vcf/freebayes_-F0_3-n10-m30_-q20_mincov20_90samples_SNPs_onlyCSV.csv","\t", escape_double = FALSE, trim_ws = TRUE)

# #colnames(mincov20_90samples_CSV) == colnames(mincov10_90samples_CSV) # mêmes colnames pour deux tableaux
# #vcf = colnames(mincov20_90samples_CSV)[1:9]

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
Erythv = c(vcf,Eryth)
# # without the outgroup: AP1 (P. lutea)
#Erythro_mincov10 = mincov10_90samples_CSV [,which(colnames(mincov10_90samples_CSV) %in% Erythv)]

#Erythro_mincov20 = mincov20_90samples_CSV [,which(colnames(mincov20_90samples_CSV) %in% Erythv)]

#rm(mincov10_90samples_CSV,mincov20_90samples_CSV )

# # clean is personnal fonction (package Pedemontana)
#inform_mincov10 = clean(Erythro_mincov10) ; dim(inform_mincov10)
# #Mode   FALSE    TRUE   pourc.clean
# #logical  447716  177092  71.65657
# # dim : [1] 177092     31
#inform_mincov20 = clean(Erythro_mincov20) ; dim(inform_mincov20)
# #Mode   FALSE    TRUE  pourc.clean
# #logical   60894   24352   71.43326
# # dim : [1] 24352    31

#data = Erythro_mincov20[,10:n.col]
#levels = apply(data,1, function(data) length(levels(as.factor(substr(as.character(data),1,3)))[!levels(as.factor(substr(as.character(data),1,3))) %in% "."]))
#plot(levels)
#Erythro_mincov20[which(levels > 3),]$ALT

#rm(Erythro_mincov10,Erythro_mincov20)

#write.table(inform_mincov10,"data_vcf/freebayes_-F0_3-n10-m13_-q20_mincov10_Eryth_SNPs_onlyCSV.csv",sep = "\t", quote = F, row.names=F)
#write.table(inform_mincov20,"data_vcf/freebayes_-F0_3-n10-m30_-q20_mincov20_Eryth_SNPs_onlyCSV.csv",sep = "\t", quote = F, row.names=F)

# charger subser Eryth ####

mincov10_Eryth_CSV <- readr::read_delim("data_vcf/freebayes_-F0_3-n10-m13_-q20_mincov10_Eryth_SNPs_onlyCSV.csv","\t", escape_double = FALSE, trim_ws = TRUE)
#mincov20_Eryth_CSV <- readr::read_delim("data_vcf/freebayes_-F0_3-n10-m30_-q20_mincov20_Eryth_SNPs_onlyCSV.csv","\t", escape_double = FALSE, trim_ws = TRUE)

a = c('AMB','AML','AOL') #apenina
c = c('CS1','CP1','CP4') #cottia
d = c('DGB','DRL') #daonensis
h = c('DMB','HC1','HGL','HS2','HP1','HPB') #hirsuta
p = c('PT1','PV1','GA2','GA4') #pedemontana
v = c('VR3','VR1','VL2','VB1') #villosa

Eryth10 = subset_reorder(mincov10_Eryth_CSV, c(a,p,c,v,h,d)) ; colnames(Eryth20)

Eryth10_r = rare(Eryth10, rare  = 0.05, r= T)

par(mfrow = c(1,1))
levels = apply(Eryth10[,-c(1:9)],1, function(data) levels(as.factor(substr(as.character(data),1,3)))[!levels(as.factor(substr(as.character(data),1,3))) %in% "."])
summary(as.factor(unlist(levels)))

par(mfrow = c(1,1))
levels = apply(Eryth10_r[,-c(1:9)],1, function(data) levels(as.factor(substr(as.character(data),1,3)))[!levels(as.factor(substr(as.character(data),1,3))) %in% "."])
summary(as.factor(unlist(levels)))

#tri selon QUAL in vcf
plot(log(Eryth10_r$QUAL), cex = 0.1, ylim = c(log(0.0000001),20)) # virer des points avec pas assez bon Phred??? qu'est-ce qu'il représente??
abline(h = log(20), col = "green")
Eryth10_r = Eryth10_r[which(Eryth10_r$QUAL >= 20),]
dim(Eryth10_r[which(Eryth10_r$QUAL >= 20),])

#permet de virer les lignes avec seulement un certain nombre de variants ####
#fonction tri dans le package
Eryth10_t = tri(data = Eryth10_r,n.r=0.9,n.c = 0.8, r = T) # avec ces seuils on vire l'individu AML car trop d'info manquantes pour cet individu
colnames(Eryth10_t[,-c(1:9)])

save2vcf(Eryth10_t)

# analyse des positions le long du genome ####
#x = as.numeric(paste(substr(as.character(inform_mincov20$CHROM),7,20) ,  as.character(inform_mincov20$POS) ,sep="."))
#plot(x = log(x), y= log(inform_mincov20$QUAL) ,cex =0.1)

# adegenet ####
Eryth20_v = read.vcfR("data_vcf/Eryth20_t.vcf", checkFile = T) ; Eryth20_v
pop = c("apennina", "apennina","apennina"
        ,"pedemontana","pedemontana"
        ,"valgau","valgau"
        ,"cottia","cottia","cottia"
        ,"villosa","villosa","villosa","villosa"
        ,"hirsuta","hirsuta","hirsuta","hirsuta","hirsuta","hirsuta"
        ,"daonensis","daonensis"
)

# analyse pour genind ####
Eryth20_v = read.vcfR("data_vcf/Eryth20_t.vcf", checkFile = T) ; Eryth20_v
Eryth20_g = vcfR2genind(Eryth20_v) ; Eryth20_g
row.names(Eryth20_g$tab)

#assign.pop = function(genind) {
#  for (i in 1:length(row.names(Eryth20_g$tab))){
#  }
#  pop(Eryth20_g) <- as.factor(pop)
#  return(genind)
#}

pop(Eryth20_g) <- as.factor(pop)

toto <- summary(Eryth20_g)
names(toto)
par(mfrow=c(2,1))
barplot(toto$loc.n.all, ylab="Number of alleles",
        main="Number of alleles per locus")
barplot(toto$Hexp-toto$Hobs, main="Heterozygosity: expected-observed",
        ylab="Hexp - Hobs")

#Is mean observed H significantly lower than mean expected H ? Nope
bartlett.test(list(toto$Hexp,toto$Hobs))
t.test(toto$Hexp,toto$Hobs,pair=T,var.equal=TRUE,alter="greater")

Eryth20_g_hwt <- hw.test(Eryth20_g, B=0) # KEZAKO??

fstat(Eryth20_g)
#     pop         Ind
#Total 0.1242216 -0.09033675
#pop   0.0000000 -0.24499159

#This table provides the three F statistics F st (pop/total), F it (Ind/total), and F is
#(ind/pop). These are overall measures which take into account all genotypes and all loci.

# excés d'He dans les espèces, donc on peut supposer écart a HW

matFst <- pairwise.fst(Eryth20_g)
matFst
pop(Eryth20_g)
plot(nj(matFst))
#        1          2          3          4          5          6
# 2 0.12868562
# 3 0.14535658 0.10341028
# 4 0.11000696 0.09727770 0.10829867
# 5 0.16043729 0.14541962 0.14920292 0.14778485
# 6 0.11634573 0.09884202 0.08943198 0.10876489 0.12199082
# 7 0.20172406 0.17631443 0.18029774 0.15988093 0.12156202 0.09343348
par(mfrow =c(1,1))
X <- tab(Eryth20_g, freq = TRUE, NA.method = "mean")
pca1 <- dudi.pca(X, scale = FALSE, scannf = FALSE, nf = 3)
barplot(pca1$eig[1:20], main = "PCA eigenvalues", col = heat.colors(50))

s.class(pca1$li, pop(Eryth20_g))
title("PCA of Primula Auricula")
add.scatter.eig(pca1$eig[1:20], 3,1,2)

col <- funky(15)
s.class(pca1$li, pop(Eryth20_g),xax=1,yax=2, col=transp(col,.6), axesell=FALSE,
        cstar=0, cpoint=3, grid=FALSE)

grp = find.clusters(Eryth20_g, max.n.clust = round(nInd(Eryth20_l))-1)#, n.pca = length(pop(Eryth20_l)))
names(grp) ; grp$Kstat ; grp$stat ; grp$grp ; grp$size
table(pop(Eryth20_g), grp$grp) ; table.value(table(pop(Eryth20_g), grp$grp), col.lab=paste("inf", 1:6),row.lab=paste("ori", 1:6))
dapc1 <- dapc(Eryth20_g, grp$grp)
dapc1 ; scatter(dapc1)

dapc1 = dapc(Eryth20_g, var.contrib = TRUE, scale = FALSE, n.pca = length(pop(Eryth20_g)), n.da = nPop(Eryth20_g) - 1)

dapc1 <- dapc(Eryth20_g)
dapc1 ;

assignplot(dapc1, only.grp=NULL, subset=NULL, new.pred=NULL, cex.lab=.75,pch=3)

scatter(dapc1, posi.da="topleft", bg="white", pch=17:22)

myCol <- c("darkblue","purple","green","orange","red","blue","darkgreen","yellow","cyan")
scatter(dapc1, posi.da="topright", bg="white",
        pch=17:22, cstar=0, col=myCol, scree.pca=TRUE,
        posi.pca="topleft")


# analyse pour genligth! ####
Eryth20_l = vcfR2genlight(Eryth20_v, n.cores = 7)
pop(Eryth20_l) <- as.factor(pop)

myFreq <- glMean(Eryth20_l)
hist(myFreq, proba=TRUE, col="gold", xlab="Allele frequencies",
     main="Distribution of (second) allele frequencies")
temp <- density(myFreq)
lines(temp$x, temp$y*1.8,lwd=3)

myFreq <- c(myFreq, 1-myFreq)
hist(myFreq, proba=TRUE, col="darkseagreen3", xlab="Allele frequencies",
     main="Distribution of allele frequencies", nclass=20)
temp <- density(myFreq, bw=.05)
lines(temp$x, temp$y*2,lwd=3)

pca1 <- glPca(Eryth20_l)
scatter(pca1, posi="bottomleft")

tre <- nj(dist(as.matrix(Eryth20_l)))
tre
plot(tre, typ="fan", cex=0.7)

myCol <- colorplot(pca1$scores,pca1$scores, transp=TRUE, cex=4)
abline(h=0,v=0, col="grey")
plot(tre, typ="fan", show.tip=FALSE)
tiplabels(pch=20, col=myCol, cex=4)
title("NJ tree ")

dapc1 <- dapc(Eryth20_l, n.pca=10, n.da=2)
scatter(dapc1,scree.da=FALSE, bg="white", posi.pca="topright", legend=TRUE,
        txt.leg=paste("group", 1:7))#, col=c("red","blue"))

myCol <- c("darkblue","purple","green","orange","red","blue","darkgreen")

compoplot(dapc1, col=myCol,lab="", txt.leg=paste("group", 1:7), ncol=2) ;popNames(Eryth20_l)
# cottia et villosa sont peut-être encore très proches

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

ploidy(Eryth20_l) = 2
glPlot(Eryth20_l, posi="topleft")

#calcul de distance
Eryth20_l.pop.dist <- poppr::bitwise.dist(Eryth20_l)
Eryth20_l.pop.dist


# dapc ##
grp = find.clusters(Eryth20_l, max.n.clust = round(nInd(Eryth20_l))-1)
names(grp) ; grp$Kstat ; grp$stat ; grp$grp ; grp$size
table(pop(Eryth20_l), grp$grp) ; table.value(table(pop(Eryth20_l), grp$grp), col.lab=paste("inf", 1:6),row.lab=paste("ori", 1:6))
dapc1 <- dapc(Eryth20_l, grp$grp)
dapc1 ; scatter(dapc1)

dapc1 = dapc(Eryth20_l, var.contrib = TRUE, scale = FALSE, n.pca = length(pop(Eryth20_l)), n.da = nPop(Eryth20_l) - 1)

dapc1 <- dapc(Eryth20_l)
dapc1 ;

assignplot(dapc1, only.grp=NULL, subset=NULL, new.pred=NULL, cex.lab=.75,pch=3)

scatter(dapc1, posi.da="topleft", bg="white", pch=17:22)

myCol <- c("darkblue","purple","green","orange","red","blue","darkgreen","yellow","cyan")
scatter(dapc1, posi.da="topright", bg="white",
        pch=17:22, cstar=0, col=myCol, scree.pca=TRUE,
        posi.pca="topleft")

# PCAdapt analysis ####
#install.packages("pcadapt")

filename <- read.pcadapt("data_vcf/Eryth10_t.vcf", type = "vcf")
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


system("sed -i '1d' data_vcf/tryhard.str ") #ne marche pas car premiere ligne en trop

spider("data_vcf/Eryth20_t.vcf","VCF","data_vcf/Eryth20_t.str","STRUCTURE")
system("sed -i '1d' data_vcf/Eryth20_t.str ")

struct2geno(file = "data_vcf/Eryth20_t.str", TESS = FALSE, diploid = T, FORMAT = 2,extra.row = 0, extra.col = 2, output = "data_vcf/Eryth20_t.geno")
#permet de créer le fichier format geno 23 individuals and 25086 markers. (SNP)

obj  <- snmf("data_vcf/Eryth20_t.geno", K = 1:14, entropy = T, ploidy = 2,
             CPU = 7,repetitions = 10, project= "new", alpha=100)

source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")

struct2geno(file = "data_vcf/tryhard.str", TESS = FALSE, diploid = T, FORMAT = 2,extra.row = 0, extra.col = 2, output = "data_vcf/tryhard.geno")
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

Pop = function(K, file,ID) {
  obj.snmf = snmf(file, K = K, alpha = 100, project = "new",iterations = 2000, repetitions = 20,
                  CPU = 7)
  ce=cross.entropy(obj,K=K)
  best = which.min(ce)
  qmatrix = Q(obj.snmf, K = K, run = best)
  barplot(t(qmatrix), col = color, border = NA, space = 0,xlab = "Individuals", ylab = "Admixture coefficients",
          names.arg =ID, las = 2)}


par(mfrow = c(3,4))
for (i in 1:12) Pop(i,"data_vcf/Eryth10_t.geno",  c(a,p,c,v,h,d)) ;beep(3)

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





