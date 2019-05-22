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
library(pbapply)
library(raster)

out = c("AP1")
a = c('AMB','AOL') #apenina  j'ai enlevé AML
c = c('CS1','CP1','CP4') #cottia
d = c('DGB','DRL') #daonensis
h = c('DMB','HC1','HGL','HS2','HP1','HPB') #hirsuta
p = c('PT1','PV1') #pedemontana
va = c('GA2','GA4')
v = c('VR3','VR1','VL2','VB1') #villosa

Eryth.file = dataset(ind= c(a,p,va,c,v,h,d)
                     ,popfile= "Populations.csv"
                     ,entryfile= "data_vcf/freebayes_-F0_3-n10-m13_-q20_mincov10_Eryth_SNPs.vcf"
                     ,name = "data_vcf/Eryth"
                     ,rare= 0.05,qual= 20,missLoci= 0.95,missInd= 0,LD= 1e4)
beep(3)


PedeHirsu.file = dataset(ind= c(a,c,p,va,h)
                         ,popfile= "Populations.csv"
                         ,entryfile= "data_vcf/freebayes_-F0_3-n10-m13_-q20_mincov10_Eryth_SNPs.vcf"
                         ,name = "data_vcf/Hirsuta"
                         ,rare= 0.05,qual= 20,missLoci= 0.95,missInd= 0,LD= 1e4)
beep(3)


Pede.file = dataset(ind= c(a,p,c,h)
                         ,popfile= "Populations.csv"
                         ,entryfile= "data_vcf/freebayes_-F0_3-n10-m13_-q20_mincov10_Eryth_SNPs.vcf"
                         ,name = "data_vcf/Pedemontana"
                         ,rare= 0.05,qual= 20,missLoci= 0.95,missInd= 0.8,LD= 1e4)
beep(3)

# adegenet ####

file = Eryth.file
file = Pede.file
file = PedeHirsu.file

# analyse pour genind ####
VCF = read.vcfR(file$.vcf, checkFile = T) ; VCF
GENIND = vcfR2genind(VCF) ; GENIND
row.names(GENIND$tab)
GENLIGHT = vcfR2genlight(VCF, n.cores = 7)
pop(GENLIGHT) <- as.factor(file$.pop)
pop(GENIND) <- as.factor(file$.pop)



fstat(GENIND)
#     pop         Ind
#Total 0.1242216 -0.09033675
#pop   0.0000000 -0.24499159

#This table provides the three F statistics F st (pop/total), F it (Ind/total), and F is
#(ind/pop). These are overall measures which take into account all genotypes and all loci.

# excés d'He dans les espèces, donc on peut supposer écart a HW

matFst <- pairwise.fst(GENIND, res.type = "matrix")
matFst
pop(GENIND)
plot(nj(matFst))
plot(root(nj(matFst), 7),main = expression("Neighbor joining d'Erythrodrosum sur les " *F[st]* " par pair.") )

par(mfrow =c(1,1))
X <- tab(GENIND, freq = TRUE, NA.method = "mean")
pca1 <- dudi.pca(X, scale = FALSE, scannf = FALSE, nf = 3)
barplot(pca1$eig[1:20], main = "PCA eigenvalues", col = heat.colors(50))

s.class(pca1$li, pop(GENIND))
title("ACP d'Erythrodrosum")
add.scatter.eig(pca1$eig[1:20], 3,1,2)

col <- funky(15)
s.class(pca1$li, pop(GENIND),xax=1,yax=2, col=transp(col,.6), axesell=FALSE,
        cstar=0, cpoint=3, grid=FALSE)

grp = find.clusters(GENIND, max.n.clust = round(nInd(GENIND))-1, n.pca = length(pop(GENIND)))
names(grp) ; grp$Kstat ; grp$stat ; grp$grp ; grp$size
table(pop(GENIND), grp$grp) ; table.value(table(pop(GENIND), grp$grp), col.lab=c("pedemontana s.l.","est-alpin","hirsuta"),row.lab=levels(file$.pop))
dapc1 <- dapc(GENIND, grp$grp, n.pca = length(pop(GENIND))-1)
dapc1 ; scatter(dapc1)

plot(grp$Kstat, type ="b", col = "blue", main = "Value of BIC versus number of clusters", xlab= "number of clusters", ylab = "BIC")
abline(v=3,col= "darkgreen", lty = 3)

dapc1 = dapc(GENIND, var.contrib = TRUE, scale = FALSE, n.pca = length(pop(GENIND)), n.da = nPop(GENIND) - 1)


dapc1 <- dapc(GENIND)
scatter(dapc1)

assignplot(dapc1, only.grp=NULL, subset=NULL, new.pred=NULL, cex.lab=.75,pch=3)

scatter(dapc1, posi.da="topleft", bg="white", pch=17:22)

myCol <- c("darkblue","purple","green","orange","red","blue","darkgreen","yellow","cyan")
scatter(dapc1, posi.da="bottomright", bg="white",
        pch=17:22, cstar=0, col=myCol, scree.pca=TRUE,
        posi.pca="bottomleft")

# analyse pour genligth! ####
GENLIGHT = vcfR2genlight(VCF, n.cores = 7)
pop(GENLIGHT) <- as.factor(file$.pop)

myFreq <- glMean(GENLIGHT)
hist(myFreq, proba=TRUE, col="gold", xlab="Allele frequencies",
     main="Distribution of (second) allele frequencies")
temp <- density(myFreq)
lines(temp$x, temp$y*1.8,lwd=3)

myFreq <- c(myFreq, 1-myFreq)
hist(myFreq, proba=TRUE, col="darkseagreen3", xlab="Allele frequencies",
     main="Distribution of allele frequencies", nclass=20)
temp <- density(myFreq, bw=.05)
lines(temp$x, temp$y*2,lwd=3)

pca1 <- glPca(GENLIGHT)
scatter(pca1, posi="bottomleft")

tre <- nj(dist(as.matrix(GENLIGHT)))
tre
plot(root(tre,22), typ="fan", cex=0.7)

myCol <- colorplot(pca1$scores,pca1$scores, transp=TRUE, cex=4)
abline(h=0,v=0, col="grey")
plot(tre, typ="fan", show.tip=FALSE)
tiplabels(pch=20, col=myCol, cex=4)
title("NJ tree ")

dapc1 <- dapc(GENLIGHT, n.pca=10, n.da=2)
scatter(dapc1,scree.da=FALSE, bg="white", posi.pca="topright", legend=TRUE,
        txt.leg=paste("group", 1:7))#, col=c("red","blue"))

myCol <- c("darkblue","purple","green","orange","red","blue","darkgreen")

compoplot(dapc1, col=myCol,lab="", txt.leg=paste("group", 1:7), ncol=2) ;popNames(GENLIGHT)
# cottia et villosa sont peut-être encore très proches

#While a large number of loci are nearly fixed (frequencies close to 0 or 1), there is an
#appreciable number of alleles with intermediate frequencies and therefore susceptible to
#contain interesting biological signal.

#Missing data should ideally be i) not too numerous and ii) randomly distributed in the dataset. In a situation like yours,
#individuals are more precious than markers, so I would discard loci with a majority of NAs, and briefly check the structure
#of the remaining missing entries.

#NAs are basically replaced to the mean allele frequency. This means individuals with NAs will tend to
#be placed closer to the origin. Also, individuals with similar patterns of NAs will be seen as more
#similar than they probably are in reality.

#If you really have a big missing value problem, and lot of NAs you cannot discard, one possibility
#would be to get a matrix of 1 and 0 where '1' indicate NAs, and do the PCA of this.
#If you obtain a structure, then this is a sign of problem - your NAs are not randomly distributed.

ploidy(GENLIGHT) = 2
glPlot(GENLIGHT, posi="topleft")
file$.ind

#calcul de distance
GENLIGHT.pop.dist <- poppr::bitwise.dist(GENLIGHT)
GENLIGHT.pop.dist

# dapc ##
grp = find.clusters(GENLIGHT, max.n.clust = round(nInd(GENLIGHT))-1)
names(grp) ; grp$Kstat ; grp$stat ; grp$grp ; grp$size
table(pop(GENLIGHT), grp$grp) ; table.value(table(pop(GENLIGHT), grp$grp), col.lab=paste("inf", 1:6),row.lab=paste("ori", 1:6))
dapc1 <- dapc(GENLIGHT, grp$grp)
dapc1 ; scatter(dapc1)

dapc1 = dapc(GENLIGHT, var.contrib = TRUE, scale = FALSE, n.pca = length(pop(GENLIGHT)), n.da = nPop(GENLIGHT) - 1)

dapc1 <- dapc(GENLIGHT)
dapc1 ;

assignplot(dapc1, only.grp=NULL, subset=NULL, new.pred=NULL, cex.lab=.75,pch=3)

scatter(dapc1, posi.da="topleft", bg="white", pch=17:22)

myCol <- c("darkblue","purple","green","orange","red","blue","darkgreen","yellow","cyan")
scatter(dapc1, posi.da="topright", bg="white",
        pch=17:22, cstar=0, col=myCol, scree.pca=TRUE,
        posi.pca="topleft")

# PCAdapt analysis ####
#install.packages("pcadapt")

filename <- read.pcadapt(file$.vcf, type = "vcf")
x <- pcadapt(input = filename, K = 9)
plot(x, option = "screeplot")
plot(x, option = "scores", pop = file$.pop)
plot(x, option = "scores", i = 2, j = 3, pop = file$.pop)
plot(x, option = "scores", i = 3, j = 4, pop = file$.pop)

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

obj  <- snmf(file$.geno, K = 1:14, entropy = T, ploidy = 2,
             CPU = 7,repetitions = 10, project= "new", alpha=10, iterations = 1000)
# Choix du K optimal (20 runs)


par(mfrow = c(1,1))
plot(obj, col = "blue", pch=1,cex=0.5)
beep(3)


color = c("orange","violet","lightgreen","red","blue","green","cyan","grey","black","yellow","darkgreen")
color = c("chartreuse3","cadetblue2","dodgerblue3","sienna3","coral2")

par(mfrow = c(1,1))
obj.snmf = snmf(file$.geno, K = 11, alpha = 10, project = "new")
qmatrix = Q(obj.snmf, K = 11)
barplot(t(qmatrix), col = color, border = 1, space = 0.05, xlab = "Individuals", ylab = "Admixture coefficients",
        names.arg =file$.ind, las = 2)

par(mfrow = c(4,1))
for (i in 2:5) Pop(K= i,files = file$.geno, ID= file$.ind) ;beep(3)

# taille de pop ####

# pegas

DNABIN = vcfR2DNAbin(VCF, consensus = T, extract.haps = F)

rownames(DNABIN)
tajima.test(DNABIN[c(1:3,8:10),]) # cott-app
tajima.test(DNABIN[c(4:7),]) # pede
tajima.test(DNABIN[c(15:20),]) # hirsuta
tajima.test(DNABIN[c(11:15),]) # villosa
tajima.test(DNABIN)

a = c("apennina", "apennina","apennina")
c = c("cottia","cottia","cottia")
da = c("daonensis","daonensis")
hi = c("hirsuta","hirsuta","hirsuta","hirsuta","hirsuta","hirsuta")
p = c("pedemontana","pedemontana")
v = c("villosa","villosa","villosa","villosa")
va = c("valgau","valgau")

# haplotype network ####
d = Eryth_dna
e <- dist.dna(d)
h <- pegas::haplotype(d)
h <- sort(h, what = "label")
(net <- pegas::haploNet(h))
ind.hap<-with(
  stack(setNames(attr(h, "index"), rownames(h))),
  table(hap=ind, pop=rownames(d)[values])
)
plot(net, size=attr(net, "freq"), scale.ratio=0.01, pie=ind.hap2, labels = T, lwd = 1, fast = T)
legend("topleft", colnames(ind.hap), col=rainbow(ncol(ind.hap)), pch=19, ncol=2)

wrong.pop<-c(a,p,va,c,v,hi,da)
ind.hap2<-with(
  stack(setNames(attr(h, "index"), rownames(h))),
  table(hap=ind, pop=wrong.pop[values])
)
plot(net, size=attr(net, "freq"), scale.ratio = 0.02, cex = 0.8, pie=ind.hap2)
legend("topleft", colnames(ind.hap2), col=rainbow(ncol(ind.hap2)), pch=20)


# ABBA BABA compute D ####

CalcPopD <- function(alignment = "alignment.fasta"){
  ##  Now we have eqn. 2 from page 2240
  ##  input is an alignment the can take multiple sequences from each
  ##  population of interest.  IMPORTANT MAKE SURE SEQUENCES ARE IN ORDER
  ##  P1, P2, P3, OUTGROUP!  Again we find the biallelic sites but now
  ##  those biallelic sites need not be fixed and we will calculate frequencies
  ##  of SNP for each population.  The way the function is set up we do need to
  ##  feed in an alignment where each sequence from a population has the same name:
  ##  pop1
  ##  AACCACAAGCCAGCTCAGCTACAG
  ##  pop1
  ##  TACAACAAGCGAGCTCAGCTACAG
  ##  pop1
  ##  GGCCACAAGCCAGCTCAGCTACAG
  ##  pop2
  ##  GGCCACAAGCCAGCTCAGCTACAG
  ##  pop2
  ##  GGCCACAAGCCAGCTCAGCTACAG
  ##  pop3
  ##  TACCACAAGCCAGCTCAGCTACAG
  ##  OUTGROUP
  ##  TACCAGGAGCCAGCTCTTCTACCC
  Mode <- function(x) {                                                                      #  i need a little mode function which R is lacking ugh
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  alignment<-read.alignment(alignment, format="fasta")                                       #  read in the alignment
  alignment.matrix<-matrix(,length(alignment$nam),nchar(alignment$seq[[1]])+1)               #  make a matrix for the alignment
  for(i in 1:length(alignment$nam)){
    alignment.matrix[i,2:ncol(alignment.matrix)]<-unlist(strsplit(alignment$seq[[i]],""))    #  fill in the matrix
  }
  alignment.matrix[,1]<-alignment$nam                                                        #  get those names into our matrix row names dont work :(
  groups<-unique(alignment$nam)
  p1 <- p2 <- p3 <- p4 <- 0                                                                  #  lets just set up the variable names from the durand paper
  numerator <- denominator <- 0
  useful<-0                                                                                  #  plus some of my own
  segregating<-0                                                                             #  plus some of my own
  seg.pos<-F                                                                                 #  plus some of my own
  for(i in 2:ncol(alignment.matrix)){                                                        #  run through all sites
    seg.pos<-F                                                                               #  reset this switch
    if(length(unique(alignment.matrix[,i]))==2){                                             #  unique(c(p1,p2,p3,o))==2 aka biallelic
      A <- Mode(alignment.matrix[alignment.matrix[, 1] == groups[4], i])                     #  lets treat the more common variant in the outgroup as "A"
      B <- unique(alignment.matrix[,i])[unique(alignment.matrix[, i]) != A]                  #  not purposely obfuscating... the other variant in variable "B"
      if(B %in% unique(alignment.matrix[alignment.matrix[, 1] == groups[3], i])){            #  makes sure that we have at least some indication of an ABBA/BABA pattern
        if(length(unique(alignment.matrix[alignment.matrix[, 1] %in% groups[1:2], i])) == 2){  #  makes sure that we've got some different resolutions in the ingroups
          useful <- useful + 1                                                                 #  lets just keep track of how many sites are even useful
          if(length(unique(alignment.matrix[alignment.matrix[, 1] == groups[1], i])) == 2) {seg.pos<-T}#  next 5 lines are a lame way of counting sites that are segregating
          if(length(unique(alignment.matrix[alignment.matrix[, 1] == groups[2], i])) == 2) {seg.pos<-T}#  vs those that are fixed another words is population sampling
          if(length(unique(alignment.matrix[alignment.matrix[, 1] == groups[3], i])) == 2) {seg.pos<-T}#  really of any value within the data set that we are examining
          if(length(unique(alignment.matrix[alignment.matrix[, 1] == groups[4], i])) == 2) {seg.pos<-T}
          if(seg.pos == T){segregating <- segregating + 1}
          #print(segregating)
          p1 <- (sum(alignment.matrix[alignment.matrix[, 1] == groups[1], i] == A))/length(alignment.matrix[alignment.matrix[, 1] == groups[1], i])  #  freq of A snp in first population
          p2 <- (sum(alignment.matrix[alignment.matrix[, 1] == groups[2], i] == A))/length(alignment.matrix[alignment.matrix[, 1] == groups[2], i])  #  freq of A snp in second population
          p3 <- (sum(alignment.matrix[alignment.matrix[, 1] == groups[3], i] == A))/length(alignment.matrix[alignment.matrix[, 1] == groups[3], i])  #  freq of A snp in third population
          p4 <- (sum(alignment.matrix[alignment.matrix[, 1] == groups[4], i] == A))/length(alignment.matrix[alignment.matrix[, 1] == groups[4], i])  #  freq of A snp in outgroup population
          #  Durands explanation of eqn 2 is lacking... as least to my feable mind!
          #  it appears to me that as written p hat is actually the frequency of SNP "B" so....
          #  snap...  vindicated my interpretation matches that found in the supplemental material of the
          #  heliconius genome paper supplement... too cool
          p1 <- 1-p1  #convert these over from proportion A to proportion B
          p2 <- 1-p2  #convert these over from proportion A to proportion B
          p3 <- 1-p3  #convert these over from proportion A to proportion B
          p4 <- 1-p4  #convert these over from proportion A to proportion B
          numerator <- ((1 - p1) * p2 * p3 * (1 - p4)) - (p1 * (1 - p2) * p3 * (1 - p4)) + numerator              #  build up our numerator sum
          denominator <- ((1 - p1) * p2 * p3 * (1 - p4)) + (p1 * (1 - p2) * p3 * (1 - p4)) + denominator          #  build up our denominator sum
        }
      }
    }
  }
  d <- numerator / denominator    #what its all about

  user.result <- list()
  user.result$d.stat <- d
  user.result$pval <- "HELP"
  user.result$align.length <- ncol(alignment.matrix) - 1
  user.result$useful.sites <- useful
  user.result$seg.sites <- segregating
  print(paste("Sites in alignment =", ncol(alignment.matrix) - 1))
  print(paste("Number of sites with ABBA or BABA patterns =", useful))
  print(paste("Number of ABBA or BABA sites that are still segregating in at least one population =", segregating))
  print(paste("D statistic =", d))
}

library(seqinr)

CalcPopD (alignment = "data/ABBA_BABA_FASTA.fas")

install.packages("evobiR")
library(evobiR)


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

# plot carto ####
plot.raster = function(rast,ext=extent (0,0,0,0),line = F,raw =T) {
  if(raw == T) {rast=raster(rast)}
  if (ext != extent (0,0,0,0)) {plot(crop(rast, ext), col = gray.colors(20, start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL))
    } else {plot(rast, col = gray.colors(20, start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL)) }
  if(line==T){lines(read.table("data_carto/WORLD_lowres.dat"))}
}


graphics.off()
Location <- read_delim("Location.csv", "&", escape_double = FALSE, trim_ws = TRUE)
colors = c("darkolivegreen3","lightsalmon3","goldenrod2","steelblue3","slateblue2","darkorchid3","violetred2")
par(mar = c(4.85, 2.825, 4, 4.56))
barplot(100,width=1,col= "lightblue", axes = F)
par(new=TRUE, mar = c(5,4,4,2))
layout(mat = matrix(c(1), nrow = 1, ncol = 1),
       heights = c(1),
       widths = c(1))
plot(crop(raster("Europ_clim/wc2.0_bio_30s_Europ_1.tif"), extent (3.95,16,43,47.5), scalebar=F),
     col = gray.colors(20, start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL), legend=F,
     xlab = "Longitude", ylab = "Latitude")
#Occurences_Alps <- read_csv("data_carto/Occurences_Alps.csv")[,-7]
#Repart = list(Repart.pedemontana = Occurences_Alps[Occurences_Alps$Taxon == "Primula pedemontana",],
#              Repart.villosa = Occurences_Alps[Occurences_Alps$Taxon == "Primula villosa",],
#              Repart.hirsuta = Occurences_Alps[Occurences_Alps$Taxon == "Primula hirsuta",],
#              Repart.daonensis = Occurences_Alps[Occurences_Alps$Taxon == "Primula daonensis",])
#tail = c(0.1,0.1,0.1,0.1) ; color = c("red","orange","green","purple")
#for (i in 1:length(Repart)) {points(x = Repart[[i]]$lon , y = Repart[[i]]$lat, cex = tail[i], col=color[i]) }

#hirsuta = click(raster("Europ_clim/wc2.0_bio_30s_Europ_1.tif"), n=30, id=FALSE, xy=T, type="n")
hirsuta = data.frame(x = c(9.629167, 9.070833, 8.395833, 8.179167, 7.762500, 7.562500, 7.287500, 7.070833, 6.687500, 6.412500, 6.270833, 6.037500, 5.837500,5.820833, 6.079167, 6.437500, 6.612500,
                           6.712500, 7.029167, 7.045833, 6.970833, 7.104167, 7.462500, 7.837500, 8.120833, 8.379167, 8.512500, 8.870833, 9.345833,10.020833,10.29583, 10.51250, 10.51250, 10.00417),
                       y= c(47.08750,47.08750,46.99583,46.97917,46.77083,46.32083,46.23750,46.21250,45.84583,45.65417,45.51250,45.37083,45.16250,44.91250,44.56250,44.66250,44.95417,45.15417,45.37083,
                            45.55417,45.69583,45.83750,45.86250,45.94583,46.12917,46.15417,46.05417,46.01250,45.98750,46.11250,46.25417, 46.58750, 46.95417, 47.07917))
polygon(x=hirsuta$x,y = hirsuta$y,col=colors[1],density = 75,border =colors[1], lwd = 0.5)
#hirsutap = click(raster("Europ_clim/wc2.0_bio_30s_Europ_1.tif"), n=15, id=FALSE, xy=T, type="n")
hirsutap = data.frame(x = c(0.60416667, -0.12916667, -0.80416667, -0.52916667, -0.10416667,  0.98750000,  1.47916667,  1.94583333,  2.04583333,  1.35416667,  1.03750000,0.77083333),
                       y= c(42.89583, 42.98750, 42.87917, 42.68750, 42.60417, 42.44583, 42.41250, 42.40417, 42.58750, 42.83750, 42.86250, 42.87917))
polygon(x=hirsutap$x,y = hirsutap$y,col=colors[1],density = 75,border =colors[1], lwd = 0.5)
#daonensis = click(raster("Europ_clim/wc2.0_bio_30s_Europ_1.tif"), n=12, id=FALSE, xy=T, type="n")
daonensis = data.frame(x = c(10.58750, 10.43750, 10.33750, 10.36250, 10.39583, 10.36250, 10.26250, 10.39583, 10.67917, 10.87083, 10.97917, 10.81250),
                    y= c(46.61250, 46.61250, 46.53750, 46.41250, 46.29583, 46.07083, 45.81250, 45.87083, 46.07083, 46.25417, 46.42917, 46.58750))
polygon(x=daonensis$x,y = daonensis$y,col=colors[2] ,density = 75,border =colors[2], lwd = 0.5)
#villosa = click(raster("Europ_clim/wc2.0_bio_30s_Europ_1.tif"), n=12, id=FALSE, xy=T, type="n")
villosa = data.frame(x = c(14.42083, 14.13750, 13.89583, 13.57083, 13.55417, 13.87917, 14.13750, 14.37083, 14.62917, 15.02917, 15.20417, 15.07083),
                     y= c(47.07917, 47.07917, 47.16250, 47.00417, 46.82917, 46.67917, 46.81250, 46.98750, 46.87917, 46.89583, 47.07917, 47.22917))
polygon(x=villosa$x,y = villosa$y,col=colors[3],density = 75,border =colors[3], lwd = 0.5)
#cottia = click(raster("Europ_clim/wc2.0_bio_30s_Europ_1.tif"), n=8, id=FALSE, xy=T, type="n")
cottia = data.frame(x = c(6.979167, 6.845833, 6.720833, 6.720833, 6.870833, 7.087500, 7.095833, 7.079167),
                    y= c(45.05417, 44.94583, 44.86250, 44.71250, 44.57917, 44.56250, 44.73750, 44.88750))
polygon(x=cottia$x,y = cottia$y,col=colors[4],density = 75,border =colors[4], lwd = 0.5)
#apennina = click(raster("Europ_clim/wc2.0_bio_30s_Europ_1.tif"), n=10, id=FALSE, xy=T, type="n")
apennina = data.frame(x = c(9.820833,  9.629167,  9.470833,  9.537500,  9.787500, 10.070833, 10.362500, 10.587500, 10.404167, 10.070833),
                      y= c(44.58750, 44.60417, 44.57083, 44.42917, 44.32083, 44.26250, 44.14583, 44.17083, 44.30417, 44.47083))
polygon(x=apennina$x,y = apennina$y,col=colors[5],density = 75,border =colors[5], lwd = 0.5)
#pedemontana = click(raster("Europ_clim/wc2.0_bio_30s_Europ_1.tif"), n=10, id=FALSE, xy=T, type="n")
pedemontana = data.frame(x = c(7.470833, 7.054167, 6.704167, 6.737500, 6.870833, 7.054167, 7.279167, 7.495833, 7.629167, 7.629167),
                      y= c(45.67917, 45.67917, 45.49583, 45.20417, 45.11250, 45.22083, 45.37083, 45.45417, 45.56250, 45.65417))
polygon(x=pedemontana$x,y = pedemontana$y,col=colors[6],density = 75,border =colors[6], lwd = 0.5)
#ecrins = click(raster("Europ_clim/wc2.0_bio_30s_Europ_1.tif"), n=8, id=FALSE, xy=T, type="n")
ecrins = data.frame(x = c(6.437500, 6.312500, 6.129167, 6.062500, 6.104167, 6.237500, 6.404167, 6.470833),
                    y= c(44.87917, 44.94583, 44.97083, 44.89583, 44.79583, 44.72917, 44.74583, 44.80417))
polygon(x=ecrins$x,y = ecrins$y,col=colors[7],density = 75,border =colors[7], lwd = 0.5)

text(x = Location$Longitude, y = Location$Latitude, Location$Code, cex = 0.65)

points(x = c(6.1667,9.1859243,5.724524), y = c(46.2,45.4654219,45.188529), pch = 4, lwd= 2)
text(x = c(6.1667,9.1859243,5.724524), y = c(46.2,45.4654219,45.188529), c("Geneva","Milano","Grenoble"), font=c(2), pos = c(3,3,2), cex =0.8)

text(x = c(12.5), y = c(44), c("N"), font=c(2), pos = 3, cex =1.5)
arrows(12.5, 43.1, 12.5, y1 = 44, length = 0.2, angle = 30, code = 2, lwd = 5)

legend("bottomright", c("P. hirsuta (6)","P. daonensis (2)","P. villosa (4)","P. cottia (4)","P. apennina (2)","P. pedemontana (2)","P. pedemontana Ecrins (2)"),
       fill = colors, text.font=c(3), cex = 0.7)

legend("bottomleft", c("                                  "," "," "," "," "), cex = 0.86)
par(new=TRUE)
layout(mat = matrix(c( 0, 0,0,0,1,0, 0, 0,0), nrow = 3, ncol = 3),
       heights = c(5, 2,1.7),
       widths = c(0.75,1.85, 5))
barplot(100,width=1.2,col= "lightblue", axes = F)
par(new=TRUE)
layout(mat = matrix(c( 0, 0,0,0,1,0, 0, 0,0), nrow = 3, ncol = 3),
       heights = c(5, 2,1.7),
       widths = c(0.8,1.7, 5))
plot(crop(raster("Europ_clim/wc2.0_bio_30s_Europ_1.tif"), extent (-2.1,4,41.5,44),
          scalebar=F), col = gray.colors(20, start = 0.3, end = 0.9, gamma = 2.2, alpha = NULL), box=F, legend=F)
#hirsutap = click(raster("Europ_clim/wc2.0_bio_30s_Europ_1.tif"), n=15, id=FALSE, xy=T, type="n")
#hirsutap = data.frame(x = c(0.60416667, -0.12916667, -0.80416667, -0.52916667, -0.10416667,  0.98750000,  1.47916667,  1.94583333,  2.04583333,  1.35416667,  1.03750000,0.77083333),
#                      y= c(42.89583, 42.98750, 42.87917, 42.68750, 42.60417, 42.44583, 42.41250, 42.40417, 42.58750, 42.83750, 42.86250, 42.87917))
#polygon(x=hirsutap$x,y = hirsutap$y,col="red",density = 50,border ="red", lwd = 0.5)
#polygon(x=hirsutap$x,y = hirsutap$y,col="red",density = 50,border ="red", lwd = 0.5)
par(new=TRUE)
layout(mat = matrix(c(1), nrow = 1, ncol = 1),
       heights = c(1),
       widths = c(1))
par(mar = c(4.85, 2.825, 4, 4.56))
plot(100,100,type="n", axes = F, ylab="", xlab="")

text(71,80,"Pyrénées", font= 2, cex= 0.8)



