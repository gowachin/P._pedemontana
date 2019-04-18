library(ape)
mat=read.dna(file='data/freebayes_-F0_3-n10-m30_-q20_mincov20_90samples_SNPs_only_FASTA.fasta',format='fasta')
mat
row.names(mat)

int = c('AP1','AMB','AML','AOL','CS1','CP1','CP4','DGB','DRL','DMB','HC1','HGL','HS2','HP1','HPB','PT1','PV1','GA2','GA4','VR3','VR1','VL2','VB1')
Erythro=mat[which(row.names(mat)%in%int),] # with one outgroup: AP1 (P. lutea)
Erythro

# subset to keep SNPs only
as.character(Erythro[,1])
Erythro[,1]

missing_bases=rep(NA,dim(Erythro)[2])
for (i in 1:length(missing_bases)){
  temp=c(unique(as.character(Erythro[,i])))
  missing_bases[i]=sum(temp=='?') # attention, totalement stringent!!
  #print(i)
}
table(missing_bases)
informativ=Erythro[,which(missing_bases == 0)]
informativ

N_bases=rep(NA,dim(informativ)[2])
for (i in 1:length(N_bases)){
  temp=c(unique(as.character(informativ[,i])))
  N_bases[i]=sum(temp!='?')
  #print(i)
}
table(N_bases)
#composition des sites :
#N_bases
# 1    2    3    4    5
# 3684 1457  319   12    1    ### donne environ 353 SNP binaires car sites ambigues
Inf_SNP=informativ[,which(N_bases>1)]
Inf_SNP
write.dna(informativ,file='data/freebayes_-F0_3-n10-m30_-q20_mincov20_subsection_Erythrodrosum_23samples_SNPs_informativ_5473_FASTA.fasta',format='fasta')

write.dna(Inf_SNP,file='data/freebayes_-F0_3-n10-m30_-q20_mincov20_subsection_Erythrodrosum_23samples_SNPs_informativ_1789_FASTA.fasta',format='fasta')


N_bases=rep(NA,dim(Erythro)[2])
for (i in 1:length(N_bases)){
  temp=c(unique(as.character(Erythro[,i])))
  N_bases[i]=sum(temp!='?')
  #print(i)
}

table(N_bases)
as.character(Erythro[,which(N_bases==7)])
SNP=Erythro[,which(N_bases>1)]
SNP
write.dna(SNP,file='data/freebayes_-F0_3-n10-m30_-q20_mincov20_subsection_Erythrodrosum_23samples_SNPs_only_FASTA.fasta',format='fasta')


# sous-groupe 'pedemontana s.l.'
pedemontana=c('AMB','AML','AOL','CS1','CP1','CP4','PT1','PV1','GA2','GA4')

# importation de l'appendice ####
library(readr)
Appendix_table <- read_csv("~/Bureau/Stage.l/Pedemontana/data/Appendix_table.csv")
#View(Appendix_table)
sort(Appendix_table$Code) ; length(Appendix_table$Code)
sort(row.names(mat)) ; length(row.names(mat))

# recodes les noms pour extraire infos de position
Appendix_table$Recode = Appendix_table$Code
unmatch = sort(row.names(mat)[row.names(mat) %.in% Appendix_table$Recode])
unmatch ; length(unmatch)
# PVE == Primula veris, génome de référence
for (i in 1:length(unmatch)) {
  Appendix_table[Appendix_table$Code == paste(unmatch[i],"'",sep ="") ,12] = unmatch[i]
}
Appendix_table$Recode


Occurences_Alps <- read_csv("data_carto/Occurences_Alps.csv")
#View(Occurences_Alps)
Occurences_Alps = Occurences_Alps[,-7]

Repart = list(Repart.pedemontana = Occurences_Alps[Occurences_Alps$Taxon == "Primula pedemontana",],
Repart.villosa = Occurences_Alps[Occurences_Alps$Taxon == "Primula villosa",],
Repart.hirsuta = Occurences_Alps[Occurences_Alps$Taxon == "Primula hirsuta",],
Repart.daonensis = Occurences_Alps[Occurences_Alps$Taxon == "Primula daonensis",])

summary(Repart[[3]])
par(mfrow = c(8,2))
for (i in 1:length(Repart)) {print(summary(Repart[[i]]))
  lat = hist(Repart[[i]]$lat) ; lon = hist(Repart[[i]]$lon)
  plot(lat, xlim = c(44,48), col = color[i])
  plot(lon, xlim = c(5.7,15.9), col = color[i])
  }
par(mfrow = c(1,1))
Sub_appendix = Appendix_table[Appendix_table$Recode %in% int,]
Alpine.ext = extent (5.5,15.9,44.1,47.8)
plot.raster("data_carto/GRASSLANDS_10min.tif",Alpine.ext,line = T)
tail = c(0.1,0.1,0.1,0.1) ; color = c("red","black","grey50","purple")
for (i in 1:length(Repart)) {points(x = Repart[[i]]$lon , y = Repart[[i]]$lat, cex = tail[i], col=color[i]) }
text(x = Sub_appendix$Latitude, y = Sub_appendix$Longitude, Sub_appendix$Recode, cex = 0.7) #text(x = 22, y = 30, "taille = 22 : seuil de 5%", col = "red", cex = 1.5)
#points(x = Sub_appendix$Latitude, y = Sub_appendix$Longitude)
abline(h=45.6, col = "red", lty = 4)


Pedemonta.ext = extent (6.1,9.3,44.5,46.8)
Pede_appendix = Appendix_table[Appendix_table$Recode %in% pedemontana,]
plot.raster("data_carto/GRASSLANDS_10min.tif",Pedemonta.ext,line = T)
tail = c(0.1,0.1,0.1,0.1) ; color = c("red","black","grey50","purple")
for (i in 1:length(Repart)) {points(x = Repart[[i]]$lon , y = Repart[[i]]$lat, cex = tail[i], col=color[i]) }
text(x = Pede_appendix$Latitude, y = Pede_appendix$Longitude, Pede_appendix$Recode, cex = 0.7, col = "blue") #text(x = 22, y = 30, "taille = 22 : seuil de 5%", col = "red", cex = 1.5)
#points(x = Sub_appendix$Latitude, y = Sub_appendix$Longitude)

#### script Rstructure ####
library(LEA)
library(ade4)
library(vegan)

source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")
#format structure créé par PGDSpider!
system("cd ; echo ########fichier home######## ; ls ;
        cd Téléchargements/PGDSpider_2.1.1.5/ ; echo fichier PGDSpider2 ; ls ;
       ./PGDSpider2.sh")

struct2geno(file = "data/freebayes_-F0_3-n10-m30_-q20_mincov20_subsection_Erythrodrosum_23samples_SNPs_only_str_rearranged.str",
            TESS = FALSE, diploid = TRUE, FORMAT = 2,
            extra.row = 0, extra.col = 2, output = "data/freebayes_23_str_rearranged.geno")
#permet de créer le fichier format geno 23 individuals and 25086 markers. (SNP)

obj  <- snmf("data/freebayes_23_str_rearranged.geno", K = 1:14, entropy = T, ploidy = 2,
             CPU = 7,repetitions = 10, project= "new", alpha=100)
# Choix du K optimal (20 runs)
ID =rearanged$ind.names

par(mfrow = c(1,1))
plot(obj, col = "blue", pch=1,cex=0.8)
beep(3)
color = c("orange","violet","lightgreen","red","blue","green","cyan","grey","black","yellow","darkgreen")

obj.snmf = snmf("data/freebayes_23_str_rearranged.geno", K = 11, alpha = 100, project = "new")
qmatrix = Q(obj.snmf, K = 11)
barplot(t(qmatrix), col = color, border = NA, space = 0, xlab = "Individuals", ylab = "Admixture coefficients",
        names.arg =ID, las = 2)

Pop = function(K) {
obj.snmf = snmf("data/freebayes_23_str_rearranged.geno", K = K, alpha = 100, project = "new",
                CPU = 7)
qmatrix = Q(obj.snmf, K = K)
barplot(t(qmatrix), col = color, border = NA, space = 0,xlab = "Individuals", ylab = "Admixture coefficients",
        names.arg =ID, las = 2)}

par(mfrow = c(3,4))
for (i in 1:12) Pop(i) ;beep(3)
# adegenet ####
library(parallel)


adegenetTutorial("dapc")

rearanged <- fasta2genlight("data/freebayes_-F0_3-n10-m30_-q20_mincov20_subsection_Erythrodrosum_23samples_SNPs_only_rearanged.fas", chunk=10, parallel=FALSE)
rearanged
glPlot(rearanged, posi="topleft")

informativ <- fasta2genlight('data/freebayes_-F0_3-n10-m30_-q20_mincov20_subsection_Erythrodrosum_23samples_SNPs_informativ_5473_FASTA.fasta', chunk=100, n.cores = NULL)
informativ
glPlot(informativ , posi="topleft")

genet.fullmat <- fasta2genlight("data/Erythrodrosum_22samples_MAP13_COV10_fullmatrix.fasta", chunk=10, n.cores = NULL)
genet.fullmat
glPlot(genet.fullmat, posi="topleft")





# td Stephanie ####

#1561_SNP.str : Fichier au format structure contenant les données génétiques (1561 loci) pour chaque individu (88ind)
#genotype.lfmm : Fichier au format lfmm contenant les données génétiques (1561 loci) pour chaque individu
#stations.txt : Fichier contenant les caractéristiques des stations : coordonnées GPS, altitude et données climatiques
#data.env.txt : Fichier contenant les mêmes caractéristiques environnementales mais pour chaque individu
#loci.txt : Fichier contenant les noms des loci

# Convertir le fichier format lfmm en objet geno :
genotype.geno<-lfmm2geno("genotype.lfmm", force = T )
stations<-read.table("stations.txt",h=T)
env.data<-read.table("data.env.txt",h=T)

####ACP sur les données environnementales (ade4) ####
# Sélectionner uniquement les variables altitude et climatiques
data.pca<-stations[,4:8]
row.names(data.pca)<-stations$Population
# Réaliser l’ACP
env.pca<-dudi.pca(data.pca, scannf = FALSE, nf = 2)
# Pourcentages d’inertie et contributions absolues
inertia.dudi(env.pca,col.inertia=TRUE) #2axes pour inertie de 89%
# Représentation graphique
par(mfrow=c(1,1))
plot(env.pca$li, ylab="Altitude et temperature min",xlab="Longitude et temperature max",main="Segregation des pop selon tendances environnementales")
text(env.pca$li, as.character(rownames(env.pca$li)), offset = 0.7 ,pos =1, col="darkgreen") ;
s.corcircle(env.pca$co)
s.arrow(2*env.pca$co,xax=1,yax=2,sub="axe1/axe2")
s.label(env.pca$li,xax=1,yax=2,sub="axe1/axe2", add.plot = T)
#reunion differente d'un point de vue env de la reunion

####Inférence de la structure génétique (LEA) ####
# Calcul de la structure pour plusieurs K : ici 1 à 14
obj  <- snmf("genotype.geno", K = 1:14, entropy = T, repetitions = 10, project= "new", alpha=100)
obj  <- snmf("genotype.geno", K = 1:14, entropy = T, repetitions = 1, project= "new", alpha=100)

# Choix du K optimal (20 runs)
plot(obj, col = "blue", pch=1,cex=0.8,lines=TRUE)

# Récupérer la valeur de cross-entropy pour chaque run pour le meilleur K
ce2=cross.entropy(obj,K=2)
ce3=cross.entropy(obj,K=3)
ce4=cross.entropy(obj,K=4)

# Plot de la structure du meilleur run pour le meilleur K identifié
par(mfrow=c(1,2))
#K=2
best2 = which.min(ce2)
qmatrix2=Q(obj,K=2, run=best2)
env = read.table("data.env.txt", h=T)
rownames(qmatrix2) = env$Population
barplot(t(qmatrix2),col=c(2,3), xlab = "Individuals", ylab = "Admixture coefficients",las=2,cex.names=0.6)

#K=3
best3 = which.min(ce3)
qmatrix3=Q(obj,K=3, run=best3)
env = read.table("data.env.txt", h=T)
rownames(qmatrix3) = env$Population
barplot(t(qmatrix3),col=c(2,3,4), xlab = "Individuals", ylab = "Admixture coefficients",las=2,cex.names=0.6)

#K=4
best4 = which.min(ce4)
qmatrix4=Q(obj,K=4, run=best3)
env = read.table("data.env.txt", h=T)
rownames(qmatrix4) = env$Population
barplot(t(qmatrix4),col=c(2,3,4,5), xlab = "Individuals", ylab = "Admixture coefficients",las=2,cex.names=0.6)

####Recherche de loci outlier pour la variable Min_Temp_Cold (LEA) ####

# Création d’un projet qui ne contient que cette variable
attach(env.data)
write.env(Min_Temp_Cold ,"Min_Temps_Cold.env")
detach(env.data)

# Corrélation entre données génotypiques et climatique
mint.cor<-lfmm("genotype.lfmm","Min_Temps_Cold.env",K=2,repetitions=1,project="new")

# Calculer les p-valeurs ajustées
zs = z.scores(Min_Temp_Cold_project, K = 2)
zs.median = apply(zs, MARGIN = 1, median) #Combine z-scores using the median
lambda = median(zs.median^2)/qchisq(0.5, df = 1) #Recalibrate
adj.p.values = pchisq(zs.median^2/lambda, df = 1, lower = FALSE) #Compute adjusted p-values from the combined z-scores
p=adj.p.values
hist(adj.p.values, col = "grey",nclass=20,xlab="Adjusted p-values")
candidates.qv.cold.lfmm = which(qvalue(adj.p.values, fdr = .05)$signif)

candi.cold.lfmm=as.character(read.table("loci.txt")[candidates.qv.cold.lfmm,1])
write.table(candi.cold.lfmm,"candidates_cold.lfmm.txt")

# Manhattan plot
plot(-log10(p), xlab= "loci number", ylab= "-log10pvalues",cex=0.8,pch=20,main="Manhattan plot Minimum Temperature of the Coldest month")
points(candidates.qv.cold.lfmm, -log10(adj.p.values[candidates.qv.cold.lfmm]), pch = 19, col = "orange")
text(candidates.qv.cold.lfmm,-log10(adj.p.values[candidates.qv.cold.lfmm]),as.character(candi.cold.lfmm),cex=0.7, pos=1)

####Analyse de redondance####

rda(genotype~environnement+ Conditions(geographie))
anova(rda)

####Outflank ####
# Données
data = read.table("genotype.lfmm")
popnames = pop(genind.data)
locinames = read.table("loci.txt")
# Calcul du FST pour chaque locus
FstDataFrame <- MakeDiploidFSTMat(data, locinames, popnames)
FST<-as.numeric(FstDataFrame$FST)
# Scan FST avec OutFLANK
output<-OutFLANK(FstDataFrame, LeftTrimFraction = 0.05, RightTrimFraction = 0.05, Hmin=0.1, NumberOfSamples=length(levels(pop names)), qthreshold = 0.05)
results<-output$results
# Distribution p-valeurs et Manhattan plot
hist(results$pvalues)
plot(-log(results$pvalues),pch=20,cex=0.4,xlab="Loci",ylab="-log10(pvalues)")
high_FST_outlier<-output$numberHighFstOutliers
# Distribution des FST
OutFLANKResultsPlotter(output, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.005)
