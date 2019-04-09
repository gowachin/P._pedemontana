#install.packages("vcfR")
#library(vcfR)
library(readr)

#vcf.cov10.fullmat.SNP = read.vcfR("data_vcf/freebayes_-F0_3-n10-m13_-q20_mincov10_90samples_SNPs_only.vcf", checkFile = F)
#vcf.cov20.fullmat.SNP = read.vcfR("data_vcf/freebayes_-F0_3-n10-m30_-q20_mincov20_90samples_SNPs_only.vcf", checkFile = F)

mincov10_90samples_CSV <- read_delim("data_vcf/freebayes_-F0_3-n10-m13_-q20_mincov10_90samples_SNPs_onlyCSV.csv",
                                     "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(mincov10_90samples_CSV)
mincov20_90samples_CSV <- read_delim("data_vcf/freebayes_-F0_3-n10-m30_-q20_mincov20_90samples_SNPs_onlyCSV.csv",
                                     "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(mincov20_90samples_CSV)

colnames(mincov20_90samples_CSV) == colnames(mincov10_90samples_CSV) # mêmes colnames pour deux tableaux
vcf = colnames(mincov20_90samples_CSV)[1:9]

Eryth = c(vcf,'AMB','AML','AOL','CS1','CP1','CP4','DGB','DRL','DMB','HC1','HGL','HS2','HP1','HPB','PT1','PV1','GA2','GA4','VR3','VR1','VL2','VB1')
# without the outgroup: AP1 (P. lutea)
Erythro_mincov10 = mincov10_90samples_CSV [,which(colnames(mincov10_90samples_CSV) %in% Eryth)]

Erythro_mincov20 = mincov20_90samples_CSV [,which(colnames(mincov20_90samples_CSV) %in% Eryth)]

rm(mincov10_90samples_CSV,mincov20_90samples_CSV )

#permet de virer lignes vides ####
clean = function(data) {
  n.col= dim(data)[2]
  exterminate = apply(data[,10:n.col],1,paste,collapse ="")
  compte = exterminate == paste(rep(".",23),collapse = "")
  print(summary(compte))
  resum = c(pourc.clean = (1-sum(compte)/length(compte) )*100)
  print(resum)

  end = data[which(compte == F),]
  return(end)
}

inform_mincov10 = clean(Erythro_mincov10) ; dim(inform_mincov10)
#   Mode   FALSE    TRUE
#logical  608327   16481
inform_mincov20 = clean(Erythro_mincov20) ; dim(inform_mincov20)
#Mode   FALSE    TRUE
#logical   84836     410

rm(Erythro_mincov10,Erythro_mincov20)

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

#tri(data = inform_mincov10,n.r=0.8, r = F)

test = tri(data = inform_mincov20,n.r=0.8,n.c = 0.8, r = T) # avec ces seuils on vire l'individu AML car trop d'info manquantes pour cet individu

plot(log(test$QUAL), cex = 0.1, ylim = c(log(0.001),20)) # virer des points avec pas assez bon Phred??? qu'est-ce qu'il représente??

write.table(test,"data_vcf/tryhard.csv",sep = "\t", quote = F, row.names=F)
system(" ./CSV_to_VCF.sh data_vcf/tryhard.csv ; mv data_vcf/tryhard.csv data_vcf/tryhard.vcf ; ls" )

#./CSV_to_VCF.sh file.txt

#mv tryhard.csv tryhard.vcf ; ls
-----------

# PGDSpider
system("cd ; echo ########fichier home######## ; ls ;
        cd Téléchargements/PGDSpider_2.1.1.5/ ; echo fichier PGDSpider2 ; ls ;
       ./PGDSpider2.sh")


library(LEA)
library(ade4)
library(vegan)

source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")

struct2geno(file = "data_vcf/tryhard.str", TESS = FALSE, diploid = T, FORMAT = 2,extra.row = 1, extra.col = 2, output = "data_vcf/tryhard.geno")
#ne marche pas car premiere ligne en trop
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





