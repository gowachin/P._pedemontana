library(ape)
mat=read.dna(file='data/freebayes_-F0.3-n10-m30_-q20_mincov20_90samples_SNPs_only_FASTA.fasta',format='fasta')
mat
row.names(mat)

int = c('AP1','AMB','AML','AOL','CS1','CP1','CP4','DGB','DRL','DMB','HC1','HGL','HS2','HP1','HPB','PT1','PV1','GA2','GA4','VR3','VR1','VL2','VB1')
Erythro=mat[which(row.names(mat)%in%int),] # with one outgroup: AP1 (P. lutea)
Erythro

# subset to keep SNPs only
as.character(Erythro[,1])

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
write.dna(SNP,file='data/freebayes_-F0.3-n10-m30_-q20_mincov20_subsection_Erythrodrosum_23samples_SNPs_only_FASTA.fasta',format='fasta')


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

unmatch = sort(row.names(mat)[row.names(mat) %!in% Appendix_table$Recode])
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
