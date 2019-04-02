library(ape)
mat=read.dna(file='data/freebayes_-F0.3-n10-m30_-q20_mincov20_90samples_SNPs_only_FASTA.fasta',format='fasta')
mat
row.names(mat)

int = c('AP1','AMB','AML','AOL','CS1','CP1','CP4','DGB','DRL','DMB','HC1','HGL','HS2','HP1','HPB','PT1','PV1','GA2','GA4','VR3','VR1','VL2','VB1')
mat2=mat[which(row.names(mat)%in%int),] # with one outgroup: AP1 (P. lutea)
mat2

Sub_appendix = Appendix_table[Appendix_table$Recode %in% int,]
Alpine.ext = extent (-10,25,40,50)
plot.raster("data_carto/GRASSLANDS_10min.tif",Alpine.ext)
#plot.obj("data_carto/elevation1x1_new.tif")
par(new = T)
text(x = Sub_appendix$Latitude, y = Sub_appendix$Longitude, Sub_appendix$Recode, cex = 0.7) #text(x = 22, y = 30, "taille = 22 : seuil de 5%", col = "red", cex = 1.5)
#points(x = Sub_appendix$Latitude, y = Sub_appendix$Longitude)
points(x = Occurences_Alps$lon , y = Occurences_Alps$lat, cex = 0.5) # attention problemes de noms entre lon et lat, c'est l'inverse.


# subset to keep SNPs only
as.character(mat2[,1])

N_bases=rep(NA,dim(mat2)[2])
for (i in 1:length(N_bases)){
  temp=c(unique(as.character(mat2[,i])))
  N_bases[i]=sum(temp!='?')
  #print(i)
}

table(N_bases)
as.character(mat2[,which(N_bases==7)])
mat3=mat2[,which(N_bases>1)]
mat3
write.dna(mat3,file='data/freebayes_-F0.3-n10-m30_-q20_mincov20_subsection_Erythrodrosum_23samples_SNPs_only_FASTA.fasta',format='fasta')


# sous-groupe 'pedemontana s.l.'
pedemontana=c('AMB','AML','AOL','CS1','CP1','CP4','PT1','PV1','GA2','GA4')

# importation de l'appendice ####
#library(readr)
Appendix_table <- read_csv("~/Bureau/Stage.l/Pedemontana/data/Appendix_table.csv")
#View(Appendix_table)
sort(Appendix_table$Code) ; length(Appendix_table$Code)
sort(row.names(mat)) ; length(row.names(mat))

unmatch = sort(row.names(mat)[row.names(mat) %!in% Appendix_table$Recode])
unmatch ; length(unmatch)
# PVE == Primula veris, génome de référence

# recodes les noms pour extraire infos de position
Appendix_table$Recode = Appendix_table$Code
for (i in 1:length(unmatch)) {
  Appendix_table[Appendix_table$Code == paste(unmatch[i],"'",sep ="") ,12] = unmatch[i]
}
Appendix_table$Recode


Occurences_Alps <- read_csv("data_carto/Occurences_Alps.csv")
#View(Occurences_Alps)
Occurences_Alps = Occurences_Alps[,-7]

