#install.packages("vcfR")
library(vcfR)
library(readr)

#vcf.cov10.fullmat.SNP = read.vcfR("data_vcf/freebayes_-F0_3-n10-m13_-q20_mincov10_90samples_SNPs_only.vcf", checkFile = F)
#vcf.cov20.fullmat.SNP = read.vcfR("data_vcf/freebayes_-F0_3-n10-m30_-q20_mincov20_90samples_SNPs_only.vcf", checkFile = F)

mincov10_90samples_CSV <- read_delim("data_vcf/freebayes_-F0_3-n10-m13_-q20_mincov10_90samples_SNPs_only CSV.csv",
                                     "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(mincov10_90samples_CSV)
mincov20_90samples_CSV <- read_delim("data_vcf/freebayes_-F0_3-n10-m30_-q20_mincov20_90samples_SNPs_only CSV.csv",
                                     "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(mincov20_90samples_CSV)

colnames(mincov20_90samples_CSV) == colnames(mincov10_90samples_CSV) # mêmes colnames pour deux tableaux
vcf = colnames(mincov10_90samples_CSV)[1:9]

Eryth = c(vcf,'AP1','AMB','AML','AOL','CS1','CP1','CP4','DGB','DRL','DMB','HC1','HGL','HS2','HP1','HPB','PT1','PV1','GA2','GA4','VR3','VR1','VL2','VB1')

Erythro_mincov10 = mincov10_90samples_CSV [,which(colnames(mincov10_90samples_CSV) %in% Eryth)] # with one outgroup: AP1 (P. lutea)

Erythro_mincov20 = mincov20_90samples_CSV [,which(colnames(mincov20_90samples_CSV) %in% Eryth)] # with one outgroup: AP1 (P. lutea)


#permet de virer lignes vides ####
clean = function(data) {
  exterminate = apply(data[,10:32],1,paste,collapse ="")
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

#permet de virer les lignes avec seulement un certain nombre de variants ####
tri = function(data,n,quiet = F, r= F,p = F) { # data est un data frame type inform_mincov10 , n = pourcentage de présence du snp , quiet est le rendu
  col = dim(data)[2]
  matrix = data[,10:col]
  matrix = matrix == "."
  row  = apply(matrix,1,sum)
  row = 1-row/(col-9) >= n
  data = data[which(row == T),]

  resum = c(pourc.end = sum(row)/length(row)*100 ,
            n.site = sum(row))

  if(quiet == F) print(resum)
  if (r== T & p == F) return(data)
  if (r== T & p == T) return(resum)
}


tri(data = inform_mincov10,n=0.8, r = F)

tri(data = inform_mincov20,n=0.8, r = F)

par(mfrow = c(2,1))
perc = seq(0,1,by = 0.1)
res = matrix(perc,nrow = length(perc), ncol = 3) ; colnames(res) = c("stringence","info.garde","n.sites")
for( i in 1:length(perc)) {
  res[i,2:3] = tri(data = inform_mincov10,n=perc[i], r = T,p=T)
}
plot(x = res[,1], y = res[,2], xlab = "seuil de stringence", ylab = "info gardée")
abline(h = 5, col = "red") ; text(x = res[,1], y = res[,2], res[,3], pos = 3)

perc = seq(0,1,by = 0.1)
res = matrix(perc,nrow = length(perc), ncol = 3) ; colnames(res) = c("stringence","info.garde","n.sites")
for( i in 1:length(perc)) {
  res[i,2:3] = tri(data = inform_mincov20,n=perc[i], r = T,p=T)
}
plot(x = res[,1], y = res[,2], xlab = "seuil de stringence", ylab = "info gardée")
abline(h = 5, col = "red") ; text(x = res[,1], y = res[,2], res[,3], pos = 3)
par(mfrow = c(1,1))

system("cd data_vcf ; ls ;
       sed -i '1 i\##fileformat=VCFv4.0' 'freebayes_-F0_3-n10-m13_-q20_mincov10_23samples_SNPs_only' " )

tryhard = tri(data = inform_mincov20,n=1, r = T)
plot(log(tryhard$QUAL)) ; abline(h = log(0.3), col = "red")
dim(tryhard[which(tryhard$QUAL > 20),])
dim(tryhard)

write.csv(tryhard, file = "tryhardtest.csv")
