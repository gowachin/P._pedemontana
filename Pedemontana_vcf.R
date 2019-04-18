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

mincov10_Eryth_CSV <- readr::read_delim("data_vcf/freebayes_-F0_3-n10-m13_-q20_mincov10_Eryth_SNPs_onlyCSV.csv","\t", escape_double = FALSE, trim_ws = TRUE)
#mincov20_Eryth_CSV <- readr::read_delim("data_vcf/freebayes_-F0_3-n10-m30_-q20_mincov20_Eryth_SNPs_onlyCSV.csv","\t", escape_double = FALSE, trim_ws = TRUE)

a = c('AMB','AML','AOL') #apenina
c = c('CS1','CP1','CP4') #cottia
d = c('DGB','DRL') #daonensis
h = c('DMB','HC1','HGL','HS2','HP1','HPB') #hirsuta
p = c('PT1','PV1','GA2','GA4') #pedemontana
v = c('VR3','VR1','VL2','VB1') #villosa

Pedemo10 = subset_reorder(mincov10_Eryth_CSV, c(a[-2],p,c,h)) ; colnames(Pedemo10)

Pedemo10_5r = rare(Pedemo10, rare  = 0.05, r= T) ; colnames(Pedemo10_5r)

levels = apply(Pedemo10[,-c(1:9)],1, function(data) levels(as.factor(substr(as.character(data),1,3)))[!levels(as.factor(substr(as.character(data),1,3))) %in% "."])
summary(as.factor(unlist(levels)))

levels = apply(Pedemo10_5r[,-c(1:9)],1, function(data) levels(as.factor(substr(as.character(data),1,3)))[!levels(as.factor(substr(as.character(data),1,3))) %in% "."])
summary(as.factor(unlist(levels)))

#tri selon QUAL in vcf
plot(log(Pedemo10_5r$QUAL), cex = 0.1, ylim = c(log(0.0000001),20)) # virer des points avec pas assez bon Phred??? qu'est-ce qu'il représente??
abline(h = log(20), col = "green")
Pedemo10_5r = Pedemo10_5r[which(Pedemo10_5r$QUAL >= 20),]
dim(Pedemo10_5r[which(Pedemo10_5r$QUAL >= 20),])

#permet de virer les lignes avec seulement un certain nombre de variants ####
#fonction tri dans le package
Pedemo10_5r_9s_8i = tri(data = Pedemo10_5r,n.r=0.9,n.c = 0.8, r = T) # avec ces seuils on vire l'individu AML car trop d'info manquantes pour cet individu
colnames(Pedemo10_5r_9s_8i[,-c(1:9)])

save2vcf(Pedemo10_5r_9s_8i)

pop = c("apennina", "apennina"#, "apennina"
        ,"pedemontana","pedemontana"
        ,"valgau","valgau"
        ,"cottia","cottia","cottia"
        ,"hirsuta","hirsuta","hirsuta","hirsuta","hirsuta","hirsuta"
)

# analyse pour genind ####
Pedemo10_5r_9s_8i_v = read.vcfR("data_vcf/Pedemo10_5r_9s_8i.vcf", checkFile = T) ; Pedemo10_5r_9s_8i_v
Pedemo10_5r_9s_8i_g = vcfR2genind(Pedemo10_5r_9s_8i_v) ; Pedemo10_5r_9s_8i_g
row.names(Pedemo10_5r_9s_8i_g$tab)

#assign.pop = function(genind) {
#  for (i in 1:length(row.names(Pedemo10_5r_9s_8i_g$tab))){
#  }
#  pop(Pedemo10_5r_9s_8i_g) <- as.factor(pop)
#  return(genind)
#}

pop(Pedemo10_5r_9s_8i_g) <- as.factor(pop)

toto <- summary(Pedemo10_5r_9s_8i_g)
names(toto)
par(mfrow=c(2,1))
barplot(toto$loc.n.all, ylab="Number of alleles",
        main="Number of alleles per locus")
barplot(toto$Hexp-toto$Hobs, main="Heterozygosity: expected-observed",
        ylab="Hexp - Hobs")

#Is mean observed H significantly lower than mean expected H ? Nope
bartlett.test(list(toto$Hexp,toto$Hobs))
t.test(toto$Hexp,toto$Hobs,pair=T,var.equal=TRUE,alter="greater")

Pedemo10_5r_9s_8i_g_hwt <- hw.test(Pedemo10_5r_9s_8i_g, B=0) # KEZAKO??

fstat(Pedemo10_5r_9s_8i_g)
#     pop         Ind
#Total 0.04777743 -0.2038086
#pop   0.00000000 -0.2642093

#This table provides the three F statistics F st (pop/total), F it (Ind/total), and F is
#(ind/pop). These are overall measures which take into account all genotypes and all loci.

# excés d'He dans les espèces, donc on peut supposer écart a HW

matFst <- pairwise.fst(Pedemo10_5r_9s_8i_g)
matFst
pop(Pedemo10_5r_9s_8i_g)
plot(nj(matFst))
#        1          2          3
# 2 0.1294091
# 3 0.1444321 0.1160444
# 4 0.1075123 0.1013335 0.1105519

par(mfrow =c(1,1))
X <- tab(Pedemo10_5r_9s_8i_g, freq = TRUE, NA.method = "mean")
pca1 <- dudi.pca(X, scale = FALSE, scannf = FALSE, nf = 3)
barplot(pca1$eig[1:20], main = "PCA eigenvalues", col = heat.colors(50))

s.class(pca1$li, pop(Pedemo10_5r_9s_8i_g))
title("PCA of Primula west clade")
add.scatter.eig(pca1$eig[1:5],posi ="bottomright", 3,1,2)

col <- funky(15)
s.class(pca1$li, pop(Pedemo10_5r_9s_8i_g),xax=1,yax=2, col=transp(col,.6), axesell=FALSE,
        cstar=0, cpoint=3, grid=FALSE)

grp = find.clusters(Pedemo10_5r_9s_8i_g, max.n.clust = round(nInd(Pedemo10_5r_9s_8i_g))-1)#, n.pca = length(pop(Eryth20_l)))
names(grp) ; grp$Kstat ; grp$stat ; grp$grp ; grp$size
table(pop(Pedemo10_5r_9s_8i_g), grp$grp) ; table.value(table(pop(Pedemo10_5r_9s_8i_g), grp$grp), col.lab=paste("inf", 1:6),row.lab=paste("ori", 1:6))
dapc1 <- dapc(Pedemo10_5r_9s_8i_g, grp$grp)
dapc1 ; scatter(dapc1)

dapc1 = dapc(Pedemo10_5r_9s_8i_g, var.contrib = TRUE, scale = FALSE, n.pca = length(pop(Pedemo10_5r_9s_8i_g)), n.da = nPop(Pedemo10_5r_9s_8i_g) - 1)

dapc1 <- dapc(Pedemo10_5r_9s_8i_g)
dapc1 ;

assignplot(dapc1, only.grp=NULL, subset=NULL, new.pred=NULL, cex.lab=.75,pch=3)

scatter(dapc1, posi.da="topleft", bg="white", pch=17:22)

myCol <- c("darkblue","purple","green","orange","red","blue","darkgreen","yellow","cyan")
scatter(dapc1, posi.da="topright", bg="white",
        pch=17:22, cstar=0, col=myCol, scree.pca=TRUE,
        posi.pca="topleft")


# LEA analysis ####
# PGDSpider
spider = function(input,inFORM,output,outFORM) {
  if (inFORM == "VCF" & outFORM == "STRUCTURE") {spid = "Spid_VCF_STRUCTURE.spid"}
  if (inFORM == "VCF" & outFORM == "PED") {spid = "Spid_VCF_PED.spid"}

  command = paste("cd ; java -Xmx1024m -Xms512M -jar Téléchargements/PGDSpider_2.1.1.5/PGDSpider2-cli.jar -inputfile Bureau/BEE/Stage/Pedemontana/",
                  input, " -inputformat ", inFORM ," -outputfile Bureau/BEE/Stage/Pedemontana/", output, " -outputformat ", outFORM ,
                  " -spid Bureau/BEE/Stage/Pedemontana/data_vcf/",spid,sep = "")
  print(command)
  system(command)
}
source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")

spider("data_vcf/Pedemo10_5r_9s_8i.vcf","VCF","data_vcf/Pedemo10_5r_9s_8i.str","STRUCTURE")
system("sed -i '1d' data_vcf/Pedemo10_5r_9s_8i.str ") #ne marche pas car premiere ligne en trop

struct2geno(file = "data_vcf/Pedemo10_5r_9s_8i.str", TESS = FALSE, diploid = T, FORMAT = 2,extra.row = 0, extra.col = 2, output = "data_vcf/Pedemo10_5r_9s_8i.geno")

obj  <- snmf("data_vcf/Pedemo10_5r_9s_8i.geno", K = 1:8, entropy = T, ploidy = 2,
             CPU = 7,repetitions = 10, project= "new", alpha=100)

# Choix du K optimal (20 runs)
ID = row.names(Pedemo10_5r_9s_8i_g$tab)

par(mfrow = c(1,1))
plot(obj, col = "blue", pch=1,cex=0.5)
beep(3)
color = c("orange","violet","lightgreen","red","blue","green","grey","black")

obj.snmf = snmf("data_vcf/Pedemo10_5r_9s_8i.geno", K = K, alpha = 100, project = "new", repetition = 10)
ce=cross.entropy(obj,K=K)
best = which.min(ce)
qmatrix = Q(obj.snmf, K = K, run = best)
barplot(t(qmatrix), col = color, border = NA, space = 0, xlab = "Individuals", ylab = "Admixture coefficients",
        names.arg =ID, las = 2)

Pop = function(K, file,ID) {
  obj.snmf = snmf(file, K = K, alpha = 100, project = "new",iterations = 2000, repetitions = 10,
                  CPU = 7)
  ce=cross.entropy(obj,K=K)
  best = which.min(ce)
  qmatrix = Q(obj.snmf, K = K, run = best)
  barplot(t(qmatrix), col = color, border = NA, space = 0,xlab = "Individuals", ylab = "Admixture coefficients",
          names.arg =ID, las = 2)}


par(mfrow = c(3,4))
for (i in 1:11) Pop(i,"data_vcf/Pedemo10_5r_9s_8i.geno",  ID ) ;beep(3)

