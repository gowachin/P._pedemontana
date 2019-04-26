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
library(introgress)
library(Pedemontana)

mincov10_Eryth_CSV <- readr::read_delim("data_vcf/freebayes_-F0_3-n10-m13_-q20_mincov10_Eryth_SNPs_onlyCSV.csv","\t", escape_double = FALSE, trim_ws = TRUE)
#mincov20_Eryth_CSV <- readr::read_delim("data_vcf/freebayes_-F0_3-n10-m30_-q20_mincov20_Eryth_SNPs_onlyCSV.csv","\t", escape_double = FALSE, trim_ws = TRUE)

a = c('AMB','AML','AOL') #apenina
c = c('CS1','CP1','CP4') #cottia
d = c('DGB','DRL') #daonensis
h = c('DMB','HC1','HGL','HS2','HP1','HPB') #hirsuta
p = c('PT1','PV1','GA2','GA4') #pedemontana
v = c('VR3','VR1','VL2','VB1') #villosa

Pedemo10 = subset_reorder(mincov10_Eryth_CSV, c(a[-2],c,p,h)) ; colnames(Pedemo10)

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
Pedemo10_5r_9s_8i = tri(data = Pedemo10_5r,n.r=1,n.c = 0.8, r = T) # avec ces seuils on vire l'individu AML car trop d'info manquantes pour cet individu
colnames(Pedemo10_5r_9s_8i[,-c(1:9)])
dim(Pedemo10_5r_9s_8i)

# supprime les positions trop proches sur un contig
data = Pedemo10_5r_9s_8i ; n = 1e4
manus = c(1)
for (i in 2:dim(data)[1]){   if (data$CHROM[i] == data$CHROM[i-1] & (data$POS[i]-data$POS[i-1]) < n ){} else {manus = c(manus,i)} }
Pedemo10_5r_9s_8i = Pedemo10_5r_9s_8i[manus,]
dim(Pedemo10_5r_9s_8i)


v= list() ; He = c()
homo = c("0/0","1/1","2/2","3/3")
hetero = c("0/1","0/2","0/3","1/2","1/3","2/3")

for(i in 1 : dim(Pedemo10_5r_9s_8i[,-c(1:9)])[2]){
  table =Pedemo10_5r_9s_8i[,-c(1:9)]
  table = table[,i]
  #levels(as.factor(substr(as.character(table),1,3)))
  v[[i]] = summary(as.factor(substr(as.character(table),1,3)))
  He[i] = sum(v[[i]][which(names(v[[i]]) %in% hetero)])/sum(v[[i]][which(names(v[[i]]) != ".")])
}

names(He) = colnames(Pedemo10_5r_9s_8i[,-c(1:9)])
He
plot(He, ylim = c(0.3,0.5))
text(x = c(1:length(He)), y = He, names(He), pos = 3)
summary(He)



#save2vcf(Eryth10_t)
tablobj2vcf(Pedemo10_5r_9s_8i,"data_vcf/Pedemo10_5r_9s_8i_pos1e4.csv","data_vcf/freebayes_vcf.head","data_vcf/Pedemo10_5r_9s_8i_pos1e4.vcf")

pop = c("apennina", "apennina"#, "apennina"
        ,"cottia","cottia","cottia"
        ,"pedemontana","pedemontana"
        ,"valgau","valgau"
        ,"hirsuta","hirsuta","hirsuta","hirsuta","hirsuta","hirsuta"
)

# analyse pour genind ####
Pedemo10_5r_9s_8i_v = read.vcfR("data_vcf/Pedemo10_5r_9s_8i_pos1e4.vcf", checkFile = T) ; Pedemo10_5r_9s_8i_v
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

fstat(Pedemo10_5r_9s_8i_g) ; beep(3)
fstat(Pedemo10_5r_9s_8i_g[which(Pedemo10_5r_9s_8i_g$pop %in% c("pedemontana","hirsuta"))]) ; beep(3)
#     pop             Ind
#Total 0.094443865327 -0.173026489424
#pop   0.000000000000 -0.295365847031

#que pour pede et valgau
#       pop             Ind
#Total 0.035734895138 -0.331346625899
#pop   0.000000000000 -0.380685269213

# que pour pede et cottia
#      pop             Ind
#Total 0.0306791851608 -0.282959844587
#pop   0.0000000000000 -0.323565763724

# que pour pede et hirsuta
#    pop             Ind
#Total 0.108398611131 -0.139744158652
#pop   0.000000000000 -0.278311331590

#This table provides the three F statistics F st (pop/total), F it (Ind/total), and F is
#(ind/pop). These are overall measures which take into account all genotypes and all loci.

# excés d'He dans les espèces, donc on peut supposer écart a HW

matFst <- pairwise.fst(Pedemo10_5r_9s_8i_g)
matFst
pop(Pedemo10_5r_9s_8i_g)
plot(nj(matFst))
#                 apennina         cottia       pedemontana        valgau
#  cottia       0.1003278856427
# pedemontana   0.1192152777891 0.0897982419468
# valgau        0.1358507278749 0.1010421724427 0.1033583099689
# hirsuta       0.1073449595393 0.1079755465314 0.0923285931921 0.0853606023649

x = Fst(Pedemo10_5r_9s_8i_g)


par(mfrow =c(1,1))
X <- tab(Pedemo10_5r_9s_8i_g, freq = TRUE, NA.method = "mean")
pca1 <- dudi.pca(X, scale = FALSE, scannf = FALSE, nf = 3)
barplot(pca1$eig[1:20], main = "PCA eigenvalues", col = heat.colors(50))

pop = c("apennina", "apennina"#, "apennina"
        ,"cottia","cottia","cottia"
        ,"pedemontana","pedemontana"
        #,"pedemontana","pedemontana"
       ,"valgau","valgau"
        ,"hirsuta","hirsuta","hirsuta","hirsuta","hirsuta","hirsuta"
)

pop(Pedemo10_5r_9s_8i_g) <- as.factor(pop)
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


compoplot(dapc1,lab="", ncol=2)

# vcfR ###
library(reshape2)
library(ggplot2)

myDiff <- genetic_diff(Pedemo10_5r_9s_8i_v, pops = as.factor(pop), method = 'nei')
#knitr::kable(head(myDiff[,1:15]))
dpf <- melt(myDiff[,c(8,14:17)], varnames=c('Index', 'Sample'), value.name = 'Depth', na.rm=TRUE)
colnames(myDiff)

ggplot(dpf, aes(x=variable, y=Depth)) + geom_violin(fill="#2ca25f", adjust = 1.2) + xlab("") + ylab("") + theme_bw()

# DAPC ####

set.seed(999)
pramx <- xvalDapc(tab(Pedemo10_5r_9s_8i_g, NA.method = "mean"), pop(Pedemo10_5r_9s_8i_g))
# mais c'est que c'est très moche dit donc...
set.seed(999)
system.time(pramx <- xvalDapc(tab(Pedemo10_5r_9s_8i_g, NA.method = "mean"), pop(Pedemo10_5r_9s_8i_g),
                              n.pca = 1:6, n.rep = 1000,
                              parallel = "multicore", ncpus = 7))
# utilisateur     système      écoulé
# 1221.338      17.356     196.886

names(pramx)
pramx[-1]

scatter(pramx$DAPC, cex = 2, legend = TRUE,
        clabel = FALSE, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1)

# PCAdapt analysis ####
#install.packages("pcadapt")

filename <- read.pcadapt("data_vcf/Pedemo10_5r_9s_8i_pos1e4.vcf", type = "vcf")
x <- pcadapt(input = filename, K = 5)
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
source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")

spider("data_vcf/Pedemo10_5r_9s_8i_pos1e4.vcf","VCF","data_vcf/Pedemo10_5r_9s_8i_pos1e4.str","STRUCTURE")
system("sed -i '1d' data_vcf/Pedemo10_5r_9s_8i_pos1e4.str ") #ne marche pas car premiere ligne en trop

struct2geno(file = "data_vcf/Pedemo10_5r_9s_8i_pos1e4.str", TESS = FALSE, diploid = T, FORMAT = 2,extra.row = 0, extra.col = 2, output = "data_vcf/Pedemo10_5r_9s_8i.geno")

obj  <- snmf("data_vcf/Pedemo10_5r_9s_8i_pos1e4.geno", K = 1:8, entropy = T, ploidy = 2,
             CPU = 7,repetitions = 10, project= "new", alpha=100)

# Choix du K optimal (20 runs)
ID = row.names(Pedemo10_5r_9s_8i_g$tab)

par(mfrow = c(1,1))
plot(obj, col = "blue", pch=1,cex=0.5)
beep(3)
color = c("orange","violet","lightgreen","red","blue","green","grey","black")

obj.snmf = snmf("data_vcf/Pedemo10_5r_9s_8i_pos1e4.geno", K = K, alpha = 100, project = "new", repetition = 10)
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
for (i in 1:11) Pop(i,"data_vcf/Pedemo10_5r_9s_8i_pos1e4.geno",  ID ) ;beep(3)



# introgress ####
mincov10_Eryth_CSV <- readr::read_delim("data_vcf/freebayes_-F0_3-n10-m13_-q20_mincov10_Eryth_SNPs_onlyCSV.csv","\t", escape_double = FALSE, trim_ws = TRUE)

c = c('CS1','CP1','CP4') #cottia
h = c('DMB','HC1','HGL','HS2','HP1','HPB') #hirsuta
p = c('PT1','PV1','GA2','GA4') #pedemontana

PedeHirsu10 = subset_reorder(mincov10_Eryth_CSV, c(c,p,h)) ; colnames(PedeHirsu10)
PedeHirsu10_5r = rare(PedeHirsu10, rare  = 0.05, r= T) ; colnames(PedeHirsu10_5r)
PedeHirsu10_5r = PedeHirsu10_5r[which(PedeHirsu10_5r$QUAL >= 20),]
PedeHirsu10_5r_9s_8i = tri(data = PedeHirsu10_5r,n.r=1,n.c = 0.8, r = T) # avec ces seuils on vire l'individu AML car trop d'info manquantes pour cet individu
# supprime les positions trop proches sur un contig
data = PedeHirsu10_5r_9s_8i ; n = 1e4
manus = c(1)
for (i in 2:dim(data)[1]){   if (data$CHROM[i] == data$CHROM[i-1] & (data$POS[i]-data$POS[i-1]) < n ){} else {manus = c(manus,i)} }
PedeHirsu10_5r_9s_8i = PedeHirsu10_5r_9s_8i[manus,]
dim(PedeHirsu10_5r_9s_8i )
tablobj2vcf(PedeHirsu10_5r_9s_8i,"data_vcf/PedeHirsu10_5r_9s_8i.csv","data_vcf/freebayes_vcf.head","data_vcf/PedeHirsu10_5r_9s_8i.vcf")


levels = apply(PedeHirsu10_5r_9s_8i[,-c(1:9)],1, function(data) levels(as.factor(substr(as.character(data),1,3)))[!levels(as.factor(substr(as.character(data),1,3))) %in% "."])
summary(as.factor(unlist(levels)))

# python2 convert_vcf2introgress_genotypes.py data_vcf/PedeHirsu10_5r_9s_8i.vcf  > data_vcf/PedeHirsu10_5r_9s_8i.introgress
# python2 convert_vcf2introgress_loci.py data_vcf/PedeHirsu10_5r_9s_8i.vcf  > data_vcf/PedeHirsu10_5r_9s_8i.loci

vcf2csv("data_vcf/PedeHirsu10_5r_9s_8i.vcf","data_vcf/PedeHirsu10_5r_9s_8i.head","data_vcf/PedeHirsu10_5r_9s_8i.csv")
PedeHirsu10_5r_9s_8i <- read_delim("data_vcf/PedeHirsu10_5r_9s_8i.csv", "\t", escape_double = FALSE, trim_ws = TRUE)

LociData = PedeHirsu10_5r_9s_8i[,c(1,2)]
LociData$locus = paste("c",substr(LociData$CHROM,7,10),".",LociData$POS, sep = "")
LociData$type = "C"
LociData$lg = substr(LociData$CHROM,7,10)
LociData = as.data.frame(LociData[,-c(1,2)])

LociData[1,]

AdmixData = PedeHirsu10_5r_9s_8i[,-c(1:9)]
AdmixData = sapply(AdmixData, substring, 1, 3)
NAs = AdmixData == "." ; AdmixData[NAs == T] = "NA/NA" ; rm(NAs)

Pop = c("Pop1","Pop1","Pop1","Pop1","Pop1","Pop2","Pop2","Pop3","Pop3","Pop3","Pop3","Pop3","Pop3")
Ind = colnames(AdmixData)
AdmixData = rbind(Pop,Ind,AdmixData)
colnames(AdmixData) = NULL

Parent_1 = AdmixData[-c(1,2),c(1:5)] ; rownames(Parent_1) = LociData[,1]
Parent_2 = AdmixData[-c(1,2),-c(1:7)] ; rownames(Parent_2) = LociData[,1]

count.matrix = prepare.data(admix.gen = AdmixData,
                            loci.data = LociData, fixed = F,
                            parental1 = Parent_1 , parental2 = Parent_2,
                            pop.id = T, ind.id = T)
hi.index.sim = est.h(introgress.data = count.matrix,
                     loci.data = LociData)
mk.image(introgress.data = count.matrix, loci.data = LociData,
         marker.order = NULL,hi.index = hi.index.sim,ylab.image = "Individuals",
         xlab.h = "population 2 ancestry", pdf = F)
abline(h= c(5,6,7), col = "red", lty = 3)
text(x = hi.index.sim[c(5,6,7),2],y = c(5,6,7), Ind[5:7], pos = c(4,4,4))
rownames(hi.index.sim) = Ind
hi.index.sim

# meh
gen.out = genomic.clines(introgress.data = count.matrix,hi.index = hi.index.sim, loci.data = LociData, sig.test = T, method = "permutation")

#gen.out$Summary.data
beep(3)

#clines.plot(cline.data = gen.out, rplots = 3, cplots = 3, pdf = F)




