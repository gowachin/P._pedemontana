library(raster)

library(Pedemontana)
library(LEA)
library(vcfR)
library(readr)
library(ade4)
library(adegenet)
library(hierfstat)
library(ape)


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



PedeHirsu.file = dataset(ind= c(a,c,p,va,h)
                         ,popfile= "Populations.csv"
                         ,entryfile= "data_vcf/freebayes_-F0_3-n10-m13_-q20_mincov10_Eryth_SNPs.vcf"
                         ,name = "data_vcf/Hirsuta"
                         ,rare= 0.05,qual= 20,missLoci= 0.95,missInd= 0,LD= 1e4)
# adegenet ####

file = Eryth.file
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

GENIND$pop
rownames(GENIND$tab)

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
scatter(dapc1,scree.da=FALSE, bg="white", posi.pca="topright", legend=F,
        txt.leg=paste("group", 1:7), col=c("red","blue"))

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

myCol <- c("slateblue2","darkorchid3","violetred2","steelblue3","goldenrod2","darkolivegreen3","lightsalmon3")
scatter(dapc1, posi.da="topleft", bg="gray90",
        pch=17:25, cstar=0, col=myCol, scree.pca=F, clabel = 1.4)


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
obj.snmf = snmf(file$.geno, K = 3, alpha = 10, project = "new")
qmatrix = Q(obj.snmf, K = 3)
barplot(t(qmatrix), col = color, border = 1, space = 0.05, xlab = "Individuals", ylab = "Admixture coefficients",
        names.arg =file$.ind, las = 2)
par(mfrow = c(4,1))
for (i in 2:5) Pop(K= i,files = file$.geno, ID= file$.ind) ;beep(3)

legend("bottomright", c("P. hirsuta","P. cottia","P. apennina","P. pedemontana","P. pedemontana Ecrins"),
       pch = c(1,13,8), text.font=c(3), cex = 0.7)

# admixture


Data=data.frame()
P1=c("apennina","cottia","pedemontana")
P2=c("apennina","cottia","ecrins")
P3=c("hirsuta")
Root=c("daonensis")
i=1
pb <- txtProgressBar(min = 1, max = length(P1)*length(P2)*length(P3)*length(Root), style = 3)
for (k in P1){
  for (j in P2) {
    for (l in P3) {
      for (m in Root) {
        #cat("\n",c(k,j,l,m),"\n")
        Data[i,1]=k
        Data[i,2]=j
        Data[i,3]=l
        Data[i,4]=m
        D = durand_freq_D(file,P1 =k,P2=j,P3 =l,Root =m,n=1e4)
        d = D[[1]]
        inf = D[[2]]
        Data[i,5]= mean(d)
        Data[i,6] = sum(d<0)/length(d)
        Data[i,7] = sum(d>0)/length(d)
        
        z <- abs(d[1]/sd(d[-1]))
        new.pval <- 2 * (1 - pnorm(z))
        Data[i,8] = new.pval
        Data[i,9] = z
        Data[i,10] = inf
        
        setTxtProgressBar(pb, i)
        i=i+1
      }
    }
  }
}
close(pb)

colnames(Data) = c("P1","P2","P3","Root",
                   "d","p<0","p>0","pval","z","site inf")
View(Data)
beep(3)

Data = rbind(Data[7,],Data[8,],Data[9,],Data[3,],Data[6,])
Data = Data[,-c(6,7,9)]
write.table(Data,"Rendu/fig/ABBA.csv",sep=" & ")

length(summary(file$.pop))
P1 = "pedemontana" ; P2 = "ecrins" ; P3 = "hirsuta" ; Root = "daonensis"

x =durand_freq_D(file,P1,P2,P3,Root,n=1e4)
hist(x)
z <- abs(x[1]/sd(x[-1]))
2 * (1 - pnorm(z))
sum(x<0)/length(x)
sum(x>0)/length(x)

z <- abs(d[1]/sd(d[-1]))
new.pval <- 2 * (1 - pnorm(z))
new.pval
