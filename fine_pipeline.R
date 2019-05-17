out = c("AP1")
a = c('AMB','AOL') #apenina  j'ai enlevé AML
c = c('CS1','CP1','CP4') #cottia
d = c('DGB','DRL') #daonensis
h = c('DMB','HC1','HGL','HS2','HP1','HPB') #hirsuta
p = c('PT1','PV1') #pedemontana
va = c('GA2','GA4')
v = c('VR3','VR1','VL2','VB1') #villosa



Eryth.file = dataset.g(ind= c(a,p,va,c)
                     ,popfile= "Populations.csv"
                     ,entryfile= "data_vcf/freebayes_-F0.3-n10-m13_-q20_mincov10_90samples_SNPs_only.vcf"
                     ,name = "clean_space/testing"
                     ,rare= 0.05,qual= 20,missLoci= 0.95,missInd= 0,LD= 1e4)
beep(3)

Eryth.file$.ind

dataset.g = function(ind,popfile,entryfile,name,rare=0,qual=0,missLoci=0,missInd=0,LD=0) {

  head = paste(name,"_head.txt",sep="")
  csv = paste(name,"_CSV.csv",sep="")

  vcf2csv(vcf = entryfile,head = head, csv = csv)

  entryfile = csv

  cat("\nreading file \n")
  lecture <- readr::read_delim(entryfile,"\t", escape_double = FALSE, trim_ws = TRUE)

  if(identical(colnames(lecture[,-c(1:9)]) , ind) == F) {
  cat("\nreordering data \n")
  reorder = subset_reorder(lecture, ind )} else {
    cat("\nno reordering needed \n")
    reorder = lecture}

  if(rare > 0) {
    cat("\ndeleting rare allele \n")
    Rare = rare(reorder, rare  = rare, r= T)} else {
      cat("\nno rare filter \n")
      Rare = reorder}

  if(qual > 0) {
    cat("\ndeleting poor quality SNP \n")
    quality = Rare[which(Rare$QUAL >= qual),]} else {
      cat("\nno qual filter \n")
      quality = Rare}

  if(missLoci > 0 | missInd > 0 ) {
    cat("\napplying treshold to missing data \n")
    missingdata = tri(data = quality,n.r=missLoci,n.c = missInd, r = T)} else {
      cat("\nno missing data filter \n")
      missingdata = quality}

  if(LD > 1) {
    cat("\ndeleting loci with narrow positions \n")
    manus = c(1)
    pb <- txtProgressBar(min = 2, max = dim(missingdata)[1], style = 3)
    for (i in 2:dim(missingdata)[1]){
      if (missingdata[i,1] == missingdata[i-1,1] & (missingdata[i,2]-missingdata[i-1,2]) < LD ){

      } else {manus = c(manus,i)}
      setTxtProgressBar(pb, i)

      }
    close(pb)
    final = missingdata[manus,]} else {
      cat("\nno LD filter \n")
      final = missingdata}

  cat("\nreading pop file \n")
  Populations <- as.data.frame(read.csv(popfile))
  colnames(Populations) = c("ind","pop")
  final_pop = subset_ord_pop(population =Populations,sub=c(colnames(final[,-c(1:9)])))

  resum = list(individuals = colnames(final[,-c(1:9)]), file_dim = dim(final))

  cat("\n")
  print(resum)
  cat("\nsaving files \n")
  name= paste(name,"_r",rare,"_q",qual,"_mL",missLoci,"_mI",missInd,"_LD1e",log10(LD),sep="")
  fichier = files(obj =final, name=name, ind=final_pop$ind, pop=final_pop$pop, head = head)

  system(paste("rm",csv,head,sep = " "))

  fichier$.ind = as.factor(as.character(fichier$.ind))
  fichier$.pop = as.factor(as.character(fichier$.pop))

  cat("\n");cat("DONE :) \n\n")

  return(fichier)
}
