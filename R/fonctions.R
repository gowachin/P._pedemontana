#' clean
#'
#' delete the monomorphic loci of a vcf file imported as a dataframe
#'
#' delete row when there is only one genotype exept the missing data. Use an apply function by row, on the data.frame[,10:n.col], because it expect the file to be a vcf originally, with the first 9 column to be informativ about the loci.
#' It only use the first 3 letters from each cell, expecting "0/0" format, and compute for levels(as.factor("the_row_itself")). In "0/0", each 0 is an allele, with the first allele before the slash and the second after.
#'
#' @param data a data frame from a vcf file (without the hashtag lignes)
#'
#' @return same data.frame as beginning, whithout the monomorphic rows
#' @return print some data about the percentage of data deleted by this function
#'
#' @author JAUNATRE Maxime
#'
#' @export
clean = function(data) {
  n.col= dim(data)[2]
  compte = apply(data[,10:n.col],1,function(data) length(levels(as.factor(substr(as.character(data),1,3)))[!levels(as.factor(substr(as.character(data),1,3))) %in% "."] ) >1)
  #exterminate = apply(data[,10:n.col],1,paste,collapse ="")
  #compte = exterminate == paste(rep(".",23),collapse = "")
  print(summary(compte))
  resum = c(pourc.clean = (1-sum(compte)/length(compte) )*100)
  print(resum)

  end = data[which(compte == T),]
  #end = data[which(compte == F),]
  return(end)
}

#' subset_reorder
#'
#' subset the main information with the list of individual and clean the dataset for monomorphic loci
#'
#' bind columns in the same order that the list put them, with a for() loop. Begin with the same 9 column because the file is expected to be a vcf.
#'
#' @param data a data frame from a vcf file (without the hashtag lignes)
#' @param list vector of individual to maintain in this subset (individual ID must be the same as the colnames in the dataframe)
#'
#' @return same data.frame, but columns are not in the same order, after the 9th column.
#'
#' @author JAUNATRE Maxime
#'
#' @export
subset_reorder = function(data,list) {
  manus = data[,1:9]
  for (i in 1:length(list)) {
    manus = cbind(manus,data[,which(colnames(data) == list[i])])
  }
  manus = clean(manus)
  return(manus)
}

#' rare
#'
#' delete the alleles with less than a "rare" percentage of presence in the dataset
#' delete the other variants, apart the two main allele (O and 1)
#'
#' @param data a data frame from a vcf file (without the hashtag lignes)
#' @param rare percentage of presence needed for keeping the allele in the set (expressed as "0.05" for 5\%)
#' @param r if the function return something : the data or the resume (dimension of the final dataframe, \% of deletion)
#' @param p if the function return the resume only
#'
#'
#' @author JAUNATRE Maxime
#'
#' @export
rare = function(data,rare = 0,quiet = T, r= F,p = F) { # data est un data frame type inform_mincov10 , quiet est le rendu
  n.col = dim(data)[2] ; n.row = dim(data)[1]
  matrix = as.matrix(data[,10:n.col])

  H. = apply(matrix,1, function(matrix) sum(substr(as.character(matrix),1,1) != ".")*2 )
  H0 = apply(matrix,1, function(matrix) sum(substr(as.character(matrix),1,1) == "0")+sum(substr(as.character(matrix),3,3) == "0") )/H. >= rare
  H1 = apply(matrix,1, function(matrix) sum(substr(as.character(matrix),1,1) == "1")+sum(substr(as.character(matrix),3,3) == "1") )/H. >= rare
  H2 = apply(matrix,1, function(matrix) sum(substr(as.character(matrix),1,1) == "2")+sum(substr(as.character(matrix),3,3) == "2") )/H. >= rare
  H3 = apply(matrix,1, function(matrix) sum(substr(as.character(matrix),1,1) == "3")+sum(substr(as.character(matrix),3,3) == "3") )/H. >= rare

  Haplo1 = matrix(substr(as.character(matrix),1,1), nrow = dim(matrix)[1], ncol = dim(matrix)[2])
  Haplo2 = matrix(substr(as.character(matrix),3,3), nrow = dim(matrix)[1], ncol = dim(matrix)[2])

  matrix[Haplo1 == "0"& H0 == F | Haplo2 == "0" & H0 == F] = "."
  matrix[Haplo1 == "1"& H1 == F | Haplo2 == "1" & H1 == F] = "."
  # matrix[Haplo1 == "2"& H2 == F | Haplo2 == "2" & H2 == F] = "."
  # matrix[Haplo1 == "3"& H3 == F | Haplo2 == "3" & H3 == F] = "."
  #tri considerant allele > 1 comme erreurs
  matrix[Haplo1 == "2" | Haplo2 == "2"] = "."
  matrix[Haplo1 == "3" | Haplo2 == "3"] = "."

  data = cbind(data[,1:9],matrix)
  data = clean(data)
  final.row = dim(data)[1]

  resum = c(n.site  =  final.row ,
            discar.site = n.row - final.row,
            p.discard = (n.row - final.row)/n.row )

  if(quiet == F) print(resum)
  if (r== T & p == F) return(data)
  if (r== T & p == T) return(resum)
}

#' tri
#'
#' delete in a first time the loci with a certain percentage of missing data, and in a second time individuals with too many missing data
#'
#' @param data a data frame from a vcf file (without the hashtag lignes)
#' @param n.r percentage of row containing information (1-missing data)
#' @param c.r idem but for individuals
#' @param quiet if the function print information at the end about percentage of data discarded
#' @param r if the function return the data or the resume
#' @param p if the function return the resume only
#'
#' @author JAUNATRE Maxime
#'
#' @export
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

  data = clean(data)

  resum = c(r.pourc.end = sum(row)/length(row)*100 ,
            r.n.site = sum(row),
            c.pourc.end = sum(col)/length(col)*100 ,
            c.n.ind = sum(col))

  if(quiet == F) print(resum)
  if (r== T & p == F) return(data)
  if (r== T & p == T) return(resum)
}

#' vcf2csv
#'
#' save as a csv file a vcf file, allowing to apply treshold on the data
#'
#' @param vcf the original vcf file
#' @param head the txt file where to stock the head of the vcf
#' @param csv the name of the final csv
#'
#' @author JAUNATRE Maxime
#'
#' @export
vcf2csv = function(vcf, head, csv) {
system(paste("./csv_makup.sh",vcf,head,csv, sep = " "))
}


#' tablobj2vcf
#'
#' save a dataframe as a vcf file, concatening the head of a vcf
#'
#' @param obj name of the data frame to save
#' @param name name for save
#' @param head the txt file where is stock the head of the vcf
#' @param vcf name of the final file
#'
#'
#' @author JAUNATRE Maxime
#'
#' @export
tablobj2vcf = function(obj,name,head,vcf) {
  write.table(obj, name ,sep = "\t", quote = F, row.names=F)

  system(paste(" cp ",name," ",name,".temp",sep = ""))
  system(paste(" sed -i '1s/.*/#&/' ",name,".temp",sep = ""))
  system(paste("cat ",head," ",name,".temp"," > ",vcf, sep = ""))
  system(paste("rm ",name,".temp",sep = ""))
}

#' Pop
#'
#' str barplot
#'
#' @param k the K to test
#' @param the file use in cross entropy (.geno)
#' @param ID  individual tags
#'
#' barplot the structure of the geno file
#'
#' @author JAUNATRE Maxime
#'
#' @export
Pop = function(K, files,ID) {
  obj.snmf = snmf(files, K = K, alpha = 10, project = "new",iterations = 2000, repetitions = 20,
                  CPU = 7, entropy = T)
  ce=cross.entropy(obj.snmf,K=K)
  best = which.min(ce)
  qmatrix = Q(obj.snmf, K = K, run = best)
  barplot(t(qmatrix), col = color, border = 1, space = 0.05,
          xlab = "Individus", ylab = "Coefficient d'admixture", main = paste("K= ",K,sep = ""),
          names.arg =ID, las = 2)}

#' spider
#'
#' push file from a type to another using pgdspider.
#' For now, just working with VCF to structure, because it need a .spid done for the software
#'
#' @param input the name of the entry file
#' @param inFORM the type of entry file
#' @param output the name of the exit file
#' @param outFORM the type of exit file
#'
#' @author JAUNATRE Maxime
#'
#' @export
spider = function(input,inFORM,output,outFORM) {
  if (inFORM == "VCF" & outFORM == "STRUCTURE") {spid = "Spid_VCF_STRUCTURE.spid"}
  if (inFORM == "VCF" & outFORM == "PED") {spid = "Spid_VCF_PED.spid"}

  command = paste("cd ; java -Xmx1024m -Xms512M -jar Téléchargements/PGDSpider_2.1.1.5/PGDSpider2-cli.jar -inputfile Bureau/BEE/Stage/Pedemontana/",
                  input, " -inputformat ", inFORM ," -outputfile Bureau/BEE/Stage/Pedemontana/", output, " -outputformat ", outFORM ,
                  " -spid Bureau/BEE/Stage/Pedemontana/data_vcf/",spid,sep = "")
  print(command)
  system(command)
}

#' files
#'
#' create files for the population genetic analysis
#' .vcf
#' .str (structure)
#' .geno (genotype for LEA)
#' .snp (for DIYABC)
#'
#' @param obj the data frame originaly made from a vcf file
#' @param name the name of the file to output
#' @param ind the list of individuals
#' @param pop the list of population assignations
#' @param head a txt file with the vcf head to take from
#'
#'
#' @author JAUNATRE Maxime
#'
#' @export
files = function(obj,name,ind,pop, head= NULL) {
.csv = paste(name,".csv",sep="")
.vcf = paste(name,".vcf",sep="")
.str = paste(name,".str",sep="")
.geno = paste(name,".geno",sep="")
.snp = paste(name,".snp",sep="")

if( is.null(head) ) {head="data_vcf/freebayes_vcf.head"} else {}

cat("\n");cat("writing VCF \n")
tablobj2vcf(obj, .csv , head , .vcf)
cat("\n");cat("VCF done, writing STR \n")
spider( .vcf ,"VCF", .str ,"STRUCTURE")
system(paste("sed -i '1d' ", .str , sep=""))

source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
suppressWarnings(source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R"))

cat("\n");cat("STR done, writing GENO \n")
struct2geno(file = .str , TESS = FALSE, diploid = T, FORMAT = 2,extra.row = 0, extra.col = 2, output =  .geno)

x =c("barplotCoeff", "barplotFromPops", "correlation",
     "correlationFromPops", "createGrid", "createGridFromAsciiRaster",
     "defaultPalette", "displayLegend", "fst",
     "getConstraintsFromAsciiRaster", "helpPops", "lColorGradients",
     "loadPkg", "maps", "mapsFromPops",
     "mapsMethodMax", "struct2geno")
rm(list = x, envir = .GlobalEnv)

cat("\n");cat("GENO done, writing SNP \n")
temp.geno = read.geno( .geno )
temp.geno= rbind(rep("A",dim(temp.geno)[2]),temp.geno)
temp.geno = cbind(c("IND",as.character(ind)),c("SEX",rep("9",dim(temp.geno)[1]-1)),c("POP",as.character(pop)),temp.geno)
colnames(temp.geno) = NULL
write.table(temp.geno, .snp ,sep = "\t", quote = F, row.names=F, col.names = F)
system(paste("sed -i '1 i\ <NM=1NF>' ", .snp ,sep=""))

fichiers = list(  .csv = paste(name,".csv",sep=""),
                  .vcf = paste(name,".vcf",sep=""),
                  .str = paste(name,".str",sep=""),
                  .geno = paste(name,".geno",sep=""),
                  .snp = paste(name,".snp",sep=""),
                  .ind = ind,
                  .pop = pop)
return(fichiers)
}

#' subset_ord_pop
#'
#' reorder and subset the population assignement
#'
#' @param population a dataframe with two column ind and pop, row is an individual and its population assignement
#' @param sub a vector with the individuals to subset
#'
#'
#' @return a data.frame
#'
#' @author JAUNATRE Maxime
#'
#' @export
subset_ord_pop = function(population,sub) {
  manus = data.frame()
  data = population[which(population$ind %in% sub),]
  for(i in 1:length(sub)){
    manus = rbind(manus,data[which(data$ind == sub[i]),])
  }
  return(as.data.frame(manus))
}

#' dataset
#'
#' create all thes files for population genetic
#'
#' @param ind vector of individuals
#' @param popfile csv file with population assignement, first column with individual ID, second column with population ID
#' @param entryfile a vcf file
#' @param name basic name and location of the files to (character chain)
#' @param rare numeric for the percentage of presence minimal for an allele to be retain (ex 0.05)
#' @param qual numeric for the quality of a loci to be retain (ex 20)
#' @param missLoci numeric for the percentage of missing data of an allele to be retain (ex 0.8)
#' @param missInd numeric for the percentage of missing data of an individual to be retain (ex 0.8)
#' @param LD numeric, minimal distance between two loci on a contig (ex 1e4)
#'
#'
#' @return a list with all the files created for this analysis
#'
#' @author JAUNATRE Maxime
#'
#' @export
dataset = function(ind,popfile,entryfile,name,rare=0,qual=0,missLoci=0,missInd=0,LD=1) {

  head = paste(name,"_head.txt",sep="")
  csv = paste(name,"_CSV.csv",sep="")

  vcf2csv(vcf = entryfile,head = head, csv = csv)

  entryfile = csv

cat("\n");cat("reading file \n")
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
  cat("\n");cat("deleting loci with narrow positions \n")
  manus = c(1)
  pb <- txtProgressBar(min = 2, max = dim(missingdata)[1], style = 3)
  for (i in 2:dim(missingdata)[1]){   if (missingdata$CHROM[i] == missingdata$CHROM[i-1] & (missingdata$POS[i]-missingdata$POS[i-1]) < LD ){} else {manus = c(manus,i)}
    setTxtProgressBar(pb, i)}
  close(pb)
  final = missingdata[manus,]} else {
    cat("\nno LD filter \n")
    final = missingdata}

cat("\n");cat("reading pop file \n")
Populations <- as.data.frame(read.csv(popfile))
colnames(Populations) = c("ind","pop")
final_pop = subset_ord_pop(population =Populations,sub=c(colnames(final[,-c(1:9)])))

resum = list(individuals = colnames(final[,-c(1:9)]), file_dim = dim(final))

cat("\n")
print(resum)
cat("\n");cat("saving files \n")
name= paste(name,"_r",rare,"_q",qual,"_mL",missLoci,"_mI",missInd,"_LD1e",log10(LD),sep="")
fichier = files(obj =final, name=name, ind=final_pop$ind, pop=final_pop$pop)

system(paste("rm",csv,head,sep = " "))

fichier$.ind = as.factor(as.character(fichier$.ind))

fichier$.pop = factor(as.character(fichier$.pop),unique(final_pop$pop))

cat("\n");cat("DONE :) \n\n")

return(fichier)
}


#' ABBA-BABA test
#'
#' compute the durand's D for admixture test between different populations
#'
#' @param file a list containning an element named .csv (a character chain with the location of a .csv file) and a vector of factor .pop containing the assignation to populations for each individual. usually made by the dataset function
#' @param P1,P2,P3,Root a character chain for the populations to compare in this test
#' @param n number of bootstrap.
#'
#'
#' @return a list with a vector (all the D value, first one is without bootstrap) and a value for the number of informativ loci
#'
#' @author JAUNATRE Maxime
#'
#' @export
durand_freq_D = function(file,P1,P2,P3,Root,n) {
  
  table = suppressMessages(substr(as.matrix(readr::read_delim(file$.csv,"\t", escape_double = FALSE, trim_ws = TRUE)[,-c(1:9)]),1,3))
  table = rbind(as.character(file$.pop),table)
  #verif
  if (P1 %in% file$.pop == F ) {message("P1 is not in file$.pop")}
  if (P2 %in% file$.pop == F ) {message("P2 is not in file$.pop")}
  if (P3 %in% file$.pop == F ) {message("P3 is not in file$.pop")}
  if (Root %in% file$.pop == F ) {message("Root is not in file$.pop")}
  
  #if ( length(unique.default(c(P1,P2,P3,Root)))<4 ) {
  #  forward = menu(c("yes","no"), title = "2 populations are the same, wish to continue?")
  #  if (forward == 1) {} else {stop("Assign other populations")}
  #  }
  
  echanti = c(P1,P2,P3,Root)
  manus = matrix(ncol = 4, nrow = dim(table)[1])
  for(i in 1:4) {
    matrix = table[,which(table[1,]==echanti[i])]
    H. = apply(matrix,1, function(matrix) sum(substr(as.character(matrix),1,1) != ".")*2 )
    H1 = apply(matrix,1, function(matrix) sum(substr(as.character(matrix),1,1) == "1")+sum(substr(as.character(matrix),3,3) == "1") )/H.
    manus[,i] = H1
  }
  
  x = apply(manus,1,sum) == 0
  manus = manus[which(x == F),]
  
  abba = (1-manus[,1])*manus[,2]*manus[,3]*(1-manus[,4])
  baba = manus[,1]*(1-manus[,2])*manus[,3]*(1-manus[,4])
  
  a_b = abba != 0 ; b_a = baba != 0
  inf = sum(a_b == T & b_a == T)
  #print(c("sites inform" = inf))
  D= c(sum(abba-baba)/sum(abba+baba))
  
  #pb <- txtProgressBar(min = 1, max = n, style = 3)
  for (i in 1: n+1){
    tmp.manus = manus[sample(c(1:dim(manus)[1]),dim(manus)[1], replace = T ),] ; tmp.manus
    
    abba = (1-tmp.manus[,1])*tmp.manus[,2]*tmp.manus[,3]*(1-tmp.manus[,4])
    baba = tmp.manus[,1]*(1-tmp.manus[,2])*tmp.manus[,3]*(1-tmp.manus[,4])
    
    D[i] = sum(abba-baba)/sum(abba+baba)
    
    # setTxtProgressBar(pb, i)
  }
  #close(pb)
  
  return(list(D,inf))
}
