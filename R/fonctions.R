#' not in
#'
#' @param x one object
#' @param y second object
#'
#' The negate of the %in% function :
#' match returns a vector of the positions of (first) matches of its first argument NOT in its second.
#'
#' @author JAUNATRE Maxime from (https://stackoverflow.com/questions/5831794/opposite-of-in)
#'
#' @export
'%.in%' = function(x,y) {!('%in%'(x,y))}


#' create documentation
#'
#' @author JAUNATRE Maxime
#'
#' @export
manual = function() {
  library("devtools")
  document()
}

#' plot a raster
#'
#' @param rast the raster file
#' @param ext the extent of this raster
#' @param T.ext T if you want to crop with an extent
#' @param line T if you want the water limit
#'
#' @author JAUNATRE Maxime
#'
#' @export
plot.raster = function(rast,ext=extent (0,0,0,0),line = F,raw =T) {
  if(raw == T) {rast=raster(rast)}
  if (ext != extent (0,0,0,0)) {plot(crop(rast, ext)) } else {plot(rast) }
  if(line==T){lines(read.table("data_carto/WORLD_lowres.dat"))}
}


#' clean
#'
#' @param data the csv from a vcf file
#'
#' delete the monomorphic loci of a vcf file imported as a CSV#'
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

#' rare
#'
#' @param data the csv from a vcf file
#' @param rare percentage of presence needed for keeping the allele in the set
#' @param r if the function return the data or the resume
#' @param p if the function return the resume only
#'
#' delete the alleles with less than a "rare" percentage of presence in the dataset
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

#' clean
#'
#' @param data the csv from a vcf file
#' @param n.r percentage of row containing information (1-missing data)
#' @param c.r idem but for individuals
#' @param quiet if the function print information at the end about percentage of data discarded
#' @param r if the function return the data or the resume
#' @param p if the function return the resume only
#'
#' delete in a first time the loci with a certain percentage of missing data, and in a second time individuals with too many missing data
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

  resum = c(r.pourc.end = sum(row)/length(row)*100 ,
            r.n.site = sum(row),
            c.pourc.end = sum(col)/length(col)*100 ,
            c.n.ind = sum(col))

  if(quiet == F) print(resum)
  if (r== T & p == F) return(data)
  if (r== T & p == T) return(resum)
}


#' subset and reorder
#'
#' @param data the csv from a vcf file
#' @param list list of individual to maintain in this subset
#'
#' subset the main information with the list of individual and clean the dataset for monomorphic loci
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


#' save to vcf a csv
#'
#' @param csv the dataframe with the information to save inside a vcf file
#'
#' save as a vcf file but need to be a vcf troncated at the origine
#' hide because tablobj2vcf is way better
#'
#' @author JAUNATRE Maxime
#'
save2vcf = function(csv) {
  name = deparse(substitute(csv))
  write.table(csv, paste("data_vcf/",name,".csv",sep = ""),sep = "\t", quote = F, row.names=F)
  system(paste(" ./CSV_to_VCF.sh data_vcf/",name,".csv ; mv data_vcf/",name,".csv data_vcf/",name,".vcf",sep="") )
}



#' save to csv a vcf
#'
#' @param vcf the original vcf file
#' @param head the txt file where to stock the head of the vcf
#' @param csv the name of the final csv
#'
#' save as a csv file a vcf file, allowing to apply treshold on the data
#'
#' @author JAUNATRE Maxime
#'
#' @export
vcf2csv = function(vcf, head, csv) {
system(paste("./csv_makup.sh",vcf,head,csv, sep = " "))
}

#' save to vcf a csv
#'
#' @param obj name of the data frame to save
#' @param name name for save
#' @param head the txt file where is stock the head of the vcf
#' @param vcf name of the final file
#'
#' save a dataframe as a vcf file, concatening the head of a vcf
#'
#' @author JAUNATRE Maxime
#'
#' @export
tablobj2vcf = function(obj,name,head,vcf) {
  write.table(obj, name ,sep = "\t", quote = F, row.names=F)
  system(paste(" sed -i '1s/.*/#&/' ",name,sep = ""))
  system(paste("cat ",head," ",name," > ",vcf, sep = ""))
  system(paste("rm ",name,sep = ""))
}

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
Pop = function(K, file,ID) {
  obj.snmf = snmf(file, K = K, alpha = 100, project = "new",iterations = 2000, repetitions = 20,
                  CPU = 7)
  ce=cross.entropy(obj,K=K)
  best = which.min(ce)
  qmatrix = Q(obj.snmf, K = K, run = best)
  barplot(t(qmatrix), col = color, border = NA, space = 0,xlab = "Individuals", ylab = "Admixture coefficients",
          names.arg =ID, las = 2)}

#' PGDSpider from r
#'
#' @param input the name of the entry file
#' @param inFORM the type of entry file
#' @param output the name of the exit file
#' @param outFORM the type of exit file
#'
#' push file from a type to another using pgdspider.
#' For now, just working with VCF to structure, because it need a .spid done for the software
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


#' files creation
#'
#' @param obj the data frame originaly made from a vcf file
#' @param name the name of the file to output
#' @param ind the list of individuals
#' @param pop the list of population assignations
#' @param head a txt file with the vcf head to take from
#'
#' create files for the population genetic analysis
#' .vcf
#' .str (structure)
#' .geno (genotype for LEA)
#' .snp (for DIYABC)
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

cat("\n");cat("writing VCF")
tablobj2vcf(obj, .csv , head , .vcf)
cat("\n");cat("VCF done, writing STR")
spider( .vcf ,"VCF", .str ,"STRUCTURE")
system(paste("sed -i '1d' ", .str , sep=""))

source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
suppressWarnings(source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R"))

cat("\n");cat("STR done, writing GENO")
struct2geno(file = .str , TESS = FALSE, diploid = T, FORMAT = 2,extra.row = 0, extra.col = 2, output =  .geno)

x =c("barplotCoeff", "barplotFromPops", "correlation",
     "correlationFromPops", "createGrid", "createGridFromAsciiRaster",
     "defaultPalette", "displayLegend", "fst",
     "getConstraintsFromAsciiRaster", "helpPops", "lColorGradients",
     "loadPkg", "maps", "mapsFromPops",
     "mapsMethodMax", "struct2geno")
rm(list = x, envir = .GlobalEnv)

cat("\n");cat("GENO done, writing SNP")
temp.geno = read.geno( .geno )
temp.geno= rbind(rep("A",dim(temp.geno)[2]),temp.geno)
temp.geno = cbind(c("IND",as.character(ind)),c("SEX",rep("9",dim(temp.geno)[1]-1)),c("POP",as.character(pop)),temp.geno)
colnames(temp.geno) = NULL
write.table(temp.geno, .snp ,sep = "\t", quote = F, row.names=F, col.names = F)
system(paste("sed -i '1 i\ <NM=1NF>' ", .snp ,sep=""))

fichiers = list(  .vcf = paste(name,".vcf",sep=""),
                  .str = paste(name,".str",sep=""),
                  .geno = paste(name,".geno",sep=""),
                  .snp = paste(name,".snp",sep="") )
return(fichiers)
}

#' reorder and subset the population assignement
#'
#' @param population a dataframe with two column ind and pop, row is an individual and its population assignement
#' @param sub a vector with the individuals to subset
#'
#' reorder and subset the population assignement
#'
#' return a data.frame
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


#' create de dataset files
#'
#' @param ind vector of individuals
#' @param popfile csv file with population assignement
#' @param entryfile csv file originally a vcf
#' @param name basic name and location of the files to (character chain)
#' @param rare numeric for the percentage of presence minimal for an allele to be retain (ex 0.05)
#' @param qual numeric for the quality of a loci to be retain (ex 20)
#' @param missLoci numeric for the percentage of missing data of an allele to be retain (ex 0.8)
#' @param missInd numeric for the percentage of missing data of an individual to be retain (ex 0.8)
#' @param LD numeric, minimal distance between two loci on a contig (ex 1e4)
#'
#' create all thes files for population genetic
#'
#' return a list with all the files created for this analysis
#'
#' @author JAUNATRE Maxime
#'
#' @export
dataset = function(ind,popfile,entryfile,name,rare,qual,missLoci,missInd,LD) {

cat("\n");cat("reading file \n")
lecture <- readr::read_delim(entryfile,"\t", escape_double = FALSE, trim_ws = TRUE)

cat("\n");cat("reordering data \n")
reorder = subset_reorder(lecture, ind )

cat("\n");cat("deleting rare allele \n")
Rare = rare(reorder, rare  = rare, r= T)

cat("\n");cat("deleting poor quality SNP \n")
quality = Rare[which(Rare$QUAL >= qual),]

cat("\n");cat("applying treshold to missing data \n")
missingdata = tri(data = quality,n.r=missLoci,n.c = missInd, r = T)

cat("\n");cat("deleting loci with narrow positions \n")
manus = c(1)
for (i in 2:dim(missingdata)[1]){   if (missingdata$CHROM[i] == missingdata$CHROM[i-1] & (missingdata$POS[i]-missingdata$POS[i-1]) < LD ){} else {manus = c(manus,i)} }
final = missingdata[manus,]

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

cat("\n");cat("DONE :) \n\n")

return(fichier)
}
