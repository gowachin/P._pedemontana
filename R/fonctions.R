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
'%!in%' = function(x,y) {!('%in%'(x,y))}


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
