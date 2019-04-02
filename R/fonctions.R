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

