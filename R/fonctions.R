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
