#' generate mate plans using MC method
#'
#' @param Nrows the number of plans file's rows
#' @param path the path of plans file
#'
#' @return 
#' @export an file contains mateplans and average inbreeding
#' @useDynLib GOplan
#'
#' @examples
runMC=function(Nrows = Nrows, path = getwd()){
  .C('mc', as.integer(Nrows), paste0(path, '/plans.txt'))
}
