#' run whole breeding system's simulation
#'
#' @param nvar the num of variable combination
#'
#' @return
#' @export
#'
#' @examples
WholeSim <- function(nvar){
  var_name = paste(VAR[nvar, ], collapse = '')
  assign('e', list(var_name = var_name, nvar = nvar), env = .GlobalEnv)
  
  ##模拟核心群选育##
  source(system.file("extra_code","23_core.R", package = "GOplan"))

  ##计算核心群平均遗传进展##
  WholeRes(var_name)
  
  dtG <- final_out[['product']][Time, ]                    ##商品群表型进展
  profit <- mean(as.numeric(final_out[['PF_out']][, 'profit']))/1000000   ##平均年收益（百万）
  return(list(dtG, profit))
  rm(final_out, pos = .GlobalEnv)
}










