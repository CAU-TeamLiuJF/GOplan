#' core_breeding_process_simulation
#'
#' @param nvar the num of variable combination
#'
#' @return
#' @export
#'
#' @examples
CoreSim <- function(nvar){
  var_name = paste(VAR[nvar, ], collapse = '')
  assign('e', list(var_name = var_name, nvar = nvar), env = .GlobalEnv)
  ##模拟核心群选育##
  source(system.file("extra_code","23_core.R", package = "GOplan"))

  ##计算核心群平均遗传进展##
  CoreRes(var_name)
  
  dtG <- c(t(out_results[["index"]][Time, ] - out_results[["index"]][1, ]))   ##综合育种值进展
  dtPhe <- c(rbind(c(Phenos + Aphe$mean), out_results[["pop0_G"]][Time, seq(1, 2*length(Aphe$trait), 2)+1]))                  ##各性状预期表型
  dtVg <- c(t(out_results[["varg"]][Time, ] - out_results[["varg"]][1, ]))   ##遗传方差下降程度
  dtInb <- c(diff(range(out_results[["dfInb"]][, 1])), diff(range(out_results[["dfInb"]][, 2]))) ##近交增量
  return(list(dtG, dtInb, profit, dtPhe, dtVg))
  rm(out_results, pos = .GlobalEnv)
}
