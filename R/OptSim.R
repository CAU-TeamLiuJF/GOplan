#' simulation of breeding process for optimization
#'
#' @param Asec service time of A's secondary group
#' @param Bsec service time of B's secondary group
#' @param Csec service time of C's secondary group
#' @param YABd service time of hybrids AB 
#' @param NAd service number of nucleus A's female stocks
#' @param NBd service number of nucleus B's female stocks
#' @param NCd service number of nucleus C's female stocks
#'
#' @return the breeding sysytem's profit of variable combination
#' @export
#'
#' @examples
OptSim <- function(x){
  if(Ctype == 33){
	Asec = x[1]; Bsec = x[2]; Csec = x[3]; YABd = x[4]
	NAd = x[5]; NBd = x[6]; NCd = x[7]
	var_name <- paste(Asec, Bsec, Csec, YABd, NAd, NBd, NCd, sep = '-')
	assign('e', list(Bsec = Bsec, Csec = Csec, Asec = Asec, YABd = YABd, #ABSor = ABSor, ABCSor = ABCSor,
					NAd = NAd, NBd = NBd, NCd = NCd,var_name = var_name), env = .GlobalEnv)
  }
  if(Ctype == 32){
	Asec = x[1]; Csec = x[2]; YABd = x[3]
	NAd = x[4]; NBd = x[5]; NCd = x[6]
	var_name <- paste(Asec, Csec, YABd, NAd, NBd, NCd, sep = '-')
	assign('e', list(Csec = Csec, Asec = Asec, YABd = YABd, #ABSor = ABSor, ABCSor = ABCSor,
					NAd = NAd, NBd = NBd, NCd = NCd,var_name = var_name), env = .GlobalEnv)
  }
  if(Ctype == 31){
	Asec = x[1]; YABd = x[2]
	NAd = x[3]; NBd = x[4]; NCd = x[5]
	var_name <- paste(Asec, YABd, NAd, NBd, NCd, sep = '-')
	assign('e', list(Asec = Asec, YABd = YABd, #ABSor = ABSor, ABCSor = ABCSor,
					NAd = NAd, NBd = NBd, NCd = NCd,var_name = var_name), env = .GlobalEnv)
  }
  if(Ctype == 42){
	Asec = x[1]; Bsec = x[2]; Csec = x[3]; Dsec = x[4]; YABd = x[5]; YCDs = x[6]
	NAd = x[7]; NBd = x[8]; NCd = x[9]; NDd = x[10]
	var_name <- paste(Asec, Bsec, Csec, Dsec, YABd, YCDs, NAd, NBd, NCd, NDd, sep = '-')
	assign('e', list(Dsec = Dsec, Bsec = Bsec, Csec = Csec, Asec = Asec, YABd = YABd, YCDs = YCDs, #ABSor = ABSor, ABCSor = ABCSor,
					NAd = NAd, NBd = NBd, NCd = NCd, NDd = NDd,var_name = var_name), env = .GlobalEnv)
  }
  if(Ctype == 41){
	Asec = x[1]; Csec = x[2]; YABd = x[3]; YCDs = x[4]
	NAd = x[5]; NBd = x[6]; NCd = x[7]; NDd = x[8]
	var_name <- paste(Asec, Csec, YABd, YCDs, NAd, NBd, NCd, NDd, sep = '-')
	assign('e', list(Csec = Csec, Asec = Asec, YABd = YABd, YCDs = YCDs, #ABSor = ABSor, ABCSor = ABCSor,
					NAd = NAd, NBd = NBd, NCd = NCd, NDd = NDd,var_name = var_name), env = .GlobalEnv)
  }
  ##模拟核心群选育##
  source(system.file("extra_code","OptCore.R", package = "GOplan"))

  ##计算核心群平均遗传进展##
    if(optCode == 0){
		WholeRes(var_name)
		gc()
		profit <- final_out[['PF_out']][Time, 'profit']/N_P
    }else{
		profit <- 0
    }
	names(profit) <- 'profit'
	profit
}










