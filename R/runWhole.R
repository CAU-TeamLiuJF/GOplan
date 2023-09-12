#' run whole system's breeding programs task
#' 
#' @param prm_path the full path to the prm file
#'
#' @return
#' @export
#'
#' @examples 
runWhole <- function(prm_path = prm_path){
	library(AlphaSimR)
	library(parallel)
	library(writexl)
	library(readxl)
	library(tidyverse)
	library(foreach)
	library(stringr)
	library(mlrMBO)
	
	knitStartTime <- Sys.time()
	options (warn = -1)
	task <- 2
	#code_path = find.packages('GOplan')
	#code_path = 'C:/Users/h/Desktop/hqq/code' #debug
	#for(i in c('task', 'code_path')){
		assign('task', task, envir = .GlobalEnv)
	#}
	### load prm files(.GlobalEnv)
	readPrm(prm_path, task)
	### load functions(.GlobalEnv)
	source(system.file("extra_code","Funs.R", package = "GOplan"),encoding = 'UTF-8')
	### construct founder pop(save as RData)
	Founders(Breed)

	TRAITS <- Aphe$trait
	if(nrow(VAR) > 0){
		res <- foreach(aa=1:nrow(VAR)) %do% WholeSim(aa)
		out <- data.frame(matrix(nrow = nrow(VAR), ncol = length(TRAITS) + 1))
	}else{
		res <- WholeSim(1)
		out <- data.frame(matrix(nrow = 1, ncol = length(TRAITS) + 1))
	} 
	if(nrow(VAR) > 0){
		varnames <- c()
		for(i in 1:nrow(VAR)){
			varnames <- c(varnames, paste(VAR[i, ], collapse = '-'))
			for(n in 1:(ncol(out)-1)){
				out[i, n] <- res[[i]][[1]][n]
			}
			out[i, ncol(out)] <- res[[i]][[2]]
		}
		out <- mutate(out, Var = varnames)
		names(out) <- c(TRAITS, 'Profit', paste(names(VAR), collapse = '-'))
	}else{
		for(n in 1:(ncol(out)-1)){
			out[1, n] <- res[[1]][n]
		}
		out[1, ncol(out)] <- res[[2]]
		names(out) <- c(TRAITS, 'Profit')
	}
	write_xlsx(out, paste0(out_path, "WholeOut.xlsx"))

	for(i in c(Breed, '1.RData', '2.RData', paste0(Breed, "founder.Rdata"))){
		unlink(paste0(out_path, i),recursive = T) 
	}
	
	cat("Analyse finish!\n")
	cat("Results generated in:", out_path, "\n")
	print(Sys.time() - knitStartTime)
}




