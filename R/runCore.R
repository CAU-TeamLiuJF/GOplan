#' run core breeding programs task
#' 
#' @param prm_path the full path to the prm file
#'
#' @return
#' @export
#'
#' @examples
runCore <- function(prm_path = prm_path){
	library(AlphaSimR)
	library(parallel)
	library(writexl)
	library(readxl)
	library(tidyverse)
	#library(magrittr)
	library(foreach)
	library(stringr)
	library(mlrMBO)
	
	knitStartTime <- Sys.time()
	options (warn = -1)
	task <- 1
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
		res <- foreach(aa=1:nrow(VAR)) %do% CoreSim(aa)
		out <- data.frame(matrix(nrow = nrow(VAR), ncol = 2*length(TRAITS) + 3))
	}else{
		res <- CoreSim(1)
		out <- data.frame(matrix(nrow = 1, ncol = 2*length(TRAITS) + 3))
	} 
	if(nrow(VAR) > 0){
		for(j in 1:3){    #c('G', 'Inb', 'Profit')
			varnames <- c()
			for(i in 1:nrow(VAR)){
				varnames <- c(varnames, paste(VAR[i, ], collapse = '-'))
				out[i, j] <- res[[i]][[j]]
			}
		}
		for(j in 4){
			for(i in 1:nrow(VAR)){
				out[i, j:ncol(out)] <- c(res[[i]][[j]], res[[i]][[j+1]])
			}
		}
		out <- mutate(out, Var = varnames)
		names(out) <- c('G', 'Inb', 'Profit', TRAITS, paste0('Vg_', TRAITS), paste(names(VAR), collapse = '-'))
		write_xlsx(out, paste0(out_path, "coreOut.xlsx"))
		
		df_varg <- df_acc <- df_acc_index <- df_index <- df_dfInb <- df_pop0_G <- data.frame(matrix(NA, ncol = 0, nrow = Time))
		out_list <- list()
		for(file in c('varg', 'acc', 'acc_index', 'index', 'dfInb', 'pop0_G')){
			df <- get(paste0('df_', file))
			for(i in 1:nrow(VAR)){
				var <- gsub("-", '', varnames[i])
				df <- cbind(df, get(paste0(var, '_out'))[[file]])
			}
			if(file %in% c('acc_index', 'index', 'dfInb')){ names(df) <- varnames }else{
			names(df) <- rep(varnames, each = length(TRAITS)) }
			out_list[[file]] <- df
		}
	}else{
		for(j in 1:3){    #c('G', 'Inb', 'Profit')
			out[1, j] <- res[[j]]
		}
		out[1, 4:ncol(out)] <- c(res[[4]], res[[5]])
		names(out) <- c('G', 'Inb', 'Profit', TRAITS, paste0('Vg_', TRAITS))
		write_xlsx(out, paste0(out_path, "coreOut.xlsx"))
		df_varg <- df_acc <- df_acc_index <- df_index <- df_dfInb <- df_pop0_G <- data.frame(matrix(NA, ncol = 0, nrow = Time))
		out_list <- list()
		for(file in c('varg', 'acc', 'acc_index', 'index', 'dfInb', 'pop0_G')){
			df <- get(paste0('df_', file))
			df <- cbind(df, get('_out')[[file]])
			out_list[[file]] <- df
		}
	}
	write_xlsx(out_list, paste0(out_path, "coreOut_detail.xlsx"))

	for(i in c(Breed, '1.RData', '2.RData', paste0(Breed, "founder.Rdata"))){
		unlink(paste0(out_path, i),recursive = T) 
	}
	
	cat("Analyse finish!\n")
	cat("Results generated in:", out_path, "\n")
	print(Sys.time() - knitStartTime)
}




