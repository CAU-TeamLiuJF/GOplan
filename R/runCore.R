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
	library(foreach)
	library(stringr)
	library(mlrMBO)

  for(i in c("AlphaSimR", "parallel", "writexl", "readxl", "tidyverse", "foreach", "stringr")){
    if(require(i) == FALSE){
      print(paste0("R package ", i, " does not exist, trying to install"))
      install.packages(i)
      if(require(i)){
        print(paste0("R package ", i, " successfully installed"))
      } else {
        stop(paste0("R package ", i, " installation failed, please install manually"))
      }
    }
  }
  
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
	if(createFounder == TRUE){
	  Founders(Breed)
	}

	TRAITS <- Aphe$trait
	if(nrow(VAR) > 0){
		res <- foreach(aa=1:nrow(VAR)) %do% CoreSim(aa)
		out <- data.frame(matrix(nrow = nrow(VAR), ncol = 4*length(TRAITS) + 5))
	}else{
		res <- CoreSim(1)
		out <- data.frame(matrix(nrow = 1, ncol = 4*length(TRAITS) + 5))
	}
	if(nrow(VAR) > 0){
		varnames <- c()
		for(i in 1:nrow(VAR)){
			varnames <- c(varnames, paste(VAR[i, ], collapse = '-'))
			out[i, ] <- c(res[[i]][[1]], res[[i]][[2]], res[[i]][[3]], res[[i]][[4]], res[[i]][[5]])
		
		}
		out <- mutate(out, Var = varnames)
		names(out) <- c('mean_G','sd_G', 'mean_Inb', 'sd_Inb', 'Profit', 
						c(paste0(c("mean","sd"), rep(TRAITS, each = 2))), 
						paste0('Vg_', c(paste0(c("mean","sd"), rep(TRAITS, each = 2)))), "var")
		write_xlsx(out, paste0(out_path, "coreOut.xlsx"))

		df_varg <- df_acc <- df_acc_index <- df_index <- df_dfInb <- df_pop0_G <- data.frame(matrix(NA, ncol = 0, nrow = Time))
		out_list <- list()
		for(file in c('varg', 'acc', 'acc_index', 'index', 'dfInb', 'pop0_G')){
			df <- get(paste0('df_', file))
			for(i in 1:nrow(VAR)){
				var <- gsub("-", '', varnames[i])
				df <- cbind(df, get(paste0(var, '_out'))[[file]])
			}
			if(file %in% c('acc_index', 'index', 'dfInb')){ names(df) <- paste0(c("mean","sd"), rep(varnames, each = 2)) }else{
			names(df) <- paste0(paste0(c("mean","sd"), "-", rep(TRAITS, each = 2)), "-",rep(varnames, each = length(TRAITS)*2)) }
			out_list[[file]] <- df
		}
	}else{
		out[1, ] <- c(res[[1]], res[[2]], res[[3]], res[[4]], res[[5]])
		names(out) <- c('mean_G','sd_G', 'mean_Inb', 'sd_Inb', 'Profit', 
						c(paste0(c("mean","sd"), rep(TRAITS, each = 2))), 
						paste0('Vg_', c(paste0(c("mean","sd"), rep(TRAITS, each = 2)))))
		write_xlsx(out, paste0(out_path, "coreOut.xlsx"))
		df_varg <- df_acc <- df_acc_index <- df_index <- df_dfInb <- df_pop0_G <- data.frame(matrix(NA, ncol = 0, nrow = Time))
		out_list <- list()
		for(file in c('varg', 'acc', 'acc_index', 'index', 'dfInb', 'pop0_G')){
			df <- get(paste0('df_', file))
			df <- cbind(df, get('_out')[[file]])
			out_list[[file]] <- df
		}
		if(file %in% c('acc_index', 'index', 'dfInb')){ names(df) <- c("mean","sd") }else{
		  names(df) <- paste0(paste0(c("mean","sd"), "-", rep(TRAITS, each = 2))) }
		out_list[[file]] <- df
	}
	write_xlsx(out_list, paste0(out_path, "coreOut_detail.xlsx"))

	for(i in c(Breed, '1.RData', '2.RData')){
		unlink(paste0(out_path, i),recursive = T)
	}

	cat("Analyse finish!\n")
	cat("Results generated in:", out_path, "\n")
	print(Sys.time() - knitStartTime)
}




