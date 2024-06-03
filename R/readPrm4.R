#' read parameter file
#'
#' @param prm_path path of parameter file
#' @param task number of analyse (1,2,3)
#'
#' @return
#' @export
#'
#' @examples
readPrm <- function(prm_path = prm_path, 
					task = task){
	prm <- read.table(prm_path, fileEncoding = 'UTF-8',fill = TRUE, blank.lines.skip = TRUE, col.names = paste0('V', 1:10))
	### read and assign Analyse Parmeters to gloablenv
	for(i in 1:which(prm$V1 == 'out_path')){
		assign(prm[i, 1], prm[i, 2], envir = .GlobalEnv)
	}
	for( i in c('Lit_pro', 'maxF', 'Nrep', 'Time', 'nSnpPerChr', 'nChr', 'nQtlPerChr')){
		 assign(i, as.numeric(get(i)), envir = .GlobalEnv) 
	}
	### create output path
	if (!file.exists(out_path)) {
		dir.create(out_path)
	}
	### read variable Parmeters and generate a combination dataframe
	Variable <- list()
	if(which(prm$V1 == 'Pars') > (which(prm$V1 == 'out_path')+1)){
		for(i in (which(prm$V1 == 'out_path')+1):(which(prm$V1 == 'Pars') - 1)){
			a <- as.numeric(t(prm[i, 2:ncol(prm)]))
			a <- a[is.na(a) == F]
			Variable[[prm[i, 1]]] <- a
		}
	}
	VAR <- expand.grid(Variable)

	### read Pop Parmeters and save as '*prm'
	for(i in c(which(prm$V1 == 'dam0'):which(prm$V1 == 'nIter'))){
		j <- as.numeric(t(prm[i, 2:ncol(prm)]))
		j <- j[is.na(j)==F]
		assign(prm[i, 1], j[j != ""], envir = .GlobalEnv)
	}
	Breed <- c('A', 'B', 'C', 'D')
	nBreed <- which(prm[which(prm$V1 == 'n_female'), ] == 0)[1]
	if(is.na(nBreed) == F){
		Breed <- Breed[1:(nBreed-2)]
	}
	
	if(task != 1 & substr(Ctype,1,1) != 2){
		cull_AB <- 1-cull_AB
		ABProb <- 1
		for(i in 1:(YABd-1)){ ABProb <- ABProb + prod(cull_AB[1:i]) }
		abPP <- c(1/ABProb)
		for(i in 1:(YABd-1)){ abPP <- c(abPP, prod(cull_AB[1:i])/ABProb) }
		assign('ABProb', ABProb, envir = .GlobalEnv);assign('abPP', abPP, envir = .GlobalEnv)
		assign('cull_AB', cull_AB, envir = .GlobalEnv);
		if(substr(Ctype,1,1) == 4){
		cull_CD <- 1-cull_CD
		CDProb <- 1
		for(i in 1:(YCDs-1)){ CDProb <- CDProb + prod(cull_CD[1:i]) }
		cdPP <- c(1/CDProb)
		for(i in 1:(YCDs-1)){ cdPP <- c(cdPP, prod(cull_CD[1:i])/CDProb) }
		assign('CDProb', CDProb, envir = .GlobalEnv);assign('cdPP', cdPP, envir = .GlobalEnv)
		assign('cull_CD', cull_CD, envir = .GlobalEnv);}
	}

	c <- c()
	for(pz in Breed){
		a <- as.numeric(prm[which(prm$V1 == 'n_female'):which(prm$V1 == 'nGeno_M'), 
							which(prm[prm$V1 == 'Breed', ] == pz)])
		names(a) <- prm[which(prm$V1 == 'n_female'):which(prm$V1 == 'nGeno_M'), 1]
		assign(paste0(pz, 'prm'), a, envir = .GlobalEnv)	
		if((EarlySelect == TRUE & Method != "blup")& (a['nfam_F']>a['nGeno_F'] | a['nfam_M']>a['nGeno_M'])){
			stop("Error: When EarlySelect is TRUE, nGeno_F/nGeno_M shouldn't less than nfam_F/nfam_M!!!")
		}
		## read each breed's detail phenotype information and save as '*phe'
		n = which(Breed == pz)
		b <- list() 
		for(i in which(prm$V1 == 'cull_s')[1]:which(prm$V1 == 'cull_sec')[1]){
			j <- as.numeric(t(prm[i, 2:ncol(prm)]))
			j <- j[is.na(j)==F]
			b[[prm$V1[i]]] <- j[j != ""]
		}
		j <- c(t(prm[which(prm$V1 == 'trait')[n], 2:ncol(prm)]))  ##get all traits name of breed
		j <- j[is.na(j)==F]
		b[['trait']] <- j[j != ""]
		c <- c(c, b[['trait']])                                   ##add breed's trait to the all trait vector
		for(i in (which(prm$V1 == 'SexLimit')[n]):(which(prm$V1 == 'analyse')[n])){
			j <- as.numeric(t(prm[i, 2:ncol(prm)]))
			j <- j[is.na(j)==F]
			b[[prm$V1[i]]] <- j[j != ""]
		}
		b$var = b$var*b$heri   ##change phenotype variance to gentic variance 
		j <- (which(prm$V1 == 'analyse')[n] + 1):(which(prm$V1 == 'analyse')[n] + length(b[['trait']]))
		j <- prm[j, 2:(1 + length(b[['trait']]))]
		##if multiple traits, generate corA matrix
		if(is.null(ncol(j)) == FALSE){
			for(i in 1:ncol(j)){
				j[, i] <- as.numeric(j[, i])
			}
			rownames(j) <- names(j)
			b[['cors']] <-  as.matrix(j)
		}
		for(i in c('cull_s', 'cull_d', 'cull_sec')){
			b[[i]] <- 1-b[[i]]
		}
		b$sProb = b$dProb = b$secProb = 1    ##calculate total proportion of all age classes (a num)
		for(k in c('s','d','sec')){
			if(task == 1 & k == 'sec' | a[[paste0('Y', k)]] <= 0){next}
			for(i in 1:(a[[paste0('Y', k)]]-1)){
				b[[paste0(k, 'Prob')]] = b[[paste0(k, 'Prob')]] + 
					prod(b[[paste0('cull_', k)]][1:i])}
		}
		for(k in c('s','d','sec')){       ##calculate separate gene proportion (a vector)
		if(task == 1 & k == 'sec' | a[[paste0('Y', k)]] <= 0){next}
		b[[paste0(k, 'PP')]] = c(1/b[[paste0(k, 'Prob')]])
			for(i in 1:(a[[paste0('Y', k)]]-1)){
				b[[paste0(k, 'PP')]] = c(b[[paste0(k, 'PP')]], 
					prod(b[[paste0('cull_', k)]][1:i])/b[[paste0(k, 'Prob')]])}
		}
		assign(paste0(pz, 'phe'), b, envir = .GlobalEnv)
		
	}

	## read initial economic parameters of each item
	for(i in which(prm$V1 == 'jb_cost')[1]:which(prm$V1 == 'geno_cost')[1]){
		assign(prm[i, 1], as.numeric(prm[i, 2]), envir = .GlobalEnv)
	}
	## read bonus of each trait
	bonus <- list()
	cost_items <- c(t(prm[which(prm$V1 == 'jb_cost'):which(prm$V1 == 'dam_sale'), 1]))
	c <- unique(c)
	for(i in c){
		j <- as.numeric(t(prm[which(prm$V1 == 'jb_cost'):which(prm$V1 == 'dam_sale'), 
									which(c == i)+2]))
		names(j) <- cost_items
		bonus[[i]] <- j
	}
	##------------------------------------------------------------------------------------##
	##read information of optimization
	if(task == 3){
		Parset <- prm[which(prm$V1 == 'Pars'):(which(prm$V1 == 'dam0') - 1), ]
		Parset <- Parset[, which(is.na(Parset[1, ]) == F)]
		names(Parset) <- c(t(Parset[1, ]))
		rownames(Parset) <- Parset$Pars
		Parset <- Parset[-1, ]
		Parset$lower <- as.integer(Parset$lower)
		Parset$upper <- as.integer(Parset$upper)
		YAB_max <- Parset['YABd', 'upper']
		YAB_min <- Parset['YABd', 'lower']
		Lit_AB = 0.5*(Aprm['n_progeny'] + Bprm['n_progeny'])
		NABd <- NABs <- N_P/PSR/Lit_AB
		if(substr(Ctype,1,1) == 3){    ##1 is the proportion of select Asec/Aprogeney
			NA_min <- NABd/(YAB_max*ABSR*(Aprm['n_progeny']*0.5)*  ## = NA2d
						Parset['Asec', 'upper']*1*                 ## = NA2d0
						(Aprm['nfam_F']-1/Aprm['Yd']))
			if(Ctype == 33){
			a <- NABd/(YAB_max*ABSR*(Aprm['n_progeny']*0.5)*   ## = NA2d
						SorAB*Parset['Bsec', 'upper'])         ## a=NBs20
			NB_min <- a/(1*(Bprm['nfam_F']-1/Bprm['Sor']/Bprm['Ys']))
			NC_min <- NABd/(SorP*Parset['Csec', 'upper']*1*
						(Cprm['nfam_F']-1/Cprm['Sor']/Cprm['Ys']))
			}else if(Ctype == 32){
			a <- NABd/(YAB_max*ABSR*(Aprm['n_progeny']*0.5)*   ## = NA2d
						SorAB*Bprm['Ys'])                      ## a=NBs20
			NB_min <- a/(1*(Bprm['nfam_F']-1/Bprm['Sor']/Bprm['Ys']))
			NC_min <- NABd/(SorP*Parset['Csec', 'upper']*1*
						(Cprm['nfam_F']-1/Cprm['Sor']/Cprm['Ys']))
			}else if(Ctype == 31){
			a <- NABd/(YAB_max*ABSR*(Aprm['n_progeny']*0.5)*   ## = NA2d
						SorAB*Bprm['Ys'])         ## a=NBs20
			NB_min <- a/(1*(Bprm['nfam_F']-1/Bprm['Sor']/Bprm['Ys']))
			NC_min <- NABd/(SorP*Cprm['Cs']*1*
						(Cprm['nfam_F']-1/Cprm['Sor']/Cprm['Ys']))
			}
			for(i in c('Parset', 'YAB_max', 'YAB_min', 'Lit_AB', 'NABd', 'NABs', 'NA_min', 'NB_min', 'NC_min')){
				assign(i, get(i), envir = .GlobalEnv)}
		}else if(substr(Ctype,1,1) == 4){
			YCD_max <- Parset['YCDs', 'upper']
			YCD_min <- Parset['YCDs', 'lower']
			NA_min <- NABd/(YAB_max*ABSR*(Aprm['n_progeny']*0.5)*  ## = NA2d
						Parset['Asec', 'upper']*1*                 ## = NA2d0
						(Aprm['nfam_F']-1/Aprm['Yd']))
			NC_min <- NCDs/(YCD_max*CDSR*(Cprm['n_progeny']*0.5)*  ## = NC2d
						Parset['Csec', 'upper']*1*                 ## = NC2d0
						(Cprm['nfam_F']-1/Cprm['Yd']))
			NCDs = NCDd = NABd/SorP
			if(Ctype == 42){
			a <- NABd/(YAB_max*ABSR*(Aprm['n_progeny']*0.5)*   ## = NA2d
						SorAB*Parset['Bsec', 'upper'])         ## a=NBs20
			NB_min <- a/(1*(Bprm['nfam_F']-1/Bprm['Sor']/Bprm['Ys']))
			
			ND_min <- NCDs/(YCD_max*CDSR*(Cprm['n_progeny']*0.5)*
							SorCD*Parset['Dsec', 'upper']*         ## a=NDs20
							1*(Dprm['nfam_F']-1/Dprm['Sor']/Dprm['Ys']))			
			}else if(Ctype == 41){
			a <- NABd/(YAB_max*ABSR*(Aprm['n_progeny']*0.5)*   ## = NA2d
						SorAB*Bprm['Ys'])        
			NB_min <- a/(1*(Bprm['nfam_F']-1/Bprm['Sor']/Bprm['Ys']))			
			D_min <- NCDs/(YCD_max*CDSR*(Cprm['n_progeny']*0.5)*
							SorCD*Dprm['Ys']*         ## a=NDs20
							1*(Dprm['nfam_F']-1/Dprm['Sor']/Dprm['Ys']))
			}
			for(i in c('Parset', 'YAB_max', 'YAB_min', 'Lit_AB', 'NABd', 'NABs', 'NA_min', 'NB_min', 'NC_min',
						'YCD_max', 'YCD_min', 'ND_min', 'NCDs', 'NCDd')){
				assign(i, get(i), envir = .GlobalEnv)}
		}
		## judge whether the input parameter is suitable or give warnings
		for(i in c('NA', 'NB', 'NC', 'ND')){
			if((substr(Ctype,1,1) == 3 & i == 'ND') == FALSE){
			cat(i, '_min：', get(paste0(i, '_min')), '\n')
			if(get(paste0(i, '_min')) >= Parset[paste0(i, 'd'), 'upper']){
				stop("The number of ", gsub("N", '', i), " pop should larger than: ", get(paste0(i, '_min')), " !\n")
			}else if(get(paste0(i, '_min')) >= 0.75*Parset[paste0(i, 'd'), 'upper']){
				cat("Warning: We suggest to expand the range of ",i, "_min!\n")
			}else if(get(paste0(i, '_min')) < Parset[paste0(i, 'd'), 'lower']){
				assign(paste0(i, '_min'), Parset[paste0(i, 'd'), 'lower'], envir = .GlobalEnv)
			}}
		}
	}
	for(i in c('Nrep', 'Time','Breed', 'VAR','cost_items','bonus')){
		assign(i, get(i), envir = .GlobalEnv)
	}
}



