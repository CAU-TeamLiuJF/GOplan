#' create founder populations for breeds
#'
#' @param Breed a vector of breeds' name
#'
#' @return
#' @export
#'
#' @examples
Founders <- function(Breed = Breed){
	for(pz in Breed){
	nfounder = ceiling(get(paste0(pz, 'prm'))[['n_female']]/1000)*200
	founderPop <- runMacs(nInd = nfounder, nChr = nChr, segSites = 2*nSnpPerChr)
	SP <- SimParam$new(founderPop)
	assign('SP', SP, envir = .GlobalEnv)
	Pars <- get(paste0(pz, 'phe'))

	if (length(Pars$trait) == 1) {
	  SP$addTraitA(
		nQtlPerChr = nQtlPerChr,
		mean =  Pars$mean,
		var = Pars$var,
		name = Pars$trait
	  )
	} else {
	  SP$addTraitA(
		nQtlPerChr = nQtlPerChr,
		mean = Pars$mean,
		var = Pars$var,
		corA = Pars$cors,
		name = Pars$trait
	  )
	}
	SP$addSnpChip(nSnpPerChr = nSnpPerChr, minSnpFreq = 0.05)
	SP$setSexes("yes_sys")

	heri <- Pars$heri
	# 建立基础群
	pop <- newPop(founderPop)
	pop <- setPheno(pop, h2 = heri)

	pop_all <- pop

	Ped <- data.frame(getPed(pop), G = 0)
	Phen <- data.frame(matrix(ncol = 3 + length(Pars$trait), nrow = 0))
	names(Phen) <- c("id", "sex", "G", Pars$trait)

	## 随机交配1-4代的构建(all in all out)
	for (g in 1:5) {
	  if (g == 1) {
		pop0 <- randCross(pop, nCrosses = length(pop@id) / 2, nProgeny = 4) # 到g=5的时候，pop0=3200male&3200female
		pop0 <- setPheno(pop0, h2 = heri)
	  } else {
		dam <- pop0[pop0@sex == "F"]
		sire <- pop0[pop0@sex == "M"]

		parents <- c(dam, sire)
		pop0 <- randCross(parents, nCrosses = length(dam@id), nProgeny = 4)
		pop0 <- setPheno(pop0, h2 = heri)

		pop_all <- c(parents, pop_all)
		Ped <- rbind(Ped, data.frame(getPed(parents), G = g - 1))
		Phen <- rbind(Phen, data.frame(
		  id = parents@id,
		  sex = parents@sex,
		  G = g - 1,
		  pheno(parents)
		))
	  }
	}
	##########################################################################
	pop00 <- pop0
	Ped00 <- Ped
	Phen00 <- Phen
	map = getSnpMap(snpChip=1, sex = "A")
	#-----------------------------------------------------------------------------------
	save(list=c('pop00', 'Ped00', 'Phen00','g', 'heri', 'Pars','SP'), file = paste0(out_path, pz, 'founder.Rdata'))
	}
}
