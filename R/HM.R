#' run HM mate
#'
#' @param NULL
#'
#' @return a list of crossplan and inbreeding values
#' @export
#'
#' @examples
HM <- function(Mate = Mate, path = 'plans.txt', maxF = maxF, Sor = Sor){
	plans = read.table(path, header = F)
	names(plans) = c('id_b', 'id_s', 'A', 'b_index', 's_index')
	if(Mate == 'HOMO'){
		plans$bs_index = 0.5*plans$b_index + 0.5*plans$s_index
		plans = plans[order(plans$bs_index, decreasing = T), ]  #reorder based Homogeneous
	}else{
		plans$bs_index = abs(plans$b_index - plans$s_index)
		plans = plans[order(plans$bs_index, decreasing = T), ]  #reorder based Heterogeneous
	}

	if(maxF == 999){
		eligible = plans
	}else{
		eligible = plans[plans$A <= 2*maxF, ]               #select eligible pairs
		re_plans = plans[plans$A > 2*maxF, ]		        #not eligible pairs
	}

	allBoar = unique(plans$id_b)
	allSow = unique(plans$id_s)
	useB = rep(0, length(unique(plans$id_b)))               #the usage of boar
	names(useB) = unique(plans$id_b)
	crossplan = matrix(allSow, nrow = length(allSow), ncol = 1)
	rownames(crossplan) = allSow
	crossplan = cbind(crossplan, rep(NA, length(allSow)))
	
	## assign boar for sow
	for(i in allSow){
		for(boar in unique(eligible[eligible$id_s == i, 'id_b'])){
			if(useB[as.character(boar)] < Sor){             #check the avaliability of boar
				crossplan[as.character(i), 2] = boar
				useB[as.character(boar)] = useB[as.character(boar)] + 1
				break
			}
		}
		if(is.na(crossplan[as.character(i), 2]) == T){      #if eligible do not have proper boar to assign, search in re_plans
			for(boar in unique(re_plans[re_plans$id_s == i, 'id_b'])){
				if(useB[as.character(boar)] < Sor){
					crossplan[as.character(i), 2] = boar
					useB[as.character(boar)] = useB[as.character(boar)] + 1
					break
				}
			}
		}
	}
	E = merge(data.frame(crossplan), plans, by.x = c('X1','X2'), by.y = c('id_s','id_b'), sort = F, all.x = T)
	list(crossplan = crossplan, Inb = mean(E$A)/2)
}
