#' calculating traits' weight for aggregate index
#'
#' @param pz the name of the breed
#'
#' @return return the weight of each trait
#' @export
#'
#' @examples
CalWeight <- function(pz){
#加载需要的参数
	Pars <- get(paste0(pz, 'phe'))
	prm <- get(paste0(pz, 'prm'))
	for(i in 1:length(prm)){assign(names(prm)[i], prm[i])}
	TRAITS <- Pars$trait

	##边际效益：性状变化一单位时收益的变化
	M_profit <- c()

	##original profit
	for(j in cost_items){
		assign(paste0(j, 1), get(j))
	}
	Lit = n_progeny

	## COST and INCOME per individual
	COST = jb_cost1 + (dam_cost1/(Lit*Pars$dProb) + sire_cost1/(Lit*Pars$sProb)) + other_cost1
	INCOME = ind_sale1 + cull_sale1/(Lit*Pars$dProb) + cull_sale1/(Lit*Pars$sProb)
	profit0 = INCOME-COST

	for(t in TRAITS){
		#获取表型变化一单位时各个经济参数的变化值
		for(j in cost_items){
			assign(paste0(j, 1), get(j) + bonus[[t]][j])
		}
		if('CZS' %in% TRAITS){
			Lit = Lit + 1
		}
		
		COST = jb_cost1 + (dam_cost1/(Lit*Pars$dProb) + sire_cost1/(Lit*Pars$sProb)) + other_cost1
		INCOME = ind_sale1 + cull_sale1/(Lit*Pars$dProb) + cull_sale1/(Lit*Pars$sProb)
		M_profit = c(M_profit, INCOME - COST - profit0) 
	}
	names(M_profit) = TRAITS

	##经济重要性=边际效益*遗传标准差
	w = M_profit*sqrt(Pars$var)

	##权重
	w = w/sum(abs(w))
	return(w)
}







