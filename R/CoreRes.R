#' cal_core_simulation_results
#'
#' @param var_name name of variable combination
#'
#' @return
#' @export
#'
#' @examples
CoreRes <- function(var_name = var_name){
	Trait = TRAITS
	file_name <- c("select", "pop0")
	out_results <- list()
	#---------------------------------------------------------------------------------------
	out_varg <- out_acc <- data.frame(matrix(ncol = 2*length(Trait), nrow = Time))
	names(out_acc) <- names(out_varg) <- c(paste0(c("mean","sd"), rep(Trait, each = 2)))
	for (phe in Trait) {
	  Data <- Data1 <- data.frame(rep(NA, Time))
	  for (j in 1:Nrep) {
			Data = cbind(Data, out_list[[j]][[paste0('VARG-pop0-', var_name)]][1:(Time), phe])
			Data1 = cbind(Data1, out_list[[j]][[paste0('acc-', var_name)]][1:Time, phe])
	  }
	  out_varg[paste0("mean", phe)] <- apply(Data[, -1], 1, mean)
	  out_varg[paste0("sd", phe)] <- apply(Data[, -1], 1, sd)
	  out_acc[paste0("mean", phe)] <- apply(Data1[, -1], 1, mean)
	  out_acc[paste0("sd", phe)] <- apply(Data1[, -1], 1, sd)
	}
	out_results[["varg"]] <- out_varg
	out_results[["acc"]] <- out_acc

	Data2 <- Data3 <- Data4 <- data.frame(rep(NA, Time))
	for (j in 1:Nrep) {
	  Data2 <- cbind(Data2, data.frame(out_list[[j]][[paste0("acc-index-", var_name)]][1:(Time), ]))
	  Data3 <- cbind(Data3, data.frame(out_list[[j]][[paste0("EBV-", var_name)]][1:(Time), 'g_index']))
	  Data4 <- cbind(Data4, out_list[[j]][[paste0("dfInb-", var_name)]][1:(Time)])
	}
	out_acc2 <- data.frame('acc_index' = apply(Data2[, -1], 1, mean))
	out_acc2$sd <- apply(Data2[, -1], 1, sd)
	out_index <- data.frame('g_index' = apply(Data3[, -1], 1, mean))
	out_index$sd <- apply(Data3[, -1], 1, sd)
	out_Inb <- data.frame('Inb' = apply(Data4[, -1], 1, mean))
	out_Inb$sd <- apply(Data4[, -1], 1, sd)

	out_results[['acc_index']] = out_acc2
	out_results[["index"]] <- out_index
	out_results[["dfInb"]] <- out_Inb

	for (name in file_name) {
	for (phe in Trait) {
	  Data <- Data0 <- Data1 <- data.frame(rep(NA, Time))
		for (j in 1:Nrep) {
		  Data <- cbind(Data, data.frame(out_list[[j]][[paste0("P_G-", name, "-", var_name)]][1:(Time), paste0("G_", phe)]))
		  Data0 <- cbind(Data0, data.frame(out_list[[j]][[paste0("P_G-", name, "-", var_name)]][1:(Time), paste0("G0_", phe)]))
		  Data1 <- cbind(Data1, data.frame(out_list[[j]][[paste0("P_G-", name, "-", var_name)]][1:(Time), paste0("G1_", phe)]))
		}
	  assign(phe, Data)
	  assign(paste0(phe, 0), Data0)
	  assign(paste0(phe, 1), Data1)
	}

	for (trait in Trait) {
	  data <- get(trait)
	  data <- data[, -1] ##
	  data$mean <- apply(data, 1, mean)
	  data$sd <- apply(data, 1, sd)
	  assign(trait, data)
	  data1 <- get(paste0(trait, 1))
	  data1 <- data1[, -1]
	  data1$mean <- apply(data1, 1, mean)
	  data1$sd <- apply(data1, 1, sd)
	  assign(paste0(trait, 1), data1)
	  data0 <- get(paste0(trait, 0))
	  data0 <- data0[, -1]
	  data0$mean <- apply(data0, 1, mean)
	  data0$sd <- apply(data0, 1, sd)
	  assign(paste0(trait, 0), data0)
	}
	mean_G <- data.frame(matrix(NA, ncol = 2*length(Trait), nrow = nrow(get(Trait[[1]]))))
	names(mean_G) <- c(paste0(c("mean","sd"), rep(Trait, each = 2)))
	n <- length(Trait)
	for (i in 1:n) {
	  a <- get(Trait[[i]])
	  mean_G[, 2*i-1] <- a$mean
	  mean_G[, 2*i] <- a$sd
	} 
	assign(paste0("mean_G_", name), mean_G)
	}
	out_results[["select_G"]] <- mean_G_select
	out_results[["pop0_G"]] <- mean_G_pop0

	#write_xlsx(out_results, paste0(out_path, var_name, "out.xlsx"))
	assign(paste0(var_name, '_out'), out_results, envir = .GlobalEnv)
	assign('out_results', out_results, envir = .GlobalEnv)
	#save.image(file = paste0(outpath, '/3.RData'))
	#############################################
	#本代码用于预估育种投入
	#############################################
	profit <- c()
	for(t in c(1:Time)){

	  #各性状第Time年预期表型进展
	  Phenos <- c(t(out_results[["pop0_G"]][t, seq(1, n, 2)] - out_results[["pop0_G"]][1, seq(1, n, 2)]))
	  names(Phenos) <- Trait
	  ##计算当前表型下各个经济参数的值
	  for(j in cost_items){
		aa <- get(j)
		for(i in names(Phenos)){
		  b <- bonus[[i]][j]
		  aa <- aa + b*Phenos[i]
		}
		assign(paste0(j, 1), aa)
	  }
	  
	  Lit = n_progeny
	  if('CZS' %in% Trait){
		Lit = Lit + Phenos['CZS']
	  }
	  
	  Npd = Nps = n_female*Lit*0.5*SR;           # Npdd/s: 核心群后代公母猪存活数目
	  #每个季节全群断奶仔猪数目
	  N_wean = n_female*Lit
	  #更新所需种猪数目=淘汰种猪数目
	  Nd0 = nSelected_F; Ns0 = nSelected_M
	  N_gxzz = Nd0+Ns0
	  #入测种猪数目（窝*每窝测定个体）
	  m_Npd = n_female*nfam_F; m_Nps = n_female*nfam_M
	  N_measure = m_Npd+m_Nps
	  if (Method == "ssgbup") {
		g_Npd = n_female*nGeno_F; g_Nps = n_female*nGeno_M
		N_genotype = g_Npd+g_Nps
	  }
	  #未留种入测种猪（售卖）
	  N_damsale = m_Npd-Nd0
	  N_siresale = m_Nps-Ns0
	  N_zzsale = N_damsale+N_siresale     #入测后未选留的
	  #做商品猪售卖的数目
	  N_production = Npd*2-N_measure                
	  #              没有入测的=所有活猪-入测猪
	  #------#
	  # COST #
	  #------# 
	  all_jbcost = N_wean*jb_cost1
	  all_sirecost = n_male*sire_cost1
	  all_damcost = n_female*dam_cost1
	  all_measurecost = N_measure*meas_cost
	  if (Method == "ssgbup") {
		all_geno_cost = N_genotype*geno_cost
	  } else {
		all_geno_cost = 0
	  }
	  COST = all_jbcost+all_sirecost+all_measurecost+all_geno_cost+all_damcost
	  #--------#
	  # INCOME #  
	  #--------#
	  in_product = N_production*ind_sale1
	  in_zzsale = N_damsale*dam_sale1+N_siresale*sire_sale1
	  in_outpig = N_gxzz*cull_sale1
	  INCOME = in_product+in_zzsale+in_outpig
	  #--------#
	  # PROFIT #
	  #--------#
	  profit = c(profit, INCOME-COST) 
	}
	profit <- mean(profit)
	names(profit) <- c('mean_profit')
	assign('profit', profit, envir = .GlobalEnv)
	assign('Phenos', Phenos, envir = .GlobalEnv)
	
	sink(paste0(out_path, "detailInfo.txt"), append=TRUE) 
	cat(var_name, "\n")
	cat("NAd:",n_female, "\t", "NAs:",n_male, "\t", "A_sor:", Sor, "\n")
	cat("YAd:",Yd, "\t", "YAs:",Ys, "\n")
	cat("NAd0:",nSelected_F, "\t", "NAs0:",nSelected_M, "\n")
	cat("nfam_F:",nfam_F, "\t", "nfam_M:",nfam_M, "\n")
	if(Method == 'ssgbup'){cat("nGeno_F:",nGeno_F, "\t", "nGeno_M:",nGeno_M, "\n")}
	sink()
}
