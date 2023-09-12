#' get whole breeding task's results
#'
#' @param var_name the var_name of variable combination
#'
#' @return
#' @export
#'
#' @examples
WholeRes <- function(var_name = var_name){
	for(pz in Breed){
	  out_list = get(paste0(pz, '_out'))
	  Trait <- get(paste0(pz,'Trait'))
	  file_name = c('select', 'pop0', 'my')
	  out_results = list()
	  #---------------------------------------------------------------------------------------
		out_varg = out_acc = data.frame(matrix(ncol = length(Trait), nrow = Time))
		names(out_acc) = names(out_varg) = Trait
		for(phe in Trait){
		  Data = Data1 = data.frame(rep(NA, Time))
		  for (j in 1:Nrep) {
			Data = cbind(Data, out_list[[j]][[paste0('VARG-pop0-', var_name)]][1:(Time), phe])
			Data1 = cbind(Data1, out_list[[j]][[paste0('acc-', var_name)]][1:Time, phe])
		  }
		  out_varg[phe] = apply(Data[, -1], 1, mean); out_acc[phe] = apply(Data1[, -1], 1, mean)
		}

		Data2 <- Data3 <- data.frame(rep(NA, Time))
		for (j in 1:Nrep) {
		  Data2 <- cbind(Data2, data.frame(out_list[[j]][[paste0("acc-index-", var_name)]][1:(Time), ]))
		  Data3 <- cbind(Data3, data.frame(out_list[[j]][[paste0("EBV-", var_name)]][1:(Time), 'g_index']))
		}
		out_acc2 <- apply(Data2[, -1], 1, mean)
		out_index <- data.frame('g_index' = apply(Data3[, -1], 1, mean))

		out_results[['varg']] = out_varg
		out_results[['acc']] = out_acc; out_results[['acc_index']] = out_acc2
		out_results[["index"]] <- out_index

	  for(name in file_name){
		for(phe in Trait){
		  Data = Data0 = Data1 = data.frame(rep(NA, Time))
		  for (j in 1:Nrep) {
			out_list = get(paste0(pz, '_out'))
			assign(paste0('P_G', j), out_list[[j]][[paste0('P_G-', name, '-', var_name)]][1:(Time), ])

			Data = cbind(Data, get(paste0('P_G', j))[, paste0('G_',phe)])
			Data0 = cbind(Data0, get(paste0('P_G', j))[, paste0('G0_',phe)])
			Data1 = cbind(Data1, get(paste0('P_G', j))[, paste0('G1_',phe)])
		  }
		  assign(phe, Data); assign(paste0(phe,0), Data0); assign(paste0(phe,1), Data1)
		}

		  for (trait in Trait) {
			data = get(trait)
			data = data[, -1] ##
			data$mean = apply(data,1,mean)
			assign(trait, data)
			data1 = get(paste0(trait, 1))
			data1 = data1[, -1]
			data1$mean = apply(data1,1,mean)
			assign(paste0(trait, 1), data1)
			data0 = get(paste0(trait, 0))
			data0 = data0[, -1]
			data0$mean = apply(data0,1,mean)
			assign(paste0(trait, 0), data0)
		  }
		  mean_G = data.frame(matrix(NA, ncol = length(Trait)*3, nrow = nrow(get(Trait[[1]]))))
		  n = length(Trait)
		  for(i in 1:n){
			a = get(Trait[[i]]); b = get(paste0(Trait[[i]], 1)); c = get(paste0(Trait[[i]], 0))
			mean_G[, i] = a$mean
			mean_G[, i+n] = b$mean
			mean_G[, i+2*n] = c$mean
		  }
		  assign(paste0("mean_G_", name), mean_G)
		}
	  out_results[[paste0(var_name, '_my_G')]] = mean_G_my
	  out_results[[paste0(var_name, '_select_G')]] = mean_G_select
	  out_results[[paste0(var_name, '_pop0_G')]] = mean_G_pop0
	  assign(paste0(pz, "_results"), out_results, envir = .GlobalEnv)
	}
	#save.image(file = paste0(out_path, '3.RData'))
	##################################################
	# predict finishers' phenotype
	##################################################
  final_out = list()
	list_SG = list()

	G = fun_gap(par1 = Ctype, par2 = Time+1)
	rrow = G[['rrow']]; rrow_P = G[['rrow_P']]
	for(i in c('A', 'B', 'C', 'D')){
	  #give an initial value (matrix 0)
	  dt = data.frame(matrix(0, ncol = ncol(A_results[[paste0(var_name, '_pop0_G')]]),
	                         nrow = nrow(A_results[[paste0(var_name, '_pop0_G')]])))
	  names(dt) = c(Trait, paste0(Trait, 1), paste0(Trait, 0))
	  assign(paste0('M', i, 'G'), dt)
	  assign(paste0(i, 'G'), dt)
	  assign(paste0(i, '0P'), dt)

	  if(exists(paste0(i, '_results'))){
	    results = get(paste0(i, '_results'))
	    a = get(paste0('M', i, 'S'))             ##the gene attribution vector of selection group
	    b = results[[paste0(var_name, '_pop0_G')]]
	    c = results[[paste0(var_name, '_my_G')]]
	    d = results[[paste0(var_name, '_select_G')]]
	    names(b) = names(c) = names(d) = c(Trait, paste0(Trait, 1), paste0(Trait, 0))
	    if(sum(a) == 0){           #assign genetic value of secondary group
	      assign(paste0('M', i, 'G'), b)
	    }else{
	      assign(paste0('M', i, 'G'), c)
	    }
	    assign(paste0(i, 'G'), d)  #assign genetic value of selection group
	    assign(paste0(i, '0P'), b)   #assign genetic value of mass
	  }}

	if(substr(Ctype, 1,1) == 3){
	  final_out[['AB']] <- data.frame(matrix(NA, ncol = 0, nrow = Time))
	  gAB = G[['gAB']]

	}else if(substr(Ctype, 1,1) == 4){
	  final_out[['AB']] <- final_out[['CD']] <- data.frame(matrix(NA, ncol = 0, nrow = Time))
	  gAB = G[['gAB']]; gCD = G[['gCD']]
	}

	 #-------------------------------------------------------------------------------
	final_out[['product']] <- data.frame(matrix(NA, ncol = 0, nrow = Time))
	final_out[['core']] = data.frame(matrix(NA, ncol = 0, nrow = Time))

	for (trait in Trait) {
	  #计算各个选择组每年的遗传优势
	  Pall = data.frame(matrix(NA, ncol = Time+2, nrow = nrow(AS)))
	  if(substr(Ctype, 1,1) == 2){
	    Phe0 = 0.5*(Bphe$mean[which(BTrait == trait)]+Aphe$mean[which(ATrait == trait)])
	    Pall[, 1] = c(rep(Aphe$mean[which(ATrait == trait)], As+1+Ad+1),
	                  rep(Bphe$mean[which(BTrait == trait)], Bs+1+Bd+1),
	                  rep(Aphe$mean[which(ATrait == trait)], Asec+1),
	                  rep(Bphe$mean[which(BTrait == trait)], Bsec+1),
	                  Phe0)
	  }else if(substr(Ctype, 1,1) == 3){
	    Phe0 = 0.25*(Bphe$mean[which(BTrait == trait)]+Aphe$mean[which(ATrait == trait)]+2*Cphe$mean[which(CTrait == trait)])
	    Pall[, 1] = c(rep(Aphe$mean[which(ATrait == trait)], As+1+Ad+1),
	                  rep(Bphe$mean[which(BTrait == trait)], Bs+1+Bd+1),
	                  rep(Cphe$mean[which(CTrait == trait)], Cs+1+Cd+1),
	                  rep(Aphe$mean[which(ATrait == trait)], Asec+1),
	                  rep(Bphe$mean[which(BTrait == trait)], Bsec+1),
	                  rep(Cphe$mean[which(CTrait == trait)], Csec+1),
	                  rep(0.5*(Bphe$mean[which(BTrait == trait)]+Aphe$mean[which(ATrait == trait)]), YABd+1),
	                  Phe0)
	  }else if(substr(Ctype, 1,1) == 4){
	    Phe0 = 0.25*(Bphe$mean[which(BTrait == trait)]+Aphe$mean[which(ATrait == trait)]+
	                   Cphe$mean[which(CTrait == trait)]+Dphe$mean[which(DTrait == trait)])
	    Pall[, 1] = c(rep(Aphe$mean[which(ATrait == trait)], As+1+Ad+1),
	                  rep(Bphe$mean[which(BTrait == trait)], Bs+1+Bd+1),
	                  rep(Cphe$mean[which(CTrait == trait)], Cs+1+Cd+1),
	                  rep(Dphe$mean[which(DTrait == trait)], Ds+1+Dd+1),
	                  rep(Aphe$mean[which(ATrait == trait)], Asec+1),
	                  rep(Bphe$mean[which(BTrait == trait)], Bsec+1),
	                  rep(Cphe$mean[which(CTrait == trait)], Csec+1),
	                  rep(Dphe$mean[which(DTrait == trait)], Dsec+1),
	                  rep(0.5*(Bphe$mean[which(BTrait == trait)]+Aphe$mean[which(ATrait == trait)]), YABd+1),
	                  rep(0.5*(Cphe$mean[which(CTrait == trait)]+Dphe$mean[which(DTrait == trait)]), YCDs+1),
	                  Phe0)
	  }

	  ##selection differential
	  ADG = BDG = CDG = DDG = ASG = BSG = CSG = DSG = A2G = B2G = C2G = D2G = c()
	  for (t in 1:Time) {                                    ##第t年的表型
	    ADG[t] = AG[t, paste0(trait, 0)] - A0P[t, paste0(trait, 0)]
	    ASG[t] = AG[t, paste0(trait, 1)] - A0P[t, paste0(trait, 1)]
	    A2G[t] = MAG[t, paste0(trait, 0)] - A0P[t, paste0(trait, 0)]
	    BDG[t] = BG[t, paste0(trait, 0)] - B0P[t, paste0(trait, 0)]
	    BSG[t] = BG[t, paste0(trait, 1)] - B0P[t, paste0(trait, 1)]
	    B2G[t] = MBG[t, paste0(trait, 1)] - B0P[t, paste0(trait, 1)]
	    CDG[t] = CG[t, paste0(trait, 0)] - C0P[t, paste0(trait, 0)]
	    CSG[t] = CG[t, paste0(trait, 1)] - C0P[t, paste0(trait, 1)]
	    if(substr(Ctype, 1,1) != 4){
	      C2G[t] = MCG[t, paste0(trait, 1)] - C0P[t, paste0(trait, 1)]
	    }else if(substr(Ctype, 1,1) == 4){
	      C2G[t] = MCG[t, paste0(trait, 0)] - C0P[t, paste0(trait, 0)]
	    }
	    DDG[t] = DG[t, paste0(trait, 0)] - D0P[t, paste0(trait, 0)]
	    DSG[t] = DG[t, paste0(trait, 1)] - D0P[t, paste0(trait, 1)]
	    D2G[t] = MDG[t, paste0(trait, 1)] - D0P[t, paste0(trait, 1)]

	    ##culmulative genetic merits
	    asg = bsg = csg = dsg = adg = bdg = cdg = ddg = mag = mbg = mcg = mdg = 0
	    for (j in 1:t) {
	      asg = asg + ASG[j]*AS[, t-j+1]
	      adg = adg + ADG[j]*AD[, t-j+1]
	      mag = mag + A2G[j]*MAS[, t-j+1]
	      bsg = bsg + BSG[j]*BS[, t-j+1]
	      bdg = bdg + BDG[j]*BD[, t-j+1]
	      mbg = mbg + B2G[j]*MBS[, t-j+1]
	      csg = csg + CSG[j]*CS[, t-j+1]
	      cdg = cdg + CDG[j]*CD[, t-j+1]
	      mcg = mcg + C2G[j]*MCS[, t-j+1]
	      ddg = ddg + DDG[j]*DD[, t-j+1]
	      dsg = dsg + DSG[j]*DS[, t-j+1]
	      mdg = mdg + D2G[j]*MDS[, t-j+1]
	    }
	    Pall[, t+1] = Pall[, 1] +
	      asg + bsg + csg + dsg + adg + bdg + cdg + ddg + mag + mbg + mcg + mdg
	  }

		#predict products' phenotype
	  final_out[['product']][, trait] <- c(t(Pall[rrow_P, 2:(Time+1)]))
	  #predict core phenotype
	  final_out[['core']][, paste0(trait, '_A')] = (unlist(Pall[rrow[1], 2:(Time+1)]) + unlist(Pall[rrow[3], 2:(Time+1)]))/2
	  final_out[['core']][, paste0(trait, '_B')] = (unlist(Pall[rrow[5], 2:(Time+1)]) + unlist(Pall[rrow[7], 2:(Time+1)]))/2
	  list_SG = list(ASG = AG-A0P, A2G = MAG-A0P, BSG = BG-B0P, B2G = MBG-B0P)
	  if(substr(Ctype, 1,1) == 3){
	    final_out[['AB']][, trait] = c(t((gAB %*% as.matrix(Pall[(G[['rrow_AB']]+1):(rrow_P-1), 2:(Time+1)]))*2))
	    final_out[['core']][, paste0(trait, '_C')] = (unlist(Pall[rrow[9], 2:(Time+1)]) + unlist(Pall[rrow[11], 2:(Time+1)]))/2
	    list_SG = list(ASG = AG-A0P, A2G = MAG-A0P, BSG = BG-B0P, B2G = MBG-B0P, CSG = CG-C0P, C2G = MCG-C0P)
	  }else if(substr(Ctype, 1,1) == 4){
	    final_out[['AB']][, trait] = c(t((gAB %*% as.matrix(Pall[(G[['rrow_AB']]+1):(G[['rrow_CD']]-1), 2:(Time+1)]))*2))
	    final_out[['CD']][, trait] = c(t((gCD %*% as.matrix(Pall[(G[['rrow_CD']]+1):(rrow_P-1), 2:(Time+1)]))*2))
	    final_out[['core']][, paste0(trait, '_C')] = (unlist(Pall[rrow[9], 2:(Time+1)]) + unlist(Pall[rrow[11], 2:(Time+1)]))/2
	    final_out[['core']][, paste0(trait, '_D')] = (unlist(Pall[rrow[13], 2:(Time+1)]) + unlist(Pall[rrow[15], 2:(Time+1)]))/2
	    list_SG = list(ASG = AG-A0P, A2G = MAG-A0P, BSG = BG-B0P, B2G = MBG-B0P,
	                   CSG = CG-C0P, C2G = MCG-C0P, DSG = DG-D0P, D2G = MDG-D0P)
	  }
	}
	#real phenotype of nuclei by simulaiton
	final_out[['A_core_out']] = A_results[[paste0(var_name, '_pop0_G')]]
	final_out[['B_core_out']] = B_results[[paste0(var_name, '_pop0_G')]]
	names(final_out[['A_core_out']]) = names(final_out[['B_core_out']]) = c(Trait, paste0(Trait, 1), paste0(Trait, 0))
	if(substr(Ctype, 1,1) == 3){
	  final_out[['C_core_out']] = C_results[[paste0(var_name, '_pop0_G')]]
	  names(final_out[['C_core_out']]) = c(Trait, paste0(Trait, 1), paste0(Trait, 0))
	}else if(substr(Ctype, 1,1) == 4){
	  final_out[['C_core_out']] = C_results[[paste0(var_name, '_pop0_G')]]
	  final_out[['D_core_out']] = D_results[[paste0(var_name, '_pop0_G')]]
	  names(final_out[['D_core_out']]) = names(final_out[['C_core_out']]) = c(Trait, paste0(Trait, 1), paste0(Trait, 0))
	}
	#############################################
	# calculate profit
	#############################################
	if(task != 3){
	  if(Ctype %in% c(31, 32, 33, 41, 42)){
	    NABd = NABs = NA2d*Aprm['n_progeny']*0.5*ABSR*ABProb;
	    NABd0 = NABd/ABProb; NA2d0 = NA2d/Aphe$secProb
	  }else{
	    NABd = NABs = NABd0 = 0
	  }
	}
	PF_out = data.frame(cost = NA, income = NA, profit = NA, N_all = NA, N_p = NA,  n_profit= NA)
	NpCd = NpCs = NpDd = NpDs = m_NpCd = m_NpCs = m_NpDd = m_NpDs = g_NpCd = g_NpCs = g_NpDd = g_NpDs = 0
	if(substr(Ctype, 1,1) != 2){
	  Lit = 0.5*(Aprm['n_progeny'] + Bprm['n_progeny'])
	}else{ Lit = Aprm['n_progeny']}
	##the number of survival offspring of breeds
	if(substr(Ctype, 1,1) == 3){
	  NpCd = NpCs = NCd*Cprm['n_progeny']*0.5*Cprm['SR'];           # NpCd/s: the number of survival offspring of breed C
	  m_NpCd = NCd*Cprm['nfam_F']; m_NpCs = NCd*Cprm['nfam_M']
	  g_NpCd = NCd*Cprm['nGeno_F']; g_NpCs = NCd*Cprm['nGeno_M']
	}else if(substr(Ctype, 1,1) == 4){
	  NpCd = NpCs = NCd*Cprm['n_progeny']*0.5*Cprm['SR'];
	  NpDd = NpDs = NDd*Dprm['n_progeny']*0.5*Dprm['SR'];
	  m_NpCd = NCd*Cprm['nfam_F']; m_NpCs = NCd*Cprm['nfam_M']
	  g_NpCd = NCd*Cprm['nGeno_F']; g_NpCs = NCd*Cprm['nGeno_M']
	  m_NpDd = NDd*Dprm['nfam_F']; m_NpDs = NDd*Dprm['nfam_M']
	  g_NpDd = NDd*Dprm['nGeno_F']; g_NpDs = NDd*Dprm['nGeno_M']
	}
	NpBd = NpBs = NBd*Bprm['n_progeny']*0.5*Bprm['SR'];           # NpBd/s: the number of survival offspring of breed B
	NpAd = NpAs = NAd*Aprm['n_progeny']*0.5*Aprm['SR'];

	for (t in 1:Time) {
	  #the phenotype progress in time t
	  for(i in c('A', 'B', 'C', 'D')){
	    if(i %in% Breed == T){
	      if('CZS' %in% names(final_out[['product']])){
	        assign(paste0(i, '_Lit'), final_out[[paste0(i, '_core_out')]][['CZS']][t])
	      }else{assign(paste0(i, '_Lit'), get(paste0(i, 'prm'))['n_progeny'])}
	    }else{
	      assign(paste0(i, '_Lit'), 0)}
	  }
	  Phenos <- c(t(final_out[['product']][t, ] - final_out[['product']][1, ]))
	  names(Phenos) <- names(final_out[['product']])
	  ##number of final survival products
	  if(substr(Ctype, 1,1) != 2){                      ##not Cytpe 2*
	    if('CZS' %in% names(final_out[['AB']])){
	      Lit = final_out[['AB']][['CZS']][t]
	    }
	    Nproduct = NABd*Lit*PSR
		N_wean = NDd*D_Lit+NCd*C_Lit+NBd*B_Lit+NAd*A_Lit+
			NA2d*A_Lit+NCd2*C_Lit+
			Nproduct/PSR
	  }else{      ##Cytpe 2*
	    if('CZS' %in% names(final_out[['product']])){
	      Lit = final_out[['A_core_out']][['CZS']][t]
	    }
	    Nproduct = NA2d*Lit*PSR
		N_wean = NDd*D_Lit+NCd*C_Lit+NBd*B_Lit+NAd*A_Lit+Nproduct/PSR
	  }
	  #number of total dams
	  N_dam = NDd+NCd+NBd+NAd+NA2d+NCd2+NABd
	  #number of stocks that need culled or repalced
	  NBd0 = NBd/Bphe$dProb; NBs0 = NBs/Bphe$sProb;
	  NAd0 = NAd/Aphe$dProb; NAs0 = NAs/Aphe$sProb
	  N_gxzz = NDd0+NDs0+NCd0+NCs0+NBd0+NBs0+NAd0+NAs0+NA2d0+NBs20+NCs20+NCd20+NDs20+NABd0+NCDs0
	  #number of mesuring or genotyping stocks
	  m_NpAd = NAd*Aprm['nfam_F']; m_NpAs = NAd*Aprm['nfam_M']
	  m_NpBd = NBd*Bprm['nfam_F']; m_NpBs = NBd*Bprm['nfam_M']
	  N_measure = m_NpAd+m_NpAs+m_NpBd+m_NpBs+m_NpCd+m_NpCs+m_NpDd+m_NpDs
	  if (Method == "ssgbup") {
	    g_NpAd = NAd*Aprm['nGeno_F']; g_NpAs = NAd*Aprm['nGeno_M']
	    g_NpBd = NBd*Bprm['nGeno_F']; g_NpBs = NBd*Bprm['nGeno_M']
	    N_genotype = g_NpAd+g_NpAs+g_NpBd+g_NpBs+g_NpCd+g_NpCs+g_NpDd+g_NpDs
	  }
	  #number of stocks that selled
	  N_damsale = (m_NpAd+m_NpBd+m_NpCd+m_NpDd) - (NDd0+NCd0+NBd0+NAd0+NA2d0+NCd20)
	  N_siresale = (m_NpAs+m_NpBs+m_NpCs+m_NpDs) - (NDs0+NCs0+NBs0+NAs0+NBs20+NCs20+NDs20)
	  N_zzsale = N_damsale+N_siresale     #measured but not selected
	  if(N_zzsale < 0){N_zzsale = 0}      #in case sometimes mutiplier need more
	  #number of individuals selled as products
	  N_production = Nproduct+NABs+NABd+NCDs+NCDd+(NpDd+NpCd+NpBd+NpAd)*2 - 
					max((N_measure+NABd0+NCDs0), N_gxzz)

	  ##calculate the bunus
	  for(j in cost_items){
	    aa <- get(j)
	    for(i in names(Phenos)){
	      b <- bonus[[i]][j]
	      aa <- aa + b*Phenos[i]
	    }
	    assign(paste0(j, 1), aa)
	  }

	  #------#
	  # COST #
	  #------#
	  all_jbcost = N_wean*jb_cost1
	  all_sirecost = (NAs+NBs+NCs+NDs+NBs2+NDs2+NCDs)*sire_cost1
	  all_damcost = N_dam*dam_cost1
	  all_measurecost = N_measure*meas_cost
	  if (Method == "ssgbup") {
	    all_geno_cost = N_genotype*geno_cost
	  } else {
	    all_geno_cost = 0
	  }
	  COST = all_jbcost+all_sirecost+all_measurecost+all_geno_cost+all_damcost+other_cost
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
	  profit = INCOME-COST
	  profit/N_production

	  # total cost、total income、total profit；total individual selled, number of products, profit per capita
	  PF_out[t, ] = c(COST, INCOME, profit, N_production, Nproduct, profit/Nproduct)
	}

	final_out[['PF_out']] = PF_out
	final_out[['product']] = as.data.frame(final_out[['product']])
	final_out[['AB']] = as.data.frame(final_out[['AB']])
	assign('final_out', final_out, envir = .GlobalEnv)
	sheets = list()
	for(i in names(final_out)){
	  sheets[[i]] = final_out[[i]]
	}
	if(task == 2){
	  write_xlsx(sheets, paste0(out_path, var_name, "out.xlsx"))
	  write_xlsx(list_SG, paste0(out_path, var_name, "SG_out.xlsx"))
	}
	sink(paste0(out_path, "detailInfo.txt"), append=TRUE)
	cat(var_name, "\n")
	cat("NAd:",NAd, "\t", "NAs:",NAs, "\t", "NAd2:", NA2d, "\t", "Asec:", Asec, "\t", "A_sor:",Aprm['Sor'], "\n")
	cat("NAd0:",NAd0, "\t", "NAs0:",NAs0, "\t", "NAd20:", NA2d0, "\n")
	cat("NBd:",NBd, "\t", "NBs:",NBs, "\t", "NBs2:", NBs2, "\t", "Bsec:", Bsec, "\t", "B_sor:",Bprm['Sor'], "\n")
	cat("NBd0:",NBd0, "\t", "NBs0:",NBs0, "\t", "NBs20:", NBs20, "\n")
	if(substr(Ctype, 1,1) != 2){
	  cat("NCd:",NCd, "\t",  "NCs:",NCs, "\t", "NCs2:", NCs2, "\t", "Csec:", Csec, "\t", "C_sor:",Cprm['Sor'], "\n")
	  cat("NCd0:",NCd0, "\t", "NCs0:",NCs0, "\t", "NCs20:", NCs20, "\n")
	  cat("NABd:",NABd, "\t", "YABd:",YABd, "\t", "SorAB:",SorAB, "\n")
	}
	if(substr(Ctype, 1,1) == 4){
	  cat("NDd:",NDd, "\t",  "NDs:",NDs, "\t", "NDs2:", NDs2, "\t", "Dsec:", Dsec, "\t", "D_sor:",Dprm['Sor'], "\n")
	  cat("NDd0:",NDd0, "\t", "NDs0:",NDs0, "\t", "NDs20:", NDs20, "\n")
	  cat("NCDs:",NCDs, "\t", "YCDs:",YCDs, "\t", "SorCD:", SorCD, "\n")
	}
	cat( "SorP:",SorP, "\n")
	cat("N_product:",Nproduct, "\n")
	sink()
}

