library(AlphaSimR); library(tidyverse);library(stringr)
  load(paste0(out_path, '1.RData'))##加载之前的工作路径
  Ped = Ped00 
  Phen = Phen00  #沿用之前的
  #-------------------------------------------------------#  
  ##一些待写出文件
  dfInb <- rep(NA, Time)
  P_G = P_G2 = P_G3 = data.frame(matrix(ncol = length(TRAITS)*6, nrow = 0))
  names(P_G) = names(P_G2) = names(P_G3) = c(paste0('P_', TRAITS), paste0('G_', TRAITS),
			   paste0('P1_', TRAITS), paste0('G1_', TRAITS),
			   paste0('P0_', TRAITS), paste0('G0_', TRAITS))				 
  VARG = VARG2 = data.frame(matrix(ncol = length(TRAITS), nrow = 0))
  names(VARG) = names(VARG2) = c(TRAITS)  
  EBV = data.frame(matrix(ncol = length(TRAITS)+1, nrow = 0))
  names(EBV) = c(paste0('g_', TRAITS), 'g_index')
  acc = data.frame(matrix(ncol = length(TRAITS), nrow = 0))
  names(acc) = c(TRAITS)
  acc2 = data.frame(matrix(ncol = 1, nrow = 0)); names(acc2) = c('index')
  ebv_all = data.frame(matrix(ncol = 2, nrow = 0))
  names(ebv_all) = c('Id', 'Index')
  #-------------------------------------------------------------------------------------------------
  rdmu_path = paste0(out_path, pz, '/', x)
  if(file.exists(rdmu_path) == FALSE){
    dir.create(rdmu_path) 
  }
  
  setwd(rdmu_path)
  #-------------------------------------------------------------------------------------------------
  
  famselect <- c(
        selectWithinFam(pop = pop00, nInd = 1, use = "rand", sex = "F"),
        selectWithinFam(pop = pop00, nInd = 1, use = "rand", sex = "M"))

      parent <- c(
        selectInd(famselect, nInd = nSelected_F, sex = "F", use = "rand"),
        selectInd(famselect, nInd = nSelected_M, sex = "M", use = "rand")
      )
      Ped <- rbind(Ped, data.frame(getPed(famselect), G = g))
      Phen <- rbind(Phen, data.frame(
        id = famselect@id,
        sex = famselect@sex,
        G = g,
        pheno(famselect)
      ))
	  famselect = famselect[famselect@id[which(famselect@id %in% parent@id == FALSE)]]

      ## 构建核心群（母4年龄，公2年龄）
      for (g in 6:(7+Yd)) {
        pop_select <- c(
          selectInd(famselect, nInd = nSelected_F, sex = "F", use = "rand"),
          selectInd(famselect, nInd = nSelected_M, sex = "M", use = "rand")
        )
		ebv_all <- rbind(ebv_all,
			data.frame(Id = pop_select@id, Index = runif(length(pop_select@id))))
        mateRes = Mate_rand()
        pop0 <- makeCross(parent, crossPlan = mateRes[['crossplan']], nProgeny = n_progeny)
        pop0 <- setPheno(pop0, h2 = heri)

        if (g == (7+Yd)) {
          famselect <- c(
            selectWithinFam(pop = pop0, nInd = nfam_F, use = "rand", sex = "F"),
            selectWithinFam(pop = pop0, nInd = nfam_M, use = "rand", sex = "M")
          )
		  dfInb[g-(6+Yd)] <- mateRes[['Inb']]
        } else {
          famselect <- c(
            selectWithinFam(pop = pop0, nInd = 1, use = "rand", sex = "F"),
            selectWithinFam(pop = pop0, nInd = 1, use = "rand", sex = "M")
          )
        }

		Ped <- rbind(Ped, data.frame(getPed(famselect), G = g))
        Phen <- rbind(Phen, data.frame(
          id = famselect@id,
          sex = famselect@sex,
          G = g,
          pheno(famselect)
        ))

        a <- Phen[Phen$id %in% parent@id, ]
		parent_org <- parent
        parent <- pop_select
		if(Yd > 1){
		for(i in 2:Yd){
			b <- a[a$sex == "F" & (a$G == g - i), ]
			b <- sample(b$id, Pars$cull_d[i-1]*length(b$id))
			parent <- c(parent, parent_org[b]) 
		}}
		if(Ys > 1){
		for(i in 2:Ys){
			b <- a[a$sex == "M" & (a$G == g - i), ]
			b <- sample(b$id, Pars$cull_s[i-1]*length(b$id))
			parent <- c(parent, parent_org[b]) 
		}}
		rm(b)

        if (Method != "blup") {
          # 基因型文件写出#
		  if (file.exists(paste0(rdmu_path, "/", g + 1)) == FALSE) {
            dir.create(paste0(rdmu_path, "/", g + 1))
          }
		  geno <- c(
            selectWithinFam(famselect, nInd = nGeno_F, sex = "F", use = "rand"),
            selectWithinFam(famselect, nInd = nGeno_M, sex = "M", use = "rand")
          )
          ID <- geno@id
          Geno <- pullSnpGeno(geno)
          if (g >= 7) {
            Geno <- rbind(
              cbind(as.data.frame(matrix(0, nrow = length(ID), ncol = 1)),
                id = ID, as.data.frame(matrix(0, nrow = length(ID), ncol = 4)), Geno
              ),
              Geno2
            )
          } else if (g == 6) {
            Geno <- cbind(as.data.frame(matrix(0, nrow = length(ID), ncol = 1)),
              id = ID, as.data.frame(matrix(0, nrow = length(ID), ncol = 4)), Geno
            )
            Geno <- Geno[Geno$id %in% Phen[Phen$G >= g - 10, "id"], ] ## 近10代基因型
            Geno2 <- Geno
          }
        }
      }
	  if (Method == "ssgblup") {
	    g_matrix_cal(snp_file = 'Geno',
                     inverse = FALSE,
					 prefix_name = paste0(rdmu_path, "/", g + 1))
	  }else if (Method == "gblup") {
	    g_matrix_cal(snp_file = 'Geno',
                     inverse = TRUE,
					 prefix_name = paste0(rdmu_path, "/", g + 1))
	  }
      P_G2[g - (6+Yd), ] <- c(
        meanP(pop0), meanG(pop0),
        meanP(pop0[pop0@sex == "M"]), meanG(pop0[pop0@sex == "M"]),
        meanP(pop0[pop0@sex == "F"]), meanG(pop0[pop0@sex == "F"]))
      VARG[g - (6+Yd), ] <- c(diag(varG(pop0)))

  g_end <- 9 + Yd + Time
  for(g in (8 + Yd):g_end){
  
    cat(paste0("run generation", g, "\n"))
	if (Method != "blup") {
	  if (file.exists(paste0(rdmu_path, '/', g+1)) == FALSE) {
	    dir.create(paste0(rdmu_path, '/', g+1))
	  }
	} else {
	  if (file.exists(paste0(rdmu_path, '/', g)) == FALSE) {
	    dir.create(paste0(rdmu_path, "/", g))
	  }
	}
    #写出系谱和表型文件到指定文件夹
	Grow = Phen[Phen$id %in% c(Ped$id, Ped$father, Ped$mother), c('id', 'G', TRAITS)]; Grow$b = 1
    write.table(Grow[(Grow$G <= g-1) & (Grow$G >= g-11), c('id', 'b', TRAITS)], paste0(rdmu_path, '/', g, "/Grow.txt"), 
                quote = F, row.names = F, col.names = F, sep = '\t')
    write.table(cbind(Ped[(Ped$G <= g-1) & (Ped$G >= g-11) , c(1, 3, 2, 4)]), paste0(rdmu_path, '/', g, "/ped.txt"), 
                quote = F, row.names = F, col.names = F, sep = '\t')
	
	# 估计育种值并依据选择指数选种
    source(system.file("extra_code","2.12_core.R", package = "GOplan"))
    
	pop_select = famselect[as.character(c(Select_F, Select_M))]
	ebv_all <- rbind(ebv_all,
			data.frame(Id = pop_select@id, Index = ebv[ebv$Id %in% pop_select@id, 'Index']))
	
	P_G[g-(7+Yd), ] = c(meanP(pop_select), meanG(pop_select), 
                       meanP(pop_select[pop_select@sex=='M']), meanG(pop_select[pop_select@sex=='M']), 
                       meanP(pop_select[pop_select@sex=='F']), meanG(pop_select[pop_select@sex=='F']))
					   
	if (Mate == 'rand') {
	  mateRes = Mate_rand()
	  dfInb[g-(6+Yd)] <- mateRes[['Inb']]      #近交系数=亲缘系数/2

	} else if (Mate == 'MC') {
	  Mate_file()
	  runMC(Nplans, getwd())
	  crossplan <- read.table('reMate.txt', header = T)
	  dfInb[g-(6+Yd)] <- crossplan[nrow(crossplan), 'id_s']
	  crossplan <- as.matrix(crossplan[-nrow(crossplan), ])
	} else {    ##调用其他软件
	  Mate_file()
	  source(Mate)
	  crossplan <- as.matrix(read.table('reMate.txt', header = T))
	}
	
    pop0 <- makeCross(parent, crossPlan = mateRes[['crossplan']], nProgeny = n_progeny)
	pop0 = setPheno(pop0, h2=heri)
    
	P_G2[g-(6+Yd), ] = c(meanP(pop0), meanG(pop0), 
                  meanP(pop0[pop0@sex=='M']), meanG(pop0[pop0@sex=='M']), 
                  meanP(pop0[pop0@sex=='F']), meanG(pop0[pop0@sex=='F']))
	VARG[g-(6+Yd), ] = c(diag(varG(pop0)))
	if(nfam_F == 0){
	  famselect = pop0   #全测
	}else{
	  famselect = c(selectWithinFam(pop = pop0, nInd = nfam_F, use = "rand", sex = "F"),
	                selectWithinFam(pop = pop0, nInd = nfam_M, use = "rand", sex = "M"))
	}			
	
	Ped = rbind(Ped, data.frame(getPed(famselect), G = g))
	Phen = rbind(Phen, data.frame(id = famselect@id, 
                                  sex = famselect@sex,
					              G = g,
                                  pheno(famselect)))
    
	a = Phen[Phen$id %in% parent@id, ]            #提取parent的性别和G信息
	parent_org <- parent
	parent <- pop_select
	if(Yd > 1){
		for(i in 2:Yd){
			b <- a[a$sex == "F" & (a$G == g - i), ]           #get remained generation
			b <- sample(b$id, Pars$cull_d[i-1]*length(b$id))  #cull parent
			parent <- c(parent, parent_org[b]) 
	}}
	if(Ys > 1){
		for(i in 2:Ys){
			b <- a[a$sex == "M" & (a$G == g - i), ]
			b <- sample(b$id, Pars$cull_s[i-1]*length(b$id))
			parent <- c(parent, parent_org[b]) 
	}}
		
    if (Method != "blup") {
      geno = c(selectWithinFam(famselect, nInd = nGeno_F, sex = "F", use = "rand"),
               selectWithinFam(famselect, nInd = nGeno_M, sex = "M", use = "rand"))	
      ID = geno@id	
      Geno = pullSnpGeno(geno)
      Geno = rbind(cbind(as.data.frame(matrix(0,nrow=length(ID),ncol=1)),
                   id = ID, as.data.frame(matrix(0,nrow=length(ID),ncol=4)),Geno),
        	       Geno2)
      Geno = Geno[Geno$id %in% Phen[Phen$G >= g-10, "id"], ]    ##近10代基因型
      Geno2 = Geno
	  if(Method == 'ssgblup'){
		g_matrix_cal(snp_file = 'Geno',
						inverse = FALSE,
						prefix_name = paste0(rdmu_path, "/", g + 1))
	  }else if (Method == "gblup") {
	    g_matrix_cal(snp_file = 'Geno',
                     inverse = TRUE,
					 prefix_name = paste0(rdmu_path, "/", g + 1))
	  }
	}	
  }

  for(i in 1:nrow(P_G2)){
    EBV[i, ] = 0
    for(j in 1:length(TRAITS)){
      EBV[i, ] <- ((P_G2[i, j] - P_G2[1, j])*w[j]/P_G2[1, j]) + EBV[i, ]
    }
  }
  out_ls1 = list()
  out_ls1[[paste0('P_G-select-', var_name)]] = P_G
  out_ls1[[paste0('P_G-pop0-', var_name)]] = P_G2
  out_ls1[[paste0('P_G-my-', var_name)]] = P_G3
  out_ls1[[paste0('VARG-pop0-', var_name)]] = VARG
  out_ls1[[paste0('EBV-', var_name)]] = EBV
  out_ls1[[paste0('acc-', var_name)]] = acc
  out_ls1[[paste0('acc-index-', var_name)]] = acc2
  out_ls1[[paste0("dfInb-", var_name)]] <- dfInb
  assign(paste0('out_', x), out_ls1)