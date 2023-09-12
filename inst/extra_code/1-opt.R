library(AlphaSimR); library(tidyverse);library(stringr)
load(paste0(out_path, '1.RData'))
rdmu_path = paste0(out_path, pz, '/', x)
if(file.exists(rdmu_path) == FALSE){
	dir.create(rdmu_path) 
}
setwd(rdmu_path)
  
Ped <- Ped00; parent <- parent00
pop0 <- pop00
famselect <- famselect00
Phen <- Phen00 # 沿用之前的
	fams <- list()
	dfInb <- rep(NA, Time)
	fams[[g-(6+Yd)]] <- pop0@id
	
	g_end <- 9 + Yd + Time
      for (g in (8 + Yd):g_end) {
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
        # 写出系谱和表型文件到指定文件夹
        Grow <- Phen[Phen$id %in% c(Ped$id, Ped$father, Ped$mother), c("id", "G", TRAITS)]
        Grow$b <- 1
        write.table(Grow[(Grow$G <= g - 1) & (Grow$G >= g - 11), c("id", "b", TRAITS)], paste0(rdmu_path, "/", g, "/Grow.txt"),
          quote = F, row.names = F, col.names = F, sep = "\t"
        )
        write.table(cbind(Ped[(Ped$G <= g - 1) & (Ped$G >= g - 11), c(1, 3, 2, 4)]), paste0(rdmu_path, "/", g, "/ped.txt"),
          quote = F, row.names = F, col.names = F, sep = "\t"
        )
		#-------------------------------------------------------------------------------------------------------------
		# 估计育种值并依据选择指数选种
        source(system.file("extra_code","2.12_core.R", package = "GOplan"))
    #-------------------------------------------------------------------------------------------------------------
        pop_select <- famselect[as.character(c(Select_F, Select_M))]

        P_G[g - (7+Yd), ] <- c(
          meanP(pop_select), meanG(pop_select),
          meanP(pop_select[pop_select@sex == "M"]), meanG(pop_select[pop_select@sex == "M"]),
          meanP(pop_select[pop_select@sex == "F"]), meanG(pop_select[pop_select@sex == "F"])
        )

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
        pop0 <- setPheno(pop0, h2 = heri)

        P_G2[g - (6+Yd), ] <- c(
          meanP(pop0), meanG(pop0),
          meanP(pop0[pop0@sex == "M"]), meanG(pop0[pop0@sex == "M"]),
          meanP(pop0[pop0@sex == "F"]), meanG(pop0[pop0@sex == "F"])
        )
        VARG[g - (6+Yd), ] <- c(diag(varG(pop0)))
		fams[[g-(6+Yd)]] <- pop0@id
        if (nfam_F == 0) {
          famselect <- pop0 # 全测
        } else {
          famselect <- c(
            selectWithinFam(pop = pop0, nInd = nfam_F, use = "rand", sex = "F"),
            selectWithinFam(pop = pop0, nInd = nfam_M, use = "rand", sex = "M")
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
		for(i in 2:Yd){
			b <- a[a$sex == "F" & (a$G == g - i), ]
			b <- sample(b$id, Pars$cull_d[i-1]*length(b$id))
			parent <- c(parent, parent_org[b]) 
		}
		for(i in 2:Ys){
			b <- a[a$sex == "M" & (a$G == g - i), ]
			b <- sample(b$id, Pars$cull_s[i-1]*length(b$id))
			parent <- c(parent, parent_org[b]) 
		}
		rm(b)
		
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
  for(i in 1:length(P_G2)){
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
