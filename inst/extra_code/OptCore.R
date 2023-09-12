fun_core = function(x){
x <<- x
source(system.file("extra_code", "1-opt.R", package = "GOplan"))
}
optCode <- 0
##judge if the number of nuclei is suitable
for(i in names(e)){
	assign(i, e[[i]])
}
NCs2 = NCs20 = NCd2 = NCd20 = NDs2 = NDs20 = 0
NBs20 = NBs2 = NA2d = NA2d0 = NCDs0 = NCDs = NCDd = 0
NDs = NCs = NDd0 = NCd0 = NDs0 = NCs0 = 0

ppp <- list()
cull_AB <- cull_AB[1]
	ABProb <- 1
	for(i in 1:(YABd-1)){ ABProb = ABProb + pow(cull_AB, i) }
	abPP <- c(1/ABProb)
	for(i in 1:(YABd-1)){ abPP = c(abPP, pow(cull_AB, i)/ABProb) }
	ppp[['abPP']] = abPP
cull_CD <- cull_CD[1]
	CDProb <- 1
	for(i in 1:(YCDs-1)){ CDProb = CDProb + pow(cull_CD, i) }
	cdPP <- c(1/CDProb)
	for(i in 1:(YCDs-1)){ cdPP = c(cdPP, pow(cull_CD, i)/CDProb) }
	ppp[['cdPP']] = cdPP
for(m in Breed){
	b = get(paste0(m, 'phe'))
	b$secProb = 1              ##total proportion of all age classes
	for(i in 1:(get(paste0(m, 'sec'))-1)){
		b$secProb = b$secProb + 
			pow(b$cull_sec[1], i)}
	b$secPP = c(1/b$secProb)
		for(i in 1:(get(paste0(m, 'sec'))-1)){
			b$secPP = c(b$secPP, pow(b$cull_sec[1], i)/b$secProb)}
	assign(paste0(m, 'phe'), b)
	ppp[[paste0(m, 'phe')]] <- b
}
	
Lit_AB = 0.5*(Aprm['n_progeny'] + Bprm['n_progeny'])
NABd <- NABs <- NABd_A <- NABd_B <- NABd_C <- N_P/PSR/Lit_AB
NABd0 = NABd/ABProb;
NA2d <- NABd0/ABSR/(Aprm['n_progeny']*0.5)
NA2d0 = NA2d/Aphe$secProb;
if(substr(Ctype,1,1) == 3){
	NDd = NCd2 = NCd20 = 0
	NCs2 = NABd/SorP; NCs20 = NCs2/Cphe$secProb
	NBs2 = NA2d/SorAB; NBs20 = NBs2/Bphe$secProb
	if((NA2d0 + NAd/Aphe$dProb) > NAd*Aprm['nfam_F']){      #if NAd can not supply enough progeny
		  NABd_A <- (NAd*Aprm['nfam_F']-NAd/Aphe$dProb)*Aphe$secProb*   #NA2d
					ABSR*(Aprm['n_progeny']*0.5)*ABProb
	}
	if(Ctype == 33){
	if((NBs20 + NBd/Bprm['Sor']/Bphe$sProb) > NBd*Bprm['nfam_M']){
		  NABd_B <- (NBd*Bprm['nfam_M']-NBd/Bprm['Sor']/Bphe$sProb)*
					Bphe$secProb*SorAB*Aphe$secProb*ABSR*(Aprm['n_progeny']*0.5)*ABProb
	}
	if((NCs20 + NCd/Cprm['Sor']/Cphe$sProb) > NCd*Cprm['nfam_M']){
		  NABd_C <- (NCd*Cprm['nfam_M']-NCd/Cprm['Sor']/Cphe$sProb)*Cphe$secProb*SorP
	}}
	if(Ctype == 32){
	NBs20 = NBs2/Bphe$sProb
	if((NBs20 + NBd/Bprm['Sor']/Bphe$sProb) > NBd*Bprm['nfam_M']){
		  NABd_B <- (NBd*Bprm['nfam_M']-NBd/Bprm['Sor']/Bphe$sProb)*Bphe$sProb*SorAB* #NA2d
					Aphe$secProb*ABSR*(Aprm['n_progeny']*0.5)*ABProb
	}
	if((NCs20 + NCd/Cprm['Sor']/Cphe$sProb) > NCd*Cprm['nfam_M']){
		  NABd_C <- (NCd*Cprm['nfam_M']-NCd/Cprm['Sor']/Cphe$sProb)*Cphe$secProb*SorP
	}}
	if(Ctype == 31){
	NBs20 = NBs2/Bphe$sProb; NCs20 = NCs2/Cphe$sProb
	if((NBs20 + NBd/Bprm['Sor']/Bphe$sProb) > NBd*Bprm['nfam_M']){
		  NABd_B <- (NBd*Bprm['nfam_M']-NBd/Bprm['Sor']/Bphe$sProb)*Bphe$sProb*SorAB* #NA2d
					Aphe$secProb*ABSR*(Aprm['n_progeny']*0.5)*ABProb
	}
	if((NCs20 + NCd/Cprm['Sor']/Cphe$sProb) > NCd*Cprm['nfam_M']){
		  NABd_C <- (NCd*Cprm['nfam_M']-NCd/Cprm['Sor']/Cphe$sProb)*Cphe$sProb*SorP
	}}
	NABd <- NABs <- min(NABd_A, NABd_B, NABd_C)
	NABd0 = NABd/ABProb;
	NA2d <- NABd0/ABSR/(Aprm['n_progeny']*0.5)
	NA2d0 = NA2d/Aphe$secProb;
	NCs2 = NABd/SorP; 
	NBs2 = NA2d/SorAB; 
	if(Ctype == 32){
	NBs20 = NBs2/Bphe$sProb;
	NCs20 = NCs2/Cphe$sProb;
	}else if(Ctype == 33){
	NBs20 = NBs2/Bphe$sProb;
	NCs20 = NCs2/Cphe$secProb;
	}else if(Ctype == 34){
	NBs20 = NBs2/Bphe$secProb;
	NCs20 = NCs2/Cphe$secProb;
	}
	cat("NABd: ", NABd, "\n")
	cat("N_P: ", NABd*PSR*Lit_AB, "\n")
	if(NABd*PSR*Lit_AB < 0.9*N_P){
		cat("This combination isn't suitable!\n")
	}
}else if(substr(Ctype,1,1) == 4){
	NCs2 = NCs20 = 0
	NCDs = NABd/SorP; NCDs0 = NCDs/CDProb
	NCd2 = NCDs0/CDSR/(Cprm['n_progeny']*0.5); NCd20 = NCd2/Cphe$secProb
	NBs2 = NA2d/SorAB; NBs20 = NBs2/Bphe$secProb
	NDs2 = NCd2/SorCD; NDs20 = NDs2/Dphe$secProb
	if((NA2d0 + NAd/Aphe$dProb) > NAd*Aprm['nfam_F']){
		  NABd_A <- (NAd*Aprm['nfam_F']-NAd/Aphe$dProb)*
					Aphe$secProb*ABSR*(Aprm['n_progeny']*0.5)*ABProb
	}
	if((NCd20 + NCd/Cphe$dProb) > NCd*Cprm['nfam_F']){
		  NABd_C <- (NCd*Cprm['nfam_F']-NCd/Cphe$dProb)*Cphe$secProb*
					(Cprm['n_progeny']*0.5)*CDSR*CDProb*SorP
	}
	if(Ctype == 42){
	if((NBs20 + NBd/Bprm['Sor']/Bphe$sProb) > NBd*Bprm['nfam_M']){
		  NABd_B <- (NBd*Bprm['nfam_M']-NBd/Bprm['Sor']/Bphe$sProb)*
					Bphe$secProb*SorAB*Aphe$secProb*ABSR*(Aprm['n_progeny']*0.5)*ABProb
	}
	if((NDs20 + NDd/Dprm['Sor']/Dphe$sProb) > NDd*Dprm['nfam_M']){
		  NABd_D <- (NDd*Dprm['nfam_M']-NDd/Dprm['Sor']/Dphe$sProb)*Dphe$secProb*SorCD*
					(Cprm['n_progeny']*0.5)*CDSR*CDProb*SorP
	}}
	if(Ctype == 41){
	NBs20 = NBs2/Bphe$sProb
	NDs20 = NDs2/Dphe$sProb
	if((NBs20 + NBd/Bprm['Sor']/Bphe$sProb) > NBd*Bprm['nfam_M']){
		  NABd_B <- (NBd*Bprm['nfam_M']-NBd/Bprm['Sor']/Bphe$sProb)*
					Bphe$sProb*SorAB*Aphe$secProb*ABSR*(Aprm['n_progeny']*0.5)*ABProb
	}
	if((NDs20 + NDd/Dprm['Sor']/Dphe$sProb) > NDd*Dprm['nfam_M']){
		  NABd_D <- (NDd*Dprm['nfam_M']-NDd/Dprm['Sor']/Dphe$sProb)*Dphe$sProb*SorCD*
					(Cprm['n_progeny']*0.5)*CDSR*CDProb*SorP
	}}
	NABd <- NABs <- min(NABd_A, NABd_B, NABd_C, NABd_D)
	NABd0 = NABd/ABProb;
	NA2d <- NABd0/ABSR/(Aprm['n_progeny']*0.5)
	NA2d0 = NA2d/Aphe$secProb;
	NBs2 = NA2d/SorAB; 
	NCDs = NCDd = NABd/SorP; NCDs0 = NCDs/CDProb
	NCd2 = NCDs0/CDSR/(Cprm['n_progeny']*0.5);
	NCd20 = NCd2/Cphe$secProb;	
	NDs2 = NCd2/SorCD; 
	if(Ctype == 41){
	NBs20 = NBs2/Bphe$sProb;
	NDs20 = NDs2/Dphe$sProb;
	}else{
	NBs20 = NBs2/Bphe$secProb;
	NDs20 = NDs2/Dphe$secProb;
	}
	cat("NABd: ", NABd, "\n")
	cat("N_P: ", NABd*PSR*Lit_AB, "\n")
	if(NABd*PSR*Lit_AB < 0.9*N_P){
		cat("This combination isn't suitable!\n")
	}
}	
cl <- makeCluster(as.numeric(Ncores))
for(pz in Breed){
	load(paste0(out_path, pz, 'founder.Rdata'))	
	if(file.exists(paste0(out_path, pz)) == FALSE){
		dir.create(paste0(out_path, pz)) 
	}
	for(i in names(ppp)){
		assign(i, ppp[[i]])
	}
	
	prm <- get(paste0(pz, 'prm'))
  	for(i in 1:length(prm)){
		assign(names(prm)[i], prm[i])
	}
	for(i in c('Ys', 'Yd', 'Ysec')){assign(gsub('Y', pz, i), get(i))}
	for(i in names(e)){
		assign(i, e[[i]])
	}
	analyse <- Pars$analyse
	w <- Pars$weigh
	TRAITS <- Pars$trait
	assign(paste0(pz, 'Trait'), TRAITS)
	
	a = names(table(analyse))
	a <- a[a != 1]
	Nclus = length(analyse[analyse==1]) + length(a)
  
	ana_phe <- TRAITS[which(analyse == 1)]
	if (length(a) != 0) {
		for(i in 1:length(a)){
			ana_phe <- append(ana_phe, paste(TRAITS[which(analyse == a[i])], collapse = "_"), after = length(ana_phe))
		}
	}
	rm(a)

  out_list <- list()
  if(pz == 'B'){	
	n_female <- NBd; NBs <- n_male <- n_female/Sor
	nSelected_M2 <- NBs20
	NBs0 = nSelected_M = round(n_male/Pars$sProb); NBd0 = nSelected_F = round(n_female/Pars$dProb)
	if(Ctype %in% c(31,32,41)){
		nSelected_M2 <- NBs20 <- NBs2 <-  0
		NBs = round((n_female + NA2d)/SorAB)
		NBs0 = nSelected_M = round(NBs/Pars$sProb)
	}
	if(n_female > 1000){
		n_female = 1000
		n_male <- n_female/Sor
		nSelected_M2 <- NBs20*(n_female/NBd)
		nSelected_M = round(n_male/Pars$sProb); NBd0 = nSelected_F = round(n_female/Pars$dProb)
		if(Ctype %in% c(31,32,41)){
			n_male = round(NBs*n_female/NBd)  ##等比例缩小
			nSelected_M = round(n_male/Pars$sProb)
		}
	}
  }else if(pz == 'C'){
    n_female <- NCd; NCs <- n_male <- n_female/Sor
	nSelected_M2 <- NCs20
	NCs0 = nSelected_M = round(n_male/Pars$sProb); NCd0 = nSelected_F = round(n_female/Pars$dProb)
	nSelected_MF <- NCd20
	if(Ctype == 31){	
		NCs0 = nSelected_M = nSelected_M + round(NABd/(SorP*Pars$sProb))
		NCs = round(nSelected_M*Pars$sProb)
	}
	if(n_female > 1000){
		n_female = 1000
		n_male <- n_female/Sor
		nSelected_M2 <- NCs20*(n_female/NCd)
		nSelected_MF <- NCd20*(n_female/NCd)
		nSelected_M = round(n_male/Pars$sProb); nSelected_F = round(n_female/Pars$dProb)
		if(Ctype %in% c(31)){
			n_male = round(NCs*n_female/NCd)  ##等比例缩小
			nSelected_M = round(n_male/Pars$sProb)
		}
	}
  }else if(pz == 'D'){	
	n_female <- NDd; NDs <- n_male <- n_female/Sor
	nSelected_M2 <- NDs20
	NDs0 = nSelected_M = round(n_male/Pars$sProb); NDd0 = nSelected_F = round(n_female/Pars$dProb)
	if(Ctype %in% c(41)){
		nSelected_M2 <- NDs20 <- NDs2 <-  0
		NDs = round(n_female + NCDs/SorCD)
		NDs0 = nSelected_M = round(NDs/Pars$sProb)
	}
	if(n_female > 1000){
		n_female = 1000
		n_male <- n_female/Sor
		nSelected_M2 <- NDs20*(n_female/NDd)
		nSelected_M = round(n_male/Pars$sProb); nSelected_F = round(n_female/Pars$dProb)
		if(Ctype %in% c(41)){
			n_male = round(NDs*n_female/NDd)  ##等比例缩小
			nSelected_M = round(n_male/Pars$sProb)
		}
	}
  }else if(pz == 'A'){
	n_female <- NAd; NAs <- n_male <- n_female/Sor
    nSelected_MF <- NA2d0
	NAs0 = nSelected_M = round(n_male/Pars$sProb); NAd0 = nSelected_F = round(n_female/Pars$dProb)
	if(n_female > 1000){
		n_female = 1000
		n_male <- n_female/Sor
		nSelected_MF <- NA2d0*(n_female/NAd)
		nSelected_M = round(n_male/Pars$sProb); nSelected_F = round(n_female/Pars$dProb)
	}
  }

  ## prepare for file dataframes
  P_G = P_G2 = P_G3 = data.frame(matrix(ncol = length(TRAITS)*6, nrow = 0))
  names(P_G) = names(P_G2) = names(P_G3) = c(paste0('P_', TRAITS), paste0('G_', TRAITS),
			   paste0('P1_', TRAITS), paste0('G1_', TRAITS),
			   paste0('P0_', TRAITS), paste0('G0_', TRAITS))				 
  VARG2 = VARG = gv_all = data.frame(matrix(ncol = length(TRAITS), nrow = 0))
  names(VARG2) = names(VARG) = c(TRAITS)  
  EBV = data.frame(matrix(ncol = length(TRAITS)+1, nrow = 0))
  names(EBV) = c(paste0('g_', TRAITS), 'g_index')
  acc = data.frame(matrix(ncol = length(TRAITS), nrow = 0))
  names(acc) = c(TRAITS)
  acc2 = data.frame(matrix(ncol = 1, nrow = 0)); names(acc2) = c('index')
  ebv_all = data.frame(matrix(ncol = 2, nrow = 0))
  names(ebv_all) = c('Id', 'Index')
  #-------------------------------------------------------------------------------------------------
famselect <- c(
	selectWithinFam(pop = pop00, nInd = 1, use = "rand", sex = "F"),
	selectWithinFam(pop = pop00, nInd = 1, use = "rand", sex = "M")
)

parent <- c(
	selectInd(famselect, nInd = nSelected_F, sex = "F", use = "rand"),
	selectInd(famselect, nInd = nSelected_M, sex = "M", use = "rand")
)
Ped <- rbind(Ped00, data.frame(getPed(famselect), G = g))
Phen <- rbind(Phen00, data.frame(
	id = famselect@id,
	sex = famselect@sex,
	G = g,
	pheno(famselect)
))
famselect = famselect[famselect@id[which(famselect@id %in% parent@id == FALSE)]]
## construct the strcture of nuclei
for (g in 6:(7+Yd)) {
pop_select <- c(
	selectInd(famselect, nInd = nSelected_F, sex = "F", use = "rand"),
	selectInd(famselect, nInd = nSelected_M, sex = "M", use = "rand"))
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
	meanP(pop0[pop0@sex == "F"]), meanG(pop0[pop0@sex == "F"])
)
VARG[g - (6+Yd), ] <- c(diag(varG(pop0)))
famselect00 <- famselect 
Ped00 <- Ped; Phen00 <- Phen; parent00 <- parent; pop00 <- pop0

save.image(file = paste0(out_path, '1.RData'))##### Nrep####
Sys.sleep(3)
clusterExport(cl = cl, list('SP', 'out_path'), envir = .GlobalEnv)
  out_list = parLapplyLB(cl = cl, as.list(c(1:Nrep)), fun_core)
  for(i in 1:Nrep){
	out_list[[i]] = out_list[[i]]$value
  }
assign(paste0(pz, "_out"), out_list)
gc()
}
stopCluster(cl)
rm(cl)



