if(QUICK == TRUE & g > 8+Yd ){
	ebv <- SimAcc()
	gv <- data.frame(Id = famselect@id, gv(famselect))
	gv$g_Index <- 0
}else{
	df = matrix(c(ana_phe, rep(g, Nclus), rep(rdmu_path, Nclus), rep(Method, Nclus)), 
				nrow = Nclus, ncol = 4)
	for(i in 1:nrow(df)){
	  if (file.exists(paste0(df[i, 3], "/", df[i, 2], "/", df[i, 1])) == FALSE) {
		  dir.create(paste0(df[i, 3], "/", df[i, 2], "/", df[i, 1]))
	  }
	  FUN_dir(df[i, 1], paste0(df[i, 3], "/", df[i, 2], "/", df[i, 1]), df[i, 4])
	}
	for(i in 1:nrow(df)){
		custom.function(c(t(df[i, ])))
	}
	#get ebv result of single trait analyse
	for(i in TRAITS[which(analyse == 1)]){
	  if(file.exists(paste0(rdmu_path, '/', g, "/", i, "/SOL"))){
		file.rename(paste0(rdmu_path, '/', g, "/", i, "/SOL"), paste0(rdmu_path, '/', g, "/", i, "/dmu.SOL"))
	  }
	  a = read.table(paste0(rdmu_path, '/', g, "/", i, "/dmu.SOL"), skip = 2)[, c(5, 8)]
	  names(a) = c('Id', paste0(i, '_EBV'))
	  assign(paste0(i, '_ebv'), a)
	}
	#get ebv result of multiple trait analyse
	double_phe = ana_phe[ana_phe %in% TRAITS[which(analyse == 1)] == FALSE]
	b <- as.data.frame(matrix(NA, ncol = 2, nrow = length(double_phe)))
	rownames(b) = double_phe
	for(i in 1:length(double_phe)){
	  b[i, ] = unlist(strsplit(double_phe[i], '_'))
	}
	for(i in TRAITS[which(analyse != 1)]){
	  i_ebvfile = rownames(which(b == i, arr.ind = TRUE))
	  if(file.exists(paste0(rdmu_path, '/', g, "/", i, "/SOL"))){
		file.rename(paste0(rdmu_path, '/', g, "/", i, "/SOL"), paste0(rdmu_path, '/', g, "/", i, "/dmu.SOL"))
	  }
	  a = read.table(paste0(rdmu_path, '/', g, "/", i_ebvfile, "/dmu.SOL"), skip = 3)[, c(3, 5, 8)]
	  a = a[a$V3 == which(b == i, arr.ind = TRUE)[2], 2:3]
	  names(a) = c('Id', paste0(i, '_EBV'))
	  assign(paste0(i, '_ebv'), a)
	  rm(a)
	}

	# calculate the index and select
	ebv1 = list()
	for(i in TRAITS){
	  ebv1[[i]] = get(paste0(i, '_ebv'))
	}
	ebv1 <- reduce(ebv1, inner_join, by = 'Id')

	gv <- data.frame(Id = famselect@id, gv(famselect))
	ebv1$Index <- gv$g_Index <- 0
	for(i in 1:length(TRAITS)){
	  ebv1$Index <- as.vector(scale(ebv1[, paste0(TRAITS[i], '_EBV')]) * w[i] + ebv1$Index)
	  gv$g_Index <- as.vector(scale(gv[, TRAITS[i]]) * w[i] + gv$g_Index)
	}
	ebv <- ebv1[ebv1$Id %in% famselect@id, ]
	ebv_gv <- merge(gv, ebv, by = "Id", sort = F)

	for(i in 1:length(TRAITS)){
	  acc[g - (7 + Yd), i] <- cor(ebv_gv[, TRAITS[i]], ebv_gv[, paste0(TRAITS[i], '_EBV')])
	}
	acc2[g - (7 + Yd), ] <- c(cor(ebv_gv$Index, ebv_gv$g_Index))
}
ebv <- ebv[order(ebv$Index, decreasing = T), ]
Select_F <- ebv[ebv$Id %in% Phen[Phen$sex == "F", "id"], ][1:nSelected_F, 'Id']
Select_M <- ebv[ebv$Id %in% Phen[Phen$sex == "M", "id"], ][1:nSelected_M, 'Id']

if(task != 1){
if (pz == "A" | (pz == "C" & substr(Ctype, 1,1) == 4)) {
  ebv2 <- ebv[ebv$Id %in% Phen[Phen$sex == "F", "id"], ]
  nSel2 <- nSelected_MF
  nSel1 <- nSelected_F
}else{
  ebv2 <- ebv[ebv$Id %in% Phen[Phen$sex == "M", "id"], ]
  nSel2 <- nSelected_M2
  nSel1 <- nSelected_M
}
if(nSel2 != 0){
	if (length(ebv2$Id) >= (nSel1 + nSel2)) {
	  pop_select2 <- famselect[as.character(ebv2[(nSel1 + 1):(nSel1 + nSel2), 'Id'])]
	  cat("select", pz, "2% = ", nSel2/(length(ebv2$Id)-nSel1), '\n')
	} else {
	  pop_select2 <- famselect[as.character(ebv2[(nSel1 + 1):length(ebv2$Id), 'Id'])]
	  cat("select", pz, "2% = 0%", '\n')
	}
	P_G3[g - (7 + Yd), ] <- c(meanP(pop_select2), meanG(pop_select2))
}
}
