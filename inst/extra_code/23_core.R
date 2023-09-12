###########################################
#本代码用于模拟核心群选育过程
#attention：建议取结果时从2开始
###########################################
##并行跑rep
fun_core = function(x){
	x <<- x
	#e <- new.env(parent = globalenv())
	#with(e, x <- x)
#env = attach(NULL, name=paste0('env',x))
	source(system.file("extra_code","1-1.R", package = "GOplan"))
}
cl <- makeCluster(as.numeric(Ncores))
clusterEvalQ(cl, expr = library(GOplan))
NCs2 = NCs20 = NCd2 = NCd20 = NDs2 = NDs20 = 0
NBs20 = NBs2 = NA2d = NA2d0 = NCDs0 = NCDs = NCDd = 0
NDd = NCd = NDs = NCs = NDd0 = NCd0 = NDs0 = NCs0 = 0
for(pz in Breed){
  load(paste0(out_path, pz, 'founder.Rdata'))
  	if(file.exists(paste0(out_path, pz)) == FALSE){
		dir.create(paste0(out_path, pz)) 
	}
	
	prm <- get(paste0(pz, 'prm'))
  	for(i in 1:length(prm)){assign(names(prm)[i], prm[i])}
	for(i in names(e)){assign(i, e[[i]])}
	
	##calculate the number of single trait analyse & mutile trait analyse
	analyse <- Pars$analyse
	w <- Pars$weigh
	TRAITS <- Pars$trait
	assign(paste0(pz, 'Trait'), TRAITS)
	
	a = names(table(analyse))
	a <- a[a != 1]
	ana_phe <- TRAITS[which(analyse == 1)]
	Nclus = length(analyse[analyse==1]) + length(a)
	if (length(a) != 0) {
		for(i in 1:length(a)){
			ana_phe <- append(ana_phe, paste(TRAITS[which(analyse == a[i])], collapse = "_"), after = length(ana_phe))
		}
	}
  rm(a)
  ### read variables combination nvar
  for(name in names(VAR)){assign(name, VAR[nvar, name])}
  for(i in c('Ys', 'Yd', 'Ysec')){assign(gsub('Y', pz, i), get(i))}
  ##if variables contain service time，need to recalculate prob----------------------------------------------------#
  if(length(intersect(names(VAR),  c('Ys', 'Yd'))) != 0){
    a =  get(paste0(pz, 'phe'))
	if('Ys' %in% names(VAR)){
	  a$sProb = 1
	  if(Ys > 1){for(i in 1:(Ys-1)){ a$sProb = a$sProb + prod(rep(a$cull_s[1], i)) }}
	  a$sPP = c(1/a$sProb)
	  if(Ys > 1){
		for(i in 1:(Ys-1)){ a$sPP = c(a$sPP, prod(rep(a$cull_s[1], i))/a$sProb) }
		a$cull_s = rep(a$cull_s[1], Ys-1)}
    }
	if('Yd' %in% names(VAR)){
      a$dProb = 1
	  if(Yd > 1){for(i in 1:(Yd-1)){ a$dProb = a$dProb + prod(rep(a$cull_d[1], i)) }}
	  a$dPP = c(1/a$dProb)
	  if(Yd > 1){for(i in 1:(Yd-1)){ a$dPP = c(a$dPP, prod(rep(a$cull_d[1], i))/a$dProb) }
	  a$cull_d = rep(a$cull_d[1], Yd-1)}
    }
	assign(paste0(pz, 'phe'), a)
	rm(a)
  }
  if(paste0(pz,  'sec') %in% names(VAR)){
	a =  get(paste0(pz, 'phe'))
    a$secProb = 1
	for(i in 1:(get(paste0(pz,  'sec'))-1)){ a$secProb = a$secProb + prod(rep(a$cull_sec[1], i)) }
	a$secPP = c(1/a$secProb)
	for(i in 1:(get(paste0(pz,  'sec'))-1)){ a$secPP = c(a$secPP, prod(rep(a$cull_sec[1], i))/a$secProb) }
	assign(paste0(pz, 'phe'), a)
	rm(a)
  }
  if('YABd' %in% names(VAR)){
    ABProb <- 1
	for(i in 1:(YABd-1)){ ABProb = ABProb + prod(rep(cull_AB[1], i)) }
	abPP <- c(1/ABProb)
	for(i in 1:(YABd-1)){ abPP = c(abPP, prod(rep(cull_AB[1], i))/ABProb) }
  }
  if('YCDs' %in% names(VAR)){
    CDProb <- 1
	for(i in 1:(YCDs-1)){ CDProb = CDProb + prod(rep(cull_CD[1], i)) }
	abPP <- c(1/CDProb)
	for(i in 1:(YCDs-1)){ cdPP = c(cdPP, prod(rep(cull_CD[1], i))/CDProb) }
  }
	Pars = get(paste0(pz, 'phe'))

  ##calculate population structure: female/male number, select number...
  n_male = n_female/Sor; assign(paste0('N',pz,'s'), n_male)
  assign(paste0('N',pz,'d'), n_female)
  nSelected_M = round(n_male/Pars$sProb)
  nSelected_F = round(n_female/Pars$dProb)
  nSelected_MF = nSelected_M2 = 0
  if(pz == 'A') {
	NA2d = n_MY = nSec
	NA2d0 = nSelected_MF = round(NA2d/Pars$secProb)
	Lit_A2 = n_progeny
  }
  if(pz == 'B'){
	NBs2 = NA2d/SorAB
	NBs20 = nSelected_M2 = NBs2/Pars$secProb
	if(Ctype %in% c(33, 42) == FALSE){
		NBs0 = nSelected_M = nSelected_M + NA2d/(SorAB*Pars$sProb)
		NBs = round(nSelected_M*Pars$sProb)}
  }
  if(pz == 'C'){
	if(Ctype %in% c(41, 42)){
		Lit_C2 = n_progeny
		NCDd = NCDs = (NA2d*Lit_A2*0.5*SR)*YABd/SorP; NCDs0 = NCDs/cdProb
		NCd2 =  NCDs0/(Lit_C2*CDSR)
		NCd20 = nSelected_MF = round(NCd2/Pars$secProb)
	}else if(Ctype %in% c(32,33)){
		NCs2 = (NA2d*Lit_A2*0.5*SR)*YABd/SorP
		NCs20 = nSelected_M2 = NCs2/Pars$secProb
	}else if(Ctype %in% c(31)){
		NCs0 = nSelected_M = nSelected_M + (NA2d*Lit_A2*0.5*SR)*YABd/(SorP*Pars$sProb)
		NCs = round(nSelected_M*Pars$sProb)
	}}
  if(pz == 'D'){
	if(Ctype != 42){
		NDs0 = nSelected_M = nSelected_M + NCDs0/(Lit_C2*SR*SorCD*Pars$sProb)
		NDs = round(nSelected_M*Pars$sProb)
	}else{
		NDs2 = NCDs0/(Lit_C2*SR*SorCD)
		NDs20 = nSelected_M2 = NDs2/Pars$sProb		
  }}
  assign(paste0('N',pz,'d0'), nSelected_F)
  assign(paste0('N',pz,'s0'), nSelected_M) 

  ##if variable contain nfam/ngeno set nfam_F/M,nGeno_F/M same as it 
  if(exists("nfam") == TRUE){
    nfam_M = nfam_F = nfam
  }
  if(exists('ngeno')){
	nGeno_F <- nGeno_M <- ngeno
  }
  save.image(file = paste0(out_path, '1.RData'))
  #####nrep####
  ##calculate n repetion in parallel
  clusterExport(cl = cl, list('SP', 'out_path'))
  out_list = parLapplyLB(cl = cl, 1:Nrep, fun_core)
  for(i in 1:Nrep){
	out_list[[i]] = out_list[[i]]$value
  }
  assign(paste0(pz, "_out"), out_list)
  gc()
}
save.image(file = paste0(out_path, '2.RData'))
stopCluster(cl)
rm(cl)
	