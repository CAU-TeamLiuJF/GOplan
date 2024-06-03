#' calculate gene flow when only A has secondary group
#'
#' @param par1 cross type of whole system: 11,12,13,21,22,23,24,31,32,33
#' @param par2 Time of whole selection years
#'
#' @return
#' @export
#'
#' @examples
fun_gap = function(par1 = Ctype, par2 = Time){
  ls_b = list()
  #----------------------------------------------------
  nrow_p = (2*sire0 - 2 + As + Bs) +
    (2*dam0 - 2 + Ad + Bd) + 1
  
  AM = 1; AF = sire0+As
  BM = AF+dam0+Ad-1; BF = BM+sire0+Bs-1
  #MB = DF+dam0+Dd-1
  #MA = MB+sire0+Bsec-1
  #MD = MA+dam0+Asec-1
  #MC = MD+sire0+Dsec-1
  #CD = MC+dam0+Csec-1
  #AB = CD+sire0+YCDs-1
  #ABCD = AB+dam0+YABd-1
  
  cAM = sire0; cAF = cAM+As+dam0-1
  cBM = cAF+Ad+sire0-1; cBF = cBM+Bs+dam0-1 
  #cMB = cDF+Dd+sire0-1; cMA = cMB+Bsec+dam0-1
  #cMD = cMA+Asec+sire0-1; cMC = cMD+Dsec+dam0-1
  #cCD = cMC+Csec+sire0-1; cAB = cCD+YCDs+dam0-1

  gAM = 0.5*Aphe$sPP; gAF = 0.5*Aphe$dPP
  gBM = 0.5*Bphe$sPP; gBF = 0.5*Bphe$dPP
  gMA = 0.5*Aphe$secPP; gMB = 0.5*Bphe$secPP
  
  if(substr(par1, 1,1) == 2){
  if(substr(par1, 2,2) == 1){
	assign(Bsec, -1)
	nrow_p = nrow_p + dam0 - 1 + Asec
	MA = BF+dam0+Bd-1; PM = MA+dam0+Asec-1; cMA = cBF+Bd+dam0-1
	rgroup = c('AM','AM','AF','AF','BM','BM', 'BF', 'BF','MA','MA', 'PM','PM')
	cggroup = c('AM','AF','AM','AF','BM', 'BF','BM', 'BF','AM','AF','BM','MA')
	fei1 = c(cAM+As-1, cAF+Ad-1, cBM+Bs-1, cBF+Bd-1, cMA+Asec-1)
  } else if(substr(par1, 2,2) == 2){ 
	nrow_p = nrow_p + dam0 + sire0 - 2 + Asec + Bsec
	MA = BF+dam0+Bd-1; MB = MA+dam0+Asec-1; PM = MB+sire0+Bsec-1
	cMA = cBF+Bd+dam0-1; cMB = cMA+Asec+sire0-1
	rgroup = c('AM','AM','AF','AF','BM','BM', 'BF', 'BF','MA','MA','MB','MB','PM','PM')
	cggroup = c('AM','AF','AM','AF','BM', 'BF','BM', 'BF','AM','AF','BM', 'BF','MA','MB')
	fei1 = c(cAM+As-1, cAF+Ad-1, cBM+Bs-1, cBF+Bd-1, cMA+Asec-1, cMB+Bsec-1)
  
  }}else if(substr(par1, 1,1) == 3){
  CM = BF+dam0+Bd-1; CF = CM+sire0+Cs-1
  cCM = cBF+Bd+sire0-1; cCF = cCM+Cs+dam0-1
  gCM = 0.5*Cphe$sPP; gCF = 0.5*Cphe$dPP; gMC = 0.5*Cphe$secPP; gAB = 0.5*abPP;
  if(substr(par1, 2,2) == 1){ 
	assign('Bsec', -1, envir = .GlobalEnv); assign('Csec', -1, envir = .GlobalEnv)
	nrow_p = nrow_p + 3*dam0 + sire0 - 4 + YABd + Asec + Cd + Cs
	MA = CF+dam0+Cd-1; AB = MA+dam0+Asec-1; PM = AB+dam0+YABd-1
	cMA = cCF+Cd+dam0-1; cAB=cMA+Asec+dam0-1
	rgroup = c('AM','AM','AF','AF','BM','BM', 'BF','BF', 'CM','CM','CF','CF','MA','MA','AB','AB','PM','PM')
	cggroup = c('AM','AF','AM','AF','BM', 'BF','BM', 'BF','CM','CF','CM','CF','AM','AF','BM','MA','CM','AB')
	fei1 = c(cAM+As-1, cAF+Ad-1, cBM+Bs-1, cBF+Bd-1, cCM+Cs-1, cCF+Cd-1, cMA+Asec-1, cAB+YABd-1)
  }else if(substr(par1, 2,2) == 2){ 
	assign('Bsec', -1, envir = .GlobalEnv)
	nrow_p = nrow_p + 3*dam0 + 2*sire0 - 5 + YABd + Asec + Csec + Cd + Cs
	MA = CF+dam0+Cd-1; MC = MA+dam0+Asec-1;AB = MC+sire0+Csec-1; PM = AB+dam0+YABd-1
	cMA = cCF+Cd+dam0-1; cMC = cMA+Asec+sire0-1; cAB=cMC+Csec+dam0-1
	rgroup = c('AM','AM','AF','AF','BM','BM', 'BF','BF', 'CM','CM','CF','CF',
			   'MA','MA','MC','MC','AB','AB','PM','PM')
	cggroup = c('AM','AF','AM','AF','BM', 'BF','BM', 'BF','CM','CF','CM','CF',
			   'AM','AF','CM','CF','BM','MA','MC','AB')
	fei1 = c(cAM+As-1, cAF+Ad-1, cBM+Bs-1, cBF+Bd-1, cCM+Cs-1, cCF+Cd-1, cMA+Asec-1, cMC+Csec-1, cAB+YABd-1)
  }else if(substr(par1, 2,2) == 3){ 
	nrow_p = nrow_p + 3*dam0 + 3*sire0 - 6 + YABd + Asec + Bsec + Csec + Cd + Cs
	MA = CF+dam0+Cd-1; MB = MA+dam0+Asec-1;MC = MB+sire0+Bsec-1;AB = MC+sire0+Csec-1; PM = AB+dam0+YABd-1
	cMA = cCF+Cd+dam0-1; cMB = cMA+Asec+sire0-1; cMC = cMB+Bsec+sire0-1; cAB=cMC+Csec+dam0-1
	rgroup = c('AM','AM','AF','AF','BM','BM', 'BF','BF', 'CM','CM','CF','CF',
			   'MA','MA','MB','MB','MC','MC','AB','AB','PM','PM')
	cggroup = c('AM','AF','AM','AF','BM', 'BF','BM', 'BF','CM','CF','CM','CF',
			   'AM','AF','BM', 'BF','CM','CF','MA','MB','MC','AB')
	fei1 = c(cAM+As-1, cAF+Ad-1, cBM+Bs-1, cBF+Bd-1, cCM+Cs-1, cCF+Cd-1, cMA+Asec-1, cMB+Bsec-1, cMC+Csec-1, cAB+YABd-1)
  
  }}else if(substr(par1, 1,1) == 4){
    nrow_p = nrow_p + 3*dam0 + 3*sire0 - 3 + YABd + YCDs + Cd + Cs + Dd + Ds
	CM = BF+dam0+Bd-1; CF = CM+sire0+Cs-1
    cCM = cBF+Bd+sire0-1; cCF = cCM+Cs+dam0-1
    gCM = 0.5*Cphe$sPP; gCF = 0.5*Cphe$dPP
	DM = CF+dam0+Cd-1; DF = DM+sire0+Ds-1
	cDM = cCF+Cd+sire0-1; cDF = cDM+Ds+dam0-1
	gDM = 0.5*Dphe$sPP; gDF = 0.5*Dphe$dPP
	gMC = 0.5*Cphe$secPP; gMD = 0.5*Dphe$secPP; gCD = 0.5*cdPP
  if(substr(par1, 2,2) == 1){ 
	assign('Bsec', -1, envir = .GlobalEnv); assign('Dsec', -1, envir = .GlobalEnv)
	nrow_p = nrow_p + dam0 - 1 + Asec
	MA = DF+dam0+Dd-1; AB = MA+dam0+Asec-1; CD = AB+dam0+YABd-1;PM = CD+sire0+YCDs-1
	cMA = cDF+Dd+dam0-1; cAB = cMA+Asec+dam0-1; cCD = cAB+YABd+sire0-1
	rgroup = c('AM','AM','AF','AF','BM','BM', 'BF','BF', 'CM','CM','CF','CF', 'DM','DM','DF','DF',
			   'MA','MA','AB','AB','CD','CD','PM','PM')
	cggroup = c('AM','AF','AM','AF','BM', 'BF','BM', 'BF','CM','CF','CM','CF','DM','DF','DM','DF',
			    'AM','AF','BM','MA','CF','DM','AB','CD')
	fei1 = c(cAM+As-1, cAF+Ad-1, cBM+Bs-1, cBF+Bd-1, cCM+Cs-1, cCF+Cd-1, cDM+Ds-1, cDF+Dd-1, 
			 cMA+Asec-1, cAB+YABd-1, cCD+YCDs-1)
  }else if(substr(par1, 2,2) == 2){ 
	nrow_p = nrow_p + 2*dam0 + 2*sire0 - 4 + Asec + Bsec + Csec + Dsec
	MA = DF+dam0+Dd-1; MB = MA+dam0+Asec-1; MC = MB+sire0+Bsec-1; MD = MC+dam0+Csec-1 
	AB = MD+sire0+Dsec-1; CD = AB+dam0+YABd-1;PM = CD+sire0+YCDs-1
	cMA = cDF+Dd+dam0-1; cMB = cMA+Asec+sire0-1; cMC = cMB+Bsec+dam0-1; cMD = cMC+Csec+sire0-1
	cAB = cMA+Dsec+dam0-1; cCD = cAB+YABd+sire0-1
	rgroup = c('AM','AM','AF','AF','BM','BM', 'BF','BF', 'CM','CM','CF','CF', 'DM','DM','DF','DF',
			   'MA','MA','MB','MB','MC','MC','MD','MD','AB','AB','CD','CD','PM','PM')
	cggroup = c('AM','AF','AM','AF','BM', 'BF','BM', 'BF','CM','CF','CM','CF','DM','DF','DM','DF',
			    'AM','AF','BM','BF','CM','CF','DM','DF','MA','MB','MC','MD','AB','CD')
	fei1 = c(cAM+As-1, cAF+Ad-1, cBM+Bs-1, cBF+Bd-1, cCM+Cs-1, cCF+Cd-1, cDM+Ds-1, cDF+Dd-1, 
			 cMA+Asec-1, cMB+Bsec-1, cMC+Csec-1, cMD+Dsec-1, cAB+YABd-1, cCD+YCDs-1)
  }}
  
	q = matrix(0, nrow_p, nrow_p)
	for(i in 1:nrow_p-1){
		if((i %in% fei1) == FALSE){
		  q[i+1, i] = 1
		}
	}
	p = q
	for (i in 1:length(rgroup)) {
    R = data.frame(matrix(0, nrow_p, nrow_p))
    
    r = get(rgroup[i])
    c = get(paste0('c',cggroup[i]))
    g = get(paste0('g',cggroup[i]))
    R[r, c:(c+length(g)-1)] = g
    assign(paste0('R', i), R)
	p = p + R
    }
	names(p) = c(1:nrow(p))
	write.csv(p, paste0(out_path, 'p.csv'), row.names = F, quote=F)
	#----------------------------------------------------
	h2 = c(rep(0,nrow_p-1),1) 
	rrow = c()
	for (j in 1:length(rgroup)) {
		n1 = c(rep(0,nrow_p))
		for(i in 1:nrow_p){
		  if((i %in% seq(1, length(rgroup), 2)) == TRUE){
			n1[get(rgroup[i])] = 1
		  }
		}
		m1 = m0 = c(rep(0,nrow_p))
		b=0;G=0;R=0;m_sum=0
		
		for (i in 2:par2) {
		  R = get(paste0('R', j))
		  m = ma_vec(data.matrix(R), n1) + ma_vec(data.matrix(p), m1)
		  n = ma_vec(data.matrix(q), n1)
		  m1 = m; n1 = n
		  b = cbind(b, m1)           
		  m_sum = m_sum + m1        
		}              
		#ls_b[[paste0('b', j)]] = b
		assign(paste0('b', j), b)
		#assign(paste0('m', j, '_sum'), m_sum)
		#assign(paste0('G', j), G)
	    rrow = c(rrow,get(rgroup[j]))
	}
	
	ls_b[['rrow']] = rrow
	ls_b[['rrow_P']] = rrow[length(rrow)]
    if(exists('gAB')){ls_b[['gAB']] = gAB; ls_b[['rrow_AB']] = AB}
	if(exists('gCD')){ls_b[['gCD']] = gCD; ls_b[['rrow_CD']] = CD}
	##extract the gene attribution vector of selection group (if not exists: 0)
	dt = data.frame(matrix(0, ncol = ncol(b1), nrow = nrow(b1)))
	for(i in c('A', 'B', 'C', 'D')){
		if(exists(paste0('M', i))){
			assign(paste0('M', i, 'S'), get(paste0('b', which(cggroup == paste0('M', i)))), envir = .GlobalEnv)
		}else{assign(paste0('M', i, 'S'), get(paste0('b', tail(which(cggroup == paste0(i, 'M')), 1))), envir = .GlobalEnv)}
		if((i %in% c('A', 'B')) | (i == 'C' & substr(par1, 1,1) == 3) | (i == 'D' & substr(par1, 1,1) == 4)){
			assign(paste0(i,'S'), get(paste0('b', which(cggroup == paste0(i,'M'))[1]))+get(paste0('b', which(cggroup == paste0(i,'M'))[2])), envir = .GlobalEnv)
			assign(paste0(i,'D'), get(paste0('b', which(cggroup == paste0(i,'F'))[1]))+get(paste0('b', which(cggroup == paste0(i,'F'))[2])), envir = .GlobalEnv)
		}else{assign(paste0(i,'S'), dt, envir = .GlobalEnv);assign(paste0(i,'D'), dt, envir = .GlobalEnv)}
    }

	return(ls_b)
}

