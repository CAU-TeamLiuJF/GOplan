#' generate dir file for running dmu
#'
#' @param phe the phenotype name
#' @param path the out path of the dir file
#' @param Method the method of estimating breeding values ('blup' or 'ssgblup')
#'
#' @return
#' @export
#'
#' @examples
FUN_dir = function(phe, path, Method = Method){
  DIR = data.frame(NA)
  DIR[1, ] = '$COMMENT'
  DIR[c(2, 4, 6, 12, 14), ] = '   '
  DIR[3, ] = '$ANALYSE 1 1 0 0'
  DIR[5, ] = paste0('$DATA  ASCII (2,', length(TRAITS), ',-9999) ../Grow.txt')
  DIR[7, ] = '$VARIABLE'
  DIR[8, ] = '#I1 I2'
  DIR[9, ] = 'Id b'
  DIR[10, ] = '#R1 R2 R3 R4'
  DIR[11, ] = paste(TRAITS, collapse=" ")                         
  DIR[13, ] = '$VAR_STR 1 PED 2 ASCII   ../ped.txt'                          
  DIR[15, ] = '$MODEL'

  if(length(strsplit(phe, "_")[[1]]) == 1){
    DIR[16, ] = '1 1 0 0 0'
	DIR[c(17, 20, 21), ] = 0
	DIR[18, ] = paste0((which(names(Grow) == phe)-2), ' 0 2 2 1')
    DIR[19, ] = '1 1'
	DIR[c(22, 25, 26), ] = ' '
	DIR[23, ] = '$RESIDUALS ASCII'
	DIR[24, ] = '$SOLUTION'
	DIR[27, ] = '$DMUAI'
	DIR[28, ] = 10
	DIR[29, ] = '1e-07'
	DIR[30, ] = '1.0d-6'
	DIR[31, ] = 1
	DIR[32, ] = 0
  } else {        
    DIR[16, ] = '2 2 0 0 0'
    DIR[c(17, 18, 23, 24, 25), ] = 0
	DIR[19, ] = paste0((which(names(Grow) == strsplit(phe, "_")[[1]][1])-2), ' 0 2 2 1')
	DIR[20, ] = paste0((which(names(Grow) == strsplit(phe, "_")[[1]][2])-2), ' 0 2 2 1')
    DIR[c(21, 22), ] = '1 1'
	DIR[c(26, 29, 30), ] = ' '
	DIR[27, ] = '$RESIDUALS ASCII'
	DIR[28, ] = '$SOLUTION'
	DIR[31, ] = '$DMUAI'
	DIR[32, ] = 10
	DIR[33, ] = '1e-07'
	DIR[34, ] = '1.0d-6'
	DIR[35, ] = 1
	DIR[36, ] = 0
  }

  if(Method == 'ssgblup'){
    DIR[13, ] = '$VAR_STR 1 PGMIX 1 ASCII  ../ped.txt    ../IND_geno.txt   ../G.txt  0.05  G-ADJUST'
  }else if(Method == 'gblup'){
	DIR[13, ] = '$VAR_STR 1 COR ASCII  ../G.txt'
  }
  
  write.table(DIR, paste0(path, '/dmu.DIR'), col.names = F, row.names = F, quote = F)
}
