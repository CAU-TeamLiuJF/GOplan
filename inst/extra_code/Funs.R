#***************************************************#
##计算基因流动法的函数
#***************************************************#
#function of multiply matrix with vecotr
ma_vec = function(matrix, vector){
	a = c(rep(0, length(vector)))
	for (i in 1:nrow(matrix)) {
		a[i] = matrix[i, ] %*% vector
	}
	return(a)
}
#***************************************************#
##funciton of generate DMU dir files
#***************************************************#

#***************************************************#
### function of calculate x^n
#***************************************************#
pow <- function(x,n){
  return(x^n)
}

#***************************************************#
### function to run_dmu
#***************************************************#
run_dmu = function(dmu_module = 'dmuai',
				   DIR_path = DIR_path){
  # 复制dmu程序和.bat或.sh脚本到工作目录
  if (Sys.info()["sysname"]=="Linux") {
	file.copy(from = system.file("extra_code","run_dmuai.sh", package = "GOplan"), to = DIR_path)
    file.copy(from = system.file("dmu","dmu1", package = "GOplan"), to = DIR_path)
    file.copy(from = system.file("dmu",dmu_module, package = "GOplan"), to = DIR_path)
	system(paste0("chmod 777 -R ", DIR_path))
  } else {
    file.copy(from = system.file("extra_code","run_dmuai.bat", package = "GOplan"), to = DIR_path)
	file.copy(from = system.file("dmu","dmu1.exe", package = "GOplan"), to = DIR_path)
    file.copy(from = system.file("dmu",paste0(dmu_module,".exe"), package = "GOplan"), to = DIR_path)
  }
  # run DMU
  if (Sys.info()["sysname"]=="Linux") {
    system(paste0(DIR_path, "/run_",dmu_module,".sh ", DIR_path, '/dmu'))
    file.remove(c("dmu1",dmu_module))
    file.remove(list.files(pattern = "*.sh"))

  } else {
    system(paste0(DIR_path, "/run_",dmu_module,".bat ", DIR_path, '/dmu'))
    file.remove(list.files(pattern = "*.exe"))
    file.remove(list.files(pattern = "*.bat"))
  }
}
#***************************************************#
### function to run_dmu in different traits
#***************************************************#
custom.function=function(x){
    if(file.exists(paste0(x[3], '/', x[2], "/", x[1]))==FALSE){
		dir.create(paste0(x[3], '/', x[2], "/", x[1]))
	}
    setwd(paste0(x[3], '/', x[2], "/", x[1]))
    run_dmu(DIR_path = paste0(x[3], '/', x[2], '/', x[1]))
}
for(i in c('ma_vec', 'custom.function', 'pow', 'run_dmu')){
	assign(i , get(i), envir = .GlobalEnv)
}

