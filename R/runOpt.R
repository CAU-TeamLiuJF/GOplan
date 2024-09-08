#' This function runs breeding programs' optimization task
#' 
#' @param prm_path the full path to the prm file
#'
#' @return
#' @export
runOpt <- function(prm_path = prm_path){
	library(AlphaSimR)
  library(parallel)
  library(writexl)
  library(readxl)
  library(tidyverse)
  library(foreach)
  library(stringr)
  library(mlrMBO)
  library(pedigree)
  library(progress)
  
	knitStartTime <- Sys.time()
	options (warn = -1)
	task <- 3
	assign('task', task, envir = .GlobalEnv)
	### load prm files(.GlobalEnv)
	readPrm(prm_path, task)
	### load functions(.GlobalEnv)
	source(system.file("extra_code","Funs.R", package = "GOplan"),encoding = 'UTF-8')
	### construct founder pop(save as RData)
	if(createFounder == TRUE){
		Founders(Breed)
	}

	optP <- list(
	  kernel = 'gauss', # kernel of the gaussian process
	  final.method = 'best.predicted', # see mlrmbo documentation (`??makeMBOControl`) and below
	  acquFunction = makeMBOInfillCritEI(),
	  filter.proposed.points = TRUE, # should bayesian opt avoid sampling parameters too close to each-other
	  filter.proposed.points.tol = 0.001, # threshold defining if 2 points are too close to each-other
	  propose.points = nPoints, # number of points evaluated at each iteration
	  nCpusProposePoints = nPoints, # number of core to use for evaluating the points at each iteration
	  totalIter = nIter, # total number of iteration
	  # time.budget = 60*30, # total time available for the optimization
	  initTrainSize = nrow(Parset)*2 # initial training data size
	  #funEval = 1, # see below
	  #nCpusFunEval = 1 # see below
	)

# mlrMBO setup ----
  ## Learner ----
  learner = makeLearner(
    "regr.km",
    predict.type = "se",
    covtype = optP$kernel,
    nugget.estim = TRUE,
    control = list(trace = FALSE),
    scaling = TRUE
  )

  ## General Control ----
  if (is.null(optP$save.on.disk.at.time)) {
    optP$save.on.disk.at.time <- Inf
  }
  control = makeMBOControl(
    n.objectives = 1,
    propose.points = optP$propose.points,  
    final.method = optP$final.method,      #best.predicted
    final.evals = 0,
    y.name = "profit",
    on.surrogate.error = "warn",
    save.on.disk.at.time = optP$save.on.disk.at.time,
    save.file.path = fp$out_path
  )

  ## Aquisition function ----
  control = setMBOControlInfill(
    control,
    crit = optP$acquFunction,
    filter.proposed.points = optP$filter.proposed.points,
    filter.proposed.points.tol = optP$filter.proposed.points.tol
  )

  ## Multi-points proposition ----
  control = setMBOControlMultiPoint(
    control,
    method = "cl",
    cl.lie = min
  )

  ## Termination condition ----
  control = setMBOControlTermination(
    control,
    iters = optP$totalIter
  )
  # Create optimized parameters: ----
  if(Ctype == 33){
  optParams <- makeParamSet(
    makeIntegerParam(id = "Asec",
                     lower = Parset['Asec', 'lower'],
                     upper = Parset['Asec', 'upper'],
                     default = round(0.5*(Parset['Asec', 'lower']+Parset['Asec', 'upper']))),
    makeIntegerParam(id = "Bsec",
                     lower = Parset['Bsec', 'lower'],
                     upper = Parset['Bsec', 'upper'],
                     default = round(0.5*(Parset['Bsec', 'lower']+Parset['Bsec', 'upper']))),
    makeIntegerParam(id = "Csec",
                     lower = Parset['Csec', 'lower'],
                     upper = Parset['Csec', 'upper'],
                     default = round(0.5*(Parset['Csec', 'lower']+Parset['Csec', 'upper']))),
    makeIntegerParam(id = "YABd",
                     lower = Parset['YABd', 'lower'],
                     upper = Parset['YABd', 'upper'],
                     default = round(0.5*(Parset['YABd', 'lower']+Parset['YABd', 'upper']))),
    makeIntegerParam(id = "NAd",
                     lower = NA_min,
                     upper = Parset['NAd', 'upper'],
                     default = round(0.5*(NA_min+Parset['NAd', 'upper']))),
	makeIntegerParam(id = "NBd",
                     lower = Parset['NBd', 'lower'],
                     upper = Parset['NBd', 'upper'],
                     default = round(0.5*(Parset['NBd', 'lower']+Parset['NBd', 'upper']))),
    makeIntegerParam(id = "NCd",
                     lower = Parset['NCd', 'lower'],
                     upper = Parset['NCd', 'upper'],
                     default = round(0.5*(Parset['NCd', 'lower']+Parset['NCd', 'upper'])))
  )}else if(Ctype == 32){
  optParams <- makeParamSet(
    makeIntegerParam(id = "Asec",
                     lower = Parset['Asec', 'lower'],
                     upper = Parset['Asec', 'upper'],
                     default = round(0.5*(Parset['Asec', 'lower']+Parset['Asec', 'upper']))),
    makeIntegerParam(id = "Csec",
                     lower = Parset['Csec', 'lower'],
                     upper = Parset['Csec', 'upper'],
                     default = round(0.5*(Parset['Csec', 'lower']+Parset['Csec', 'upper']))),
    makeIntegerParam(id = "YABd",
                     lower = Parset['YABd', 'lower'],
                     upper = Parset['YABd', 'upper'],
                     default = round(0.5*(Parset['YABd', 'lower']+Parset['YABd', 'upper']))),
    makeIntegerParam(id = "NAd",
                     lower = NA_min,
                     upper = Parset['NAd', 'upper'],
                     default = round(0.5*(NA_min+Parset['NAd', 'upper']))),
	makeIntegerParam(id = "NBd",
                     lower = Parset['NBd', 'lower'],
                     upper = Parset['NBd', 'upper'],
                     default = round(0.5*(Parset['NBd', 'lower']+Parset['NBd', 'upper']))),
    makeIntegerParam(id = "NCd",
                     lower = Parset['NCd', 'lower'],
                     upper = Parset['NCd', 'upper'],
                     default = round(0.5*(Parset['NCd', 'lower']+Parset['NCd', 'upper'])))
  )}else if(Ctype == 31){
  optParams <- makeParamSet(
    makeIntegerParam(id = "Asec",
                     lower = Parset['Asec', 'lower'],
                     upper = Parset['Asec', 'upper'],
                     default = round(0.5*(Parset['Asec', 'lower']+Parset['Asec', 'upper']))),
    makeIntegerParam(id = "YABd",
                     lower = Parset['YABd', 'lower'],
                     upper = Parset['YABd', 'upper'],
                     default = round(0.5*(Parset['YABd', 'lower']+Parset['YABd', 'upper']))),
    makeIntegerParam(id = "NAd",
                     lower = NA_min,
                     upper = Parset['NAd', 'upper'],
                     default = round(0.5*(NA_min+Parset['NAd', 'upper']))),
	makeIntegerParam(id = "NBd",
                     lower = Parset['NBd', 'lower'],
                     upper = Parset['NBd', 'upper'],
                     default = round(0.5*(Parset['NBd', 'lower']+Parset['NBd', 'upper']))),
    makeIntegerParam(id = "NCd",
                     lower = Parset['NCd', 'lower'],
                     upper = Parset['NCd', 'upper'],
                     default = round(0.5*(Parset['NCd', 'lower']+Parset['NCd', 'upper'])))
  )}else if(Ctype == 42){
  optParams <- makeParamSet(
    makeIntegerParam(id = "Asec",
                     lower = Parset['Asec', 'lower'],
                     upper = Parset['Asec', 'upper'],
                     default = round(0.5*(Parset['Asec', 'lower']+Parset['Asec', 'upper']))),
    makeIntegerParam(id = "Bsec",
                     lower = Parset['Bsec', 'lower'],
                     upper = Parset['Bsec', 'upper'],
                     default = round(0.5*(Parset['Bsec', 'lower']+Parset['Bsec', 'upper']))),
    makeIntegerParam(id = "Csec",
                     lower = Parset['Csec', 'lower'],
                     upper = Parset['Csec', 'upper'],
                     default = round(0.5*(Parset['Csec', 'lower']+Parset['Csec', 'upper']))),
	makeIntegerParam(id = "Dsec",
                     lower = Parset['Dsec', 'lower'],
                     upper = Parset['Dsec', 'upper'],
                     default = round(0.5*(Parset['Dsec', 'lower']+Parset['Dsec', 'upper']))),
    makeIntegerParam(id = "YABd",
                     lower = Parset['YABd', 'lower'],
                     upper = Parset['YABd', 'upper'],
                     default = round(0.5*(Parset['YABd', 'lower']+Parset['YABd', 'upper']))),
    makeIntegerParam(id = "YCDs",
                     lower = Parset['YCDs', 'lower'],
                     upper = Parset['YCDs', 'upper'],
                     default = round(0.5*(Parset['YCDs', 'lower']+Parset['YCDs', 'upper']))),
	makeIntegerParam(id = "NAd",
                     lower = NA_min,
                     upper = Parset['NAd', 'upper'],
                     default = round(0.5*(NA_min+Parset['NAd', 'upper']))),
	makeIntegerParam(id = "NBd",
                     lower = Parset['NBd', 'lower'],
                     upper = Parset['NBd', 'upper'],
                     default = round(0.5*(Parset['NBd', 'lower']+Parset['NBd', 'upper']))),
    makeIntegerParam(id = "NCd",
                     lower = Parset['NCd', 'lower'],
                     upper = Parset['NCd', 'upper'],
                     default = round(0.5*(Parset['NCd', 'lower']+Parset['NCd', 'upper']))),
	makeIntegerParam(id = "NDd",
                     lower = Parset['NDd', 'lower'],
                     upper = Parset['NDd', 'upper'],
                     default = round(0.5*(Parset['NDd', 'lower']+Parset['NDd', 'upper'])))
  )}else if(Ctype == 41){
  optParams <- makeParamSet(
    makeIntegerParam(id = "Asec",
                     lower = Parset['Asec', 'lower'],
                     upper = Parset['Asec', 'upper'],
                     default = round(0.5*(Parset['Asec', 'lower']+Parset['Asec', 'upper']))),
    makeIntegerParam(id = "Csec",
                     lower = Parset['Csec', 'lower'],
                     upper = Parset['Csec', 'upper'],
                     default = round(0.5*(Parset['Csec', 'lower']+Parset['Csec', 'upper']))),
    makeIntegerParam(id = "YABd",
                     lower = Parset['YABd', 'lower'],
                     upper = Parset['YABd', 'upper'],
                     default = round(0.5*(Parset['YABd', 'lower']+Parset['YABd', 'upper']))),
    makeIntegerParam(id = "YCDs",
                     lower = Parset['YCDs', 'lower'],
                     upper = Parset['YCDs', 'upper'],
                     default = round(0.5*(Parset['YCDs', 'lower']+Parset['YCDs', 'upper']))),
	makeIntegerParam(id = "NAd",
                     lower = NA_min,
                     upper = Parset['NAd', 'upper'],
                     default = round(0.5*(NA_min+Parset['NAd', 'upper']))),
	makeIntegerParam(id = "NBd",
                     lower = Parset['NBd', 'lower'],
                     upper = Parset['NBd', 'upper'],
                     default = round(0.5*(Parset['NBd', 'lower']+Parset['NBd', 'upper']))),
    makeIntegerParam(id = "NCd",
                     lower = Parset['NCd', 'lower'],
                     upper = Parset['NCd', 'upper'],
                     default = round(0.5*(Parset['NCd', 'lower']+Parset['NCd', 'upper']))),
	makeIntegerParam(id = "NDd",
                     lower = Parset['NDd', 'lower'],
                     upper = Parset['NDd', 'upper'],
                     default = round(0.5*(Parset['NDd', 'lower']+Parset['NDd', 'upper'])))
  )}

  # 初始化进度条
  pb <- progress_bar$new(
    format = "  Running [:bar] :percent eta: :eta",
    total = nIter*nPoints+optP$initTrainSize,  # 假设总共有100次调用，实际值根据需要调整
    width = 60
    )

# 包装原始的 OptSim 函数以增加进度条更新
  OptSimWithProgress <- function(x) {
    # 更新进度条
    pb$tick()
  
    # 调用原始的 OptSim 函数
    return(OptSim(x))

    if (pb$is_finished) {
      cat("\n100% Complete\n")
    }
  }

  fun <- makeSingleObjectiveFunction(
    name = "breedingSimulation",
    id = "breedingSimulation",
    description = "breedingSimulation",
    fn = OptSimWithProgress,
    vectorized = FALSE,
    par.set = optParams,
    noisy = TRUE,
    minimize = FALSE)

  # Initial design ----(Latin square)
  initDes = generateDesign(n = optP$initTrainSize,
                           par.set = optParams,
                           fun = lhs::improvedLHS)

  # Bayesian optimization  ----
  sink(paste0(out_path, "/OptLog.txt"))
  run <- mbo(
    fun = fun,
    design = initDes,
    learner = learner,
    control = control,
    show.info = TRUE
  )

# Generated outputs ----
  results <- as.data.frame(run$opt.path)
  write.csv(results[, 1:(ncol(results)-10)], paste0(out_path, 'Opt_results.csv'), row.names = F, quote = F)
  #save(results, file = paste0(out_path, 'results.RData'))
	sink()
	cat("Analyse finish!\n")
	cat("Results generated in:", out_path, "\n")
	print(Sys.time() - knitStartTime)
	log_name = paste0(out_path, "/OptLog.txt")
	#getOptRes(log_name, out_path)
}




