#' run rand mate
#'
#' @param NULL
#'
#' @return a list of crossplan and inbreeding values
#' @export
#'
#' @examples
Mate_rand <- function(maxF = maxF){
	allId <- Phen[Phen$id %in% parent@id, c('id', 'sex')]
	names(allId) <- c('Id', 'sex')
	allId = unique(allId)
	##待配母猪
	sow <- allId[allId$sex == 'F', ]
	##待配公猪
	boar <- allId[allId$sex == 'M', ]
	sor <- floor(length(sow$Id)/length(boar$Id))
	##生成所有配种可能
	plans <- data.frame(id_s = rep(sow$Id, each = nrow(boar)),
						id_b = rep(boar$Id, times = nrow(sow)))

		##匹配近交系数
		count.id <- union(sow$Id, boar$Id)
		##将需要计算亲缘关系的个体而未在子代中出现的个体加到子代中，父母设为未知
		add <- !(count.id %in% Ped[,1])
		Ped1 <- Ped[, c(1,3,2)]
		if (any(add)) {
		  Ped1 <- rbind(data.frame(id = count.id[add], father = 0, mother = 0), Ped1)
		}
		## 给系谱排序，祖先在前子代在后
		Ped1 <- Ped1[order(pedigree::orderPed(Ped1)), ]
		## 设定需要计算亲缘系数的个体
		needed <- Ped1[, 1] %in% count.id

		## 计算亲缘关系矩阵
		pedigree::makeA(Ped1, which=needed)
		## 注意，A中的V1和V2代表的是需要计算近交的个体在系谱中的位置，不是个体编号
		A <- read.table("A.txt") 
		## 换上个体编号
		A$id1 <- Ped1[A$V1, 1]
		A$id2 <- Ped1[A$V2, 1]
		## 切片
		A <- A[, c('id1', 'id2', 'V3')]
		names(A)[3] <- 'A'
		## 上面的是上三角矩阵，补成全矩阵(1、2列ID交换)
		A2 <- data.frame(id1=A$id2, id2=A$id1, A=A$A, stringsAsFactors = FALSE)
		A <- rbind(A, A2)
		A <- unique(A)

		names(A)[1:2] <- c('id_s', 'id_b')

	plans <- merge(plans, A, by=c('id_s', 'id_b'), all.x = TRUE)
	rownames(plans) <- paste0(plans$id_s, '-', plans$id_b)
	##随机配种方案
	if(maxF == 999){    #不限制近交的随机交配
		plan0 <- list()
		sow_id <- sow$Id
		Rows0 <- c()
		for(i in boar$Id){
		  if(sow_id[1] != 0){
			if(i == tail(boar$Id, 1)){
			  plan0[[i]] <- sows <- sow_id
			}else{
			  if(length(sow_id) >= sor){
				plan0[[i]] <- sows <- sample(sow_id, sor)
				sow_id <- sow_id[sow_id %in% sows == FALSE]
			  }else{
				plan0[[i]] <- sows <- sow_id
				sow_id <- 0
			  }
			}
			Rows0 <- c(Rows0, paste0(sows, '-', i))
		  }else{
			break
		  }
		}
		E0 <- sum(plans[Rows0, 'A'])
		##最终选配计划表（两列——公、母）
		crossplan <- data.frame()
		for (i in boar$Id) {
		  crossplan <- rbind(data.frame(id_b = i, id_s = plan0[[i]]),
							 crossplan)
		}
		unmated_s <- sow$Id[sow$Id %in% crossplan$id_s == FALSE]
		if(length(unmated_s) != 0){         ##如果还有未配种的母猪，随机挑选公猪与它们配种
		  plan2 <- data.frame(id_b = sample(boar$Id, length(unmated_s)), id_s = unmated_s)
		  crossplan <- rbind(crossplan, plan2)
		  plan2 <- merge(plan2, A, by=c('id_s', 'id_b'), all.x = TRUE)
		  E0 <- sum(plans[Rows0, 'A'], plan2$A)
		} 
	}else{   #set inbreeding threshold
		eligible = plans[plans$A <= 2*maxF, ]               #select eligible pairs
		re_plans = plans[plans$A > 2*maxF, ]		        #not eligible pairs
		re_plans = re_plans[order(re_plans$A, decreasing = F), ]
		allBoar = table(plans$id_b)
		allBoar = allBoar[order(allBoar, decreasing = F)] 
		allSow = unique(plans$id_s)
		useS = c()
		useB = rep(0, length(unique(plans$id_b)))               #the usage of boar
		names(useB) = names(allBoar)
		crossplan = matrix(allSow, nrow = length(allSow), ncol = 1)
		crossplan = cbind(crossplan, rep(NA, length(allSow)))
		rownames(crossplan) = allSow
		
		## assign boar for sow
		for(i in names(allBoar)){
			ava_sow = eligible[eligible$id_b == i, 'id_s']
			ava_sow = ava_sow[ava_sow%in%useS == F]           #avaliable sows for boar
			if(length(ava_sow) >= sor){
				nSow = sor
			}else{
				nSow = length(ava_sow)
			}
			useB = c(useB, nSow)
			selSow = sample(ava_sow, nSow)
			crossplan$id_b[crossplan$id_s %in% selSow] = i
			useS = c(useS, selSow)
		}
		if(length(useS)<length(allSow)){                       #check wether all sow has been assigned
			for(i in allSow[allSow %in% useS == F]){
				ava_boar = names(useB[useB<sor])
				ava_boar = re_plans[re_plans$id_s == i, 'id_b']
				if(length(ava_boar) != 0){                     #check wether there are avaliable boar
					crossplan[as.character(i), 2] = ava_boar[1]
					useB[as.character(ava_boar[1])] = useB[as.character(ava_boar[1])] + 1
				}else{
					ava_boar = eligible[eligible$id_s == i, 'id_b']
					crossplan[as.character(i), 2] = sample(ava_boar, 1)
				}
			}
		}
		E0 = merge(data.frame(crossplan), plans, by.x = c('X1','X2'), by.y = c('id_s','id_b'), sort = F, all.x = T)
		E0 = sum(E0$A)
	}
	list(crossplan = as.matrix(crossplan), Inb = E0/length(sow$Id)/2)
}

