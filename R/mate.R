#' prepare files for MC mate or others
#'
#' @param NULL
#'
#' @return
#' @export
#'
#' @examples
Mate_file <- function(){
	ebvAll <- ebv_all
	allId <- Phen[Phen$id %in% parent@id, c('id', 'sex')]
	names(allId) <- c('Id', 'sex')
	##待配母猪
	sow <- allId[allId$sex == 'F', ]
	##待配公猪
	boar <- allId[allId$sex == 'M', ]
	Sor <- floor(length(sow$Id)/length(boar$Id))
	##生成所有配种可能
	plans <- data.frame(id_s = rep(sow$Id, each = nrow(boar)),
						id_b = rep(boar$Id, times = nrow(sow)))

	if(Gmate == FALSE){
		##匹配近交系数
		count.id <- union(sow$Id, boar$Id)
		##将需要计算亲缘关系的个体而未在子代中出现的个体加到子代中，父母设为未知
		add <- !(count.id %in% Ped[,1])
		Ped1 <- Ped[, c(1,3,2)]
		if (any(add)) {
		  Ped1 <- rbind(data.frame(id = count.id[add], sire = 0, dam = 0), Ped1)
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

		names(A)[1:2] <- c('id_b', 'id_s')
	}else{
		g_matrix_cal(snp_file = 'Geno',
								inverse = FALSE,
								prefix_name = paste0(rdmu_path, "/", g))
		A <- read.table("../G.txt")
		names(A) <- c('id1', 'id2', 'A')
		## 上面的是上三角矩阵，补成全矩阵(1、2列ID交换)
		A2 <- data.frame(id1=A$id2, id2=A$id1, A=A$A, stringsAsFactors = FALSE)
		A <- rbind(A, A2)
		A <- unique(A)
		names(A) <- c('id_b', 'id_s', 'A')
	}
	
	plans <- merge(plans, A, by=c('id_b', 'id_s'), all.x = TRUE)
	plans <- plans[plans$id_b %in% boar$Id, ]
	plans <- plans[plans$id_s %in% sow$Id, ]
	if(Mate %in% c('rand', 'MC') == FALSE){
	##获得父母对后代的预估ebv
		rownames(ebvAll) <- as.character(ebvAll$Id)
		plans <- within(plans, {b_index = ebvAll[plans$id_b, 'Index']
								s_index = ebvAll[plans$id_s, 'Index']})
	}
	write.table(plans, 'plans.txt', col.names=F, row.names=F, quote=F)
	assign('Nrows', nrow(plans), envir = .GlobalEnv)
}



















