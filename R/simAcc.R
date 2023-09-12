#' ebv simulation according accuracy
#' 
#' @param NULL
#'
#' @return a dataframe of simulated ebv
#' @export
#'
#' @examples
SimAcc <- function(){
	ebv <- TBV <- as.data.frame(gv(famselect))
	ebv$Index <- 0
	ebv$Id <- famselect@id

	a <- rnorm(nrow(TBV), mean = 0, sd = 1)
	for(i in TRAITS){
		Racc <- acc[1, i]
		u1 <- mean(TBV[, i])
		s1 <- sd(TBV[, i])

		ebv[, i] <- u1 + Racc*(TBV[, i] - u1) + sqrt(1-Racc*Racc)*a*s1
		ebv$Index <- as.vector(scale(ebv[, i]) * w[which(TRAITS == i)] + ebv$Index)
	}
	return(ebv)
}
