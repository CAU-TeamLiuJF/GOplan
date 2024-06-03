#' extract opt outlog results
#'
#' @param res the output results of runOpt()
#' @param out_path path of output results
#'
#' @return
#' @export
#'
#' @examples
getOptRes <- function(res = res, 
					  out_path = out_path){
	res <- res[substr(res$V1, 1, 5) == '[mbo]', ]

	a <- list()
	for(i in 1:length(res)){
		b <- unlist(strsplit(res[i], ' '))
		b <- gsub(":", '', b)
		b <- gsub(";", '', b)
		b <- gsub("Asec=", '', b)
		b <- gsub("Bsec=", '', b)
		b <- gsub("Csec=", '', b)
		b <- gsub("YABd=", '', b)
		b <- gsub("NAd=", '', b)
		b <- gsub("NBd=", '', b)
		b <- gsub("NCd=", '', b)
		b <- b[is.null(b) == FALSE]
		a[[i]] <- as.data.frame(t(data.frame(b[c(2:9, 13)])))
	}
	b <- data.table::rbindlist(a)
	names(b) <- c('Nmbo', 'Asec', 'Bsec', 'Csec', 'YABd', 'NAd', 'NBd', 'NCd', 'profit')
	write.csv(b, paste0(out_path, 'mboRes.csv'), row.names = F, quote = F)
}
