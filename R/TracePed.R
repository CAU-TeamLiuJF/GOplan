#' Tracing ped for n generation
#'
#' @param ped pedigree file
#' @param n  the number of generation want to trace
#' @param id  the vector of individuals that wanted to trace pedigree
#' 
#' @return return a pedigree file contain exact generations
#' @export
#'
#' @examples
TracePed <- function(ped = Ped, n = 3, id = famselect@id){
	Pedi = data.frame(matrix(nrow = 0, ncol = 4))
	names(Pedi) = c("id", "father", "mother", "G")
	names(ped) = c("id", "mother", "father", "G")
	for(i in 1:n){
		Pedi = rbind(Pedi, ped[ped$id %in% id, c(1,3,2,4)])                   #加入这一世代的个体ped
		id = c(ped[ped$id %in% id, "father"], ped[ped$id %in% id, "mother"])  #将id更新为个体的父母亲
	}
	return(Pedi)
}


