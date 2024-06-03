#' construct G matrix
#'
#' @param snp_file which plink file need to calculate g matrix, .raw format
#' @param inverse need to calculate inverse g matrix or not
#' @param prefix_name the prefix name of output file's path
#'
#' @return
#' @export
#'
#' @examples
g_matrix_cal <- function(snp_file = 'Geno',
                         inverse = FALSE,
                         prefix_name = NULL) {

  suppressPackageStartupMessages(require(dplyr))

  # 读取raw格式基因型文件
  geno_raw <- get(snp_file)
  ind_id <- data.frame(id = geno_raw$id, num = 1:nrow(geno_raw))
  geno_raw <- geno_raw[, 7:ncol(geno_raw)]
  geno_raw <- as.matrix(geno_raw)

  if (inverse == FALSE) {
    gmat <- G_mat_cpp(geno_raw, inverse = inverse)$G
  } else if(inverse == TRUE) {
    gmat <- G_mat_cpp(geno_raw, inverse = inverse)$G_inv
  } else { stop("  [inverse] input parameter error!!!  \n\n") }

  # lower.tri函数获取下三角矩阵的逻辑值
  ind <- lower.tri(gmat, diag = TRUE)

  gmat_col_3_form <- data.frame(row = row(gmat)[ind], col = col(gmat)[ind], val = c(gmat)[ind]) %>%
    dplyr::left_join(ind_id, by = c('row' = 'num')) %>%
    dplyr::left_join(ind_id, by = c('col' = 'num')) %>%
    dplyr::select(4,5,3) %>%
    dplyr::rename(id1 = 1, id2 = 2)

    # 输出G矩阵
    G <- dplyr::filter(gmat_col_3_form, !is.na(id1) & !is.na(id2))
    if(is.null(prefix_name) == F){
	data.table::fwrite(x = G, file = paste0(prefix_name, "/G.txt"), quote = F, sep = " ", col.names = F, row.names = F)
    write.table(unique(G$id1), file = paste0(prefix_name, "/IND_geno.txt"), quote = F, sep = " ", col.names = F, row.names = F)
	}else{
	data.table::fwrite(x = G, file = paste0("G.txt"), quote = F, sep = " ", col.names = F, row.names = F)
    write.table(unique(G$id1), file = paste0("IND_geno.txt"), quote = F, sep = " ", col.names = F, row.names = F)
	}
}
