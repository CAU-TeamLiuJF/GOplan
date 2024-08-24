#define ARMA_USE_BLAS
#define ARMA_USE_LAPACK
#define ARMA_USE_OPENMP
#define ARMA_DONT_USE_WRAPPER

//[[Rcpp::depends(RcppArmadillo)]]


#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' @export
// [[Rcpp::export]]
List G_mat_cpp(arma::Mat<int> & M, bool inverse=false) {
  //M: ind_num * snp_num, row_name = ind_id, col_name = snp_id
  // VanRaden (2008)
  // G = ZZ' / (2 \sum p_i (1-p_i))

  arma::Mat<double> M1=arma::conv_to<arma::Mat<double>>::from(M);
  int ind_num=M1.n_rows;
  int snp_num=M1.n_cols;

  Rcpp::Rcout<<"ind_num: "<<ind_num<<"  snp_num: "<<snp_num<<endl;

  // Initialisation
  arma::Row<double> p(snp_num),q(snp_num);
  arma::Mat<double> G,G_inv;
  Rcpp::List result;

  for(int i=0; i < snp_num; i++) {
    arma::vec col=M1.col(i);
    arma::uvec non_missing=arma::find(col != -9);
    // indices of non-missing elements
    double col_mean;
    if(non_missing.n_elem > 0) {
      col_mean=arma::mean(col(non_missing));
      // column mean excluding -9 elements
    } else {
      col_mean=0;
      // if all elements are -9, set column mean to 0
    } col.replace(-9,col_mean);
    // replace -9 with column mean
    p[i]=col_mean;
    M1.col(i)=col-col_mean;
  }

  p=0.5*p;
  q=1-p;

  // Formulate
  G = M1 * M1.t() / sum(2*p%q);

  // Inverse   positive-definite matrix RcppEigen
  if(inverse==true){
    arma::vec eigval = arma::eig_sym(G);
    if(all(eigval > 0)){
      Rcpp::Rcout<<"[G]: positive-definite matrix"<<endl;
      G_inv=inv(G);
    }else{
      Rcpp::Rcout<<"[G]: not positive-definite matrix, add 0.0001 to diag"<<endl;
      G.diag()=G.diag()+0.0001;
      G_inv = inv(G);
      G.diag()=G.diag()-0.0001;
    }
  }else{
    Rcpp::Rcout<<"[G]: do not inverse matrix"<<endl;
  }

  result = Rcpp::List::create(Rcpp::Named("G") = G, Rcpp::Named("G_inv") = G_inv);

  return result;
}

// 包装函数的声明
extern "C" SEXP G_mat_cpp(SEXP MSEXP, SEXP inverseSEXP);

// 包装函数的实现
extern "C" SEXP G_mat_cpp(SEXP MSEXP, SEXP inverseSEXP) {
    BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    
    // 将 SEXP 转换为 arma::Mat<int>
    arma::Mat<int> M = Rcpp::as<arma::Mat<int>>(MSEXP);
    bool inverse = Rcpp::as<bool>(inverseSEXP);
    
    // 调用实际的 C++ 函数
    rcpp_result_gen = Rcpp::wrap(G_mat_cpp(M, inverse));
    return rcpp_result_gen;
    END_RCPP
}
