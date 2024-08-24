#include <Rcpp.h>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include <ctime>
#include "fun.h"  // 包含fun.h

const double DOUBLE_EPS = 1e-15;

// 使用Rcpp进行自动注册
using namespace Rcpp;

// [[Rcpp::export]]
void mc_cpp(int nplan, std::string fPath) {
  int nrows = nplan;
  time_t cur_time;

  std::ifstream plans(fPath); // 读取文件
  if (!plans.is_open()) {
    Rcpp::Rcerr << "Error opening file: " << strerror(errno) << std::endl;
    return;
  }

  std::vector<std::vector<int>> id_bs(nrows, std::vector<int>(2));
  std::vector<double> coA(nrows);
  for (int i = 0; i < nrows; i++) {
    plans >> id_bs[i][0] >> id_bs[i][1] >> coA[i];
  }
  plans.close();

  std::vector<int> id_s(nrows), id_b(nrows);
  for (int i = 0; i < nrows; i++) {
    id_b[i] = id_bs[i][0];
    id_s[i] = id_bs[i][1];
  }

  auto remove_duplicates = [](std::vector<int>& vec) {
    std::sort(vec.begin(), vec.end());
    vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
  };

  remove_duplicates(id_s);
  remove_duplicates(id_b);

  int ns = id_s.size(), nb = id_b.size();
  std::vector<std::vector<double>> reA(ns, std::vector<double>(nb));
  for (int i = 0; i < nrows; ++i) {
    reA[search(id_bs[i][1], id_s)][search(id_bs[i][0], id_b)] = coA[i];
  }

  int sor = ns / nb;
  std::vector<std::vector<int>> bestMates(nb, std::vector<int>(2 * sor));
  std::vector<std::vector<int>> mates(nb, std::vector<int>(2 * sor));
  std::vector<int> mate_rows(nb);
  if ((ns % nb) * 1.0 / nb > 0.5) ++sor;

  for (int j = 0; j < nb - 1; ++j) {
    auto tp2 = N_array(j * sor, (j + 1) * sor);
    for (int i = 0; i < sor; i++) {
      mates[j][i] = tp2[i];
    }
    mate_rows[j] = sor;
  }
  auto tp2 = N_array((nb - 1) * sor, ns);
  for (int i = 0; i < ns - (nb - 1) * sor; i++) {
    mates[nb - 1][i] = tp2[i];
  }
  mate_rows[nb - 1] = ns - (nb - 1) * sor;

  double E0 = 0;
  for (int i = 0; i < nb; i++) {
    for (int j = 0; j < mate_rows[i]; j++) {
      E0 += reA[mates[i][j]][i];
    }
  }

  Rcpp::Rcout << E0 << std::endl;
  double bestE = E0;

  int maxEva = 100 * ns;
  int maxRe = 10 * ns;
  int Neva = 0, round = 0;
  int Nre = 1;
  int rd_s1, rd_s2, rand_s1, rand_s2;
  double delta, p, rand_p, Tem = 1.0;

  time(&cur_time);
  Rcpp::Rcout << "start annealing: " << ctime(&cur_time) << std::endl;

  while (Nre > 0) {
    Neva = 0;
    Nre = 0;
    while (Neva < maxEva && Nre < maxRe) {
      auto rand_b = rand_array(nb);
      rd_s1 = rand() % mate_rows[rand_b[0]];
      rd_s2 = rand() % mate_rows[rand_b[1]];
      rand_s1 = mates[rand_b[0]][rd_s1];
      rand_s2 = mates[rand_b[1]][rd_s2];
      delta = reA[rand_s1][rand_b[1]] + reA[rand_s2][rand_b[0]] - reA[rand_s1][rand_b[0]] - reA[rand_s2][rand_b[1]];
      p = exp(-delta / Tem);
      rand_p = static_cast<double>(rand()) / (RAND_MAX + 1.0);
      if (delta < 0.0 || (delta > 0.0 && rand_p < p)) {
        mates[rand_b[0]][rd_s1] = rand_s2;
        mates[rand_b[1]][rd_s2] = rand_s1;

        E0 += delta;
        ++Nre;
        ++Neva;
      } else {
        ++Neva;
      }
      if (E0 < bestE) {
        bestE = E0;
        bestMates = mates;
      }
      ++round;
    }
    Tem *= 0.9;
    if (fabs(Tem) <= DOUBLE_EPS) break;
  }

  time(&cur_time);
  Rcpp::Rcout << "round: " << round << " finish time: " << ctime(&cur_time) << std::endl;

  std::ofstream outf("reMate.txt");
  outf << "id_b id_s\n";
  for (int i = 0; i < nb; i++) {
    for (int j = 0; j < mate_rows[i]; j++) {
      outf << id_b[i] << " " << id_s[bestMates[i][j]] << "\n";
    }
  }
  outf << "averageInb " << bestE / (ns * 2.0) << "\n";
  outf.close();
}
