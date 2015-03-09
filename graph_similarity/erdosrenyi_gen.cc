#include <cmath>
#include "erdosrenyi_gen.h"

std::vector< std::vector<int> > ErdosRenyi_Gen::gen_unweighted_graph(const int n, const std::vector<double>& params) {
  return gen_unweighted_graph(n, params[0]);
}

std::vector< std::vector<int> > ErdosRenyi_Gen::gen_unweighted_graph(const int n, const double p) {
  const double tol = 1e-6;
  if(p < tol) {
    std::vector< std::vector<int> > A(n, std::vector<int>(n, 0));
    return A;
  }
  else if(p < 1 + tol && p > 1 - tol) {
    std::vector< std::vector<int> > A(n, std::vector<int>(n, 1));
    return A;
  }
  else {
    std::vector< std::vector<int> > A(n, std::vector<int>(n, 0));
    //allocate and init arrays to zero
    int v = 1;
    int w = -1;
    while(v <= n) {
      double rv1 = gen_uniform_rn();
      // previously added some ROUND_CONST
      // but I believe it was unecessary
      w = (int) (w + 1 + floor(log(1-rv1)/log(1-p)));
      while((w >= v) && (v <= n)) {
	w -= v;
	v++;
      }
      if(v < n) {
	A[v][w] = 1;
	A[w][v] = 1;
      }
    }
    return A;
  }
}

std::vector< std::vector<double> > ErdosRenyi_Gen::gen_graph(const int n, const std::vector<double>& params) {
  return gen_graph(n, params[0]);
}

std::vector< std::vector<double> > ErdosRenyi_Gen::gen_graph(const int n, const double p) {
  std::vector< std::vector<int> > A_temp = gen_unweighted_graph(n, p);
  std::vector< std::vector<double> > A(n);
  for(int i = 0; i < n; i++) {
    A[i] = std::vector<double>(A_temp[i].begin(), A_temp[i].end());
  }
  return A;
}
