#include <cmath>
#include "chunglu_gen.h"

std::vector< std::vector<int> > ChungLu_Gen::gen_unweighted_graph(const int n, const std::vector<double>& params) {
  return gen_unweighted_graph(n, params[0], params[1]);
}

std::vector< std::vector<int> > ChungLu_Gen::gen_unweighted_graph(const int n, const double p, const double r) {
  std::vector<double> weights(n);
  double normalization = 0;
  for(int i = 1; i <= n; i++) {
    normalization += pow(i, r);
  }
  normalization *= n*p/pow(n, r);
  normalization = sqrt(normalization);

  // normalize weights so that the division by the sum
  // is unecessary when assigning A_{ij} probabilities
  for(int i = 0; i < n; i++) {
    weights[i] = n*p*pow((i+1.0)/n, r)/normalization;
  }

  std::vector< std::vector<double> > P(n, std::vector<double>(n));
  for(int i = 0; i < n; i++) {
    for(int j = i; j < n; j++) {
      if(weights[i]*weights[j] < 1) {
	P[i][j] = weights[i]*weights[j];
      }
      else {
	P[i][j] = 1;
      }
    }
  }

  std::vector< std::vector<int> > A(n, std::vector<int>(n));
  for(int i = 0; i < n; i++) {
    for(int j = i; j < n; j++) {
      if(gen_uniform_rn() < P[i][j]) {
	A[i][j] = 1;
	A[j][i] = 1;
      }
      else {
	A[i][j] = 0;
	A[j][i] = 0;
      }
    }
  }
  
  return A;
}


std::vector< std::vector<double> > ChungLu_Gen::gen_graph(const int n, const std::vector<double>& params) {
  return gen_graph(n, params[0], params[1]);
}

std::vector< std::vector<double> > ChungLu_Gen::gen_graph(const int n, const double p, const double r) {
  std::vector< std::vector<int> > A_temp = gen_unweighted_graph(n, p, r);
  std::vector< std::vector<double> > A(n);
  for(int i = 0; i < n; i++) {
    A[i] = std::vector<double>(A_temp[i].begin(), A_temp[i].end());
  }
  return A;
}
