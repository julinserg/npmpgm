#ifndef GRAPH_GEN_ALG_
#define GRAPH_GEN_ALG_

#include <vector>
#include <random>

class Graph_Gen_Alg {
 public:
  Graph_Gen_Alg();
  virtual ~Graph_Gen_Alg();
  virtual std::vector< std::vector<int> > gen_unweighted_graph(const int n, const std::vector<double>& params) = 0;
  virtual std::vector< std::vector<double> > gen_graph(const int n, const std::vector<double>& params) = 0;
 protected:
  double gen_uniform_rn() const;
 private:
  std::mt19937 *mt;
  double normalization;
};

#endif
