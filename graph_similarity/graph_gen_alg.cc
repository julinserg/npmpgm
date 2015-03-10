#include <chrono>
#include "graph_gen_alg.h"

Graph_Gen_Alg::Graph_Gen_Alg() {
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  mt = new std::mt19937(seed);
  normalization = (double) mt->max()+1;
}

Graph_Gen_Alg::~Graph_Gen_Alg() {
  delete mt;
}

double Graph_Gen_Alg::gen_uniform_rn() const {
  return (*mt)()/(normalization);
}



