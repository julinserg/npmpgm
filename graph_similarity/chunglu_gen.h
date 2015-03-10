#ifndef CHUNGLU_GRAPH_
#define CHUNGLU_GRAPH_

#include "graph_gen_alg.h"

class ChungLu_Gen : public Graph_Gen_Alg {
public:
ChungLu_Gen() {}
~ChungLu_Gen() {}
 std::vector< std::vector<int> > gen_unweighted_graph(const int n, const std::vector<double>& params);
 std::vector< std::vector<int> > gen_unweighted_graph(const int n, const double p, const double r);
 std::vector< std::vector<double> > gen_graph(const int n, const std::vector<double>& params);
 std::vector< std::vector<double> > gen_graph(const int n, const double p, const double r);
};

#endif
