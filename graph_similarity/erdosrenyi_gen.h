#ifndef ER_GRAPH_
#define ER_GRAPH_

#include "graph_gen_alg.h"

class ErdosRenyi_Gen : public Graph_Gen_Alg {
public:
ErdosRenyi_Gen() {}
~ErdosRenyi_Gen() {}
 std::vector< std::vector<int> > gen_unweighted_graph(const int n, const std::vector<double>& params);
 std::vector< std::vector<int> > gen_unweighted_graph(const int n, const double p);
 std::vector< std::vector<double> > gen_graph(const int n, const std::vector<double>& params);
 std::vector< std::vector<double> > gen_graph(const int n, const double p);
};

#endif
