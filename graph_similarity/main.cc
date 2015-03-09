#include <sstream>
#include <iostream>
#include "graph_gen_alg.h"
#include "chunglu_gen.h"
#include "erdosrenyi_gen.h"
#include "embeddings.h"
// #include "kernels.h"
#include "gaussian_kernel.h"
#include "dmaps.h"
#include "util_fns.h"

// // for testing arnoldi
// #include "dmaps_EIG.h"

int main(int argc, char** argv) {

  std::cout << "<---------------------------------------->" << std::endl;

  const int graph_size = atoi(argv[2]);
  std::string graph_type;
  Graph_Gen_Alg* alg;
  std::vector< std::vector<double> > graph_params;
  int n_pts;

  double epsilon = 1e-3;

  if(std::string(argv[1]) == "-er") {
    graph_type = "ER";
    alg = new ErdosRenyi_Gen;
    // ErdosRenyi_Gen alg;
    const int n_pvals = atoi(argv[3]);
    n_pts = n_pvals;
    const double pmin = 0;
    const double pmax = 1;
    const double dp = (pmax - pmin)/(n_pvals - 1);
    graph_params = std::vector< std::vector<double> >(n_pts, std::vector<double>(1));
    for(int i = 0; i < n_pvals; i++) {
      graph_params[i][0] = pmin + i*dp;
    }

    std::cout << "--> Generating " << n_pts << " Erdos-Renyi graphs" << std::endl;
    std::cout << "---> " << "graph size is: " << graph_size << std::endl;
    std::cout << "---> " << n_pvals << " p values between " << pmin << " and " << pmax << std::endl;
  }
  else if(std::string(argv[1]) == "-cl") {
    // ChungLu_Gen alg;
    graph_type = "CL";
    alg = new ChungLu_Gen;

    const int n_pvals = atoi(argv[3]);
    double pmin = 0.5;
    double pmax = 1.0;
    double dp = (pmax - pmin)/(n_pvals - 1);

    const int n_rvals = atoi(argv[4]);
    double rmin = 0.0;
    double rmax = 0.5;
    double dr = (rmax - rmin)/(n_rvals - 1);
    
    epsilon = atof(argv[5]);

    n_pts = n_pvals*n_rvals;
    graph_params = std::vector< std::vector<double> >(n_pts, std::vector<double>(2));
    int counter = 0;
    for(int i = 0; i < n_pvals; i++) {
      for(int j = 0; j < n_rvals; j++) {
	graph_params[counter][0] = pmin + i*dp;
	graph_params[counter++][1] = rmin + j*dr;
      }
    }

    std::cout << "--> Generating " << n_pts << " Chung-Lu graphs" << std::endl;
    std::cout << "----> " << "graph size is: " << graph_size << std::endl;
    std::cout << "----> " << n_pvals << " p values between " << pmin << " and " << pmax << std::endl;
    std::cout << "----> " << n_rvals << " r values between " << rmin << " and " << rmax << std::endl;
  }

  const int n_spec_params = 100;
  std::vector<double> spectral_params(n_spec_params);
  const double spec_param_max = 0.01;
  const double spec_param_min = 0.0001;
  const double dspec_param = (spec_param_max - spec_param_min)/(n_spec_params - 1);
  for(int i = 0; i < n_spec_params; i++) {
    spectral_params[i] = spec_param_min + i*dspec_param;
  }
    

    
  std::vector< std::vector<double> > graph_embedding(n_pts);
  for(int i = 0; i < n_pts; i++) {
    graph_embedding[i] = spectral_embedding(alg->gen_unweighted_graph(graph_size, graph_params[i]), spectral_params);
  }
  std::cout << "--> Graphs generated and embedded" << std::endl;
  const int k = 5;
  std::vector<double> eigvals(k);
  std::vector< std::vector<double> > eigvects(k);
  std::vector< std::vector<double> > W(k);
  double median = get_median(get_squared_distances(graph_embedding));
  std::cout << "--> Median squared distance: " << median << std::endl;
  Gaussian_Kernel* gk = new Gaussian_Kernel(epsilon);
  const int dmaps_success = dmaps::map(graph_embedding, gk, eigvals, eigvects, W, k, 1e-12);
  if(dmaps_success != 1) {
    std::cout << "dmaps encountered an error" << std::endl;
  }
  std::cout << "--> DMAP computed" << std::endl;
  
  std::ofstream output_eigvals("./output_data/" + graph_type + "_embedding_eigvals.csv");
  std::ofstream output_eigvects("./output_data/" + graph_type + "_embedding_eigvects.csv");
  // std::ofstream output_W("./output_data/dmaps_W.csv");
  std::ofstream output_graphparams("./output_data/" + graph_type + "_embedding_graphparams.csv");

  save_vector(eigvals, output_eigvals);
  save_matrix(eigvects, output_eigvects);
  // save_matrix(W, output_W);
  save_matrix(graph_params, output_graphparams);
  output_eigvals.close();
  output_eigvects.close();
  // output_W.close();
  output_graphparams.close();

  std::cout << "<---------------------------------------->" << std::endl;

  // std::cout << "--> Testing" << std::endl;

  // std::vector<double> eigvals_EIG(n_pts);
  // std::vector< std::vector<double> > eigvects_EIG(n_pts);
  // std::vector< std::vector<double> > W_EIG(n_pts);
  // dmaps_EIG::map(graph_embedding, gk, eigvals_EIG, eigvects_EIG, W_EIG, 1e-12);

  // std::ofstream output_eigvals_EIG("./output_data/" + graph_type + "_embedding_eigvals_EIG.csv");
  // std::ofstream output_eigvects_EIG("./output_data/" + graph_type + "_embedding_eigvects_EIG.csv");
  // std::ofstream output_graphparams_EIG("./output_data/" + graph_type + "_embedding_graphparams_EIG.csv");

  // save_vector(eigvals_EIG, output_eigvals_EIG);
  // save_matrix(eigvects_EIG, output_eigvects_EIG);
  // save_matrix(graph_params, output_graphparams_EIG);
  // output_eigvals_EIG.close();
  // output_eigvects_EIG.close();
  // output_graphparams_EIG.close();

  // std::cout << "<---------------------------------------->" << std::endl;

  return 1;
}
