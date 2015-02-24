#include "mainwindow.h"
#include <QApplication>
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
int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    std::cout << "<---------------------------------------->" << std::endl;

    const int graph_size = atoi(argv[2]);
    std::string graph_type;
    Graph_Gen_Alg* alg;
    std::vector< std::vector<double> > graph_params;
    int n_pts;

    double epsilon = 1e-3;

    int argcL = 1;
    QString strComLine = "-np %1";
    strComLine = strComLine.arg(10);
    char *param = strComLine.toLatin1().data();
    char **argvL = &param;
    int m_thred = 4;
    int m_rank = 0;
   #ifdef HAVE_MPI
       ///
       /// MPI init
       ///
       MPI_Init(&argcL, &argvL);
       MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
       MPI_Comm_size(MPI_COMM_WORLD, &m_thred);
       MPI_Barrier(MPI_COMM_WORLD);
   #endif

    MainWindow w;
    w.show();

    return a.exec();
}
