#include "mainwindow.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
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
