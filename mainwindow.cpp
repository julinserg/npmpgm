#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "csv.h"
#include "armadillo"
#include <stdlib.h>
using namespace arma;
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    m_somthread = new QSOMThread();
    m_somthread->start();
    m_somthread->moveToThread(m_somthread);
    mat B;
   /* mat A = randn(2,3);
    mat B = randn(4,5);

    field<mat> F(2,1);
    F(0,0) = A;
    F(1,0) = B;

    F.print("F:");

    F.save("mat_field");*/
   /* csv_parser csv("data/characterTrainLabel.csv");
    mat B;
    while(file_parser.has_more_rows())
    {
            unsigned int i = 0;
            csv_row row = file_parser.get_row();
            B.set_size(row_count,row.size());
            for (i = 0; i < row.size(); i++)
            {
                double val =  atof(row[i].c_str());

                B(row_count-1,i) = val;
            }
            row_count++;
    }*/
    QList<QStringList> listfile = CSV::parseFromFile("data/characterTrainLabel.csv");
    B.set_size(listfile.size(),listfile.first().size());
    int i = 0;
    int j = 0;
    foreach (QStringList rowList, listfile) {
        j = 0;
        foreach (QString valStr, rowList) {
             bool ok;
             B(i,j) = valStr.toDouble(&ok);
             j++;
        }
        i++;
    }

}

MainWindow::~MainWindow()
{
    delete ui;
}
