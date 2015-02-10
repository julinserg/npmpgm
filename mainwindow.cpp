#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "armadillo"
#include <stdlib.h>
#include "cgetdata.h"
#include <QFile>
#include <QTextStream>
using namespace arma;
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);  
    QString str1("data/characterTrainData.csv");
    QString str2("data/characterTrainLabel.csv");
    field<mat> data;
    mat label;
    CGetData::getCellFromFile(str1,str2,data,label);   
    field<mat> dataForSOM;
    CGetData::formingTrainDataForSOM(data,label,dataForSOM);
    QString namefile = "./som_train_class_%1.csv";
    for(int k = 0; k < dataForSOM.n_rows; ++k)
    {
        QString namefilecur = namefile.arg(k+1);
        QFile file(namefilecur);
        mat MATRIX = dataForSOM(k,0);
        if (file.open(QFile::WriteOnly|QFile::Truncate))
        {
          QTextStream stream(&file);
          for (int i = 0; i < MATRIX.n_rows; ++i)
          {
              for (int j = 0; j < MATRIX.n_cols; ++j)
              {
                  double val = MATRIX(i,j);
                  stream << val << " ";
              }
              stream << "\n";
          }
          file.close();
        }
    }
    int t =0;
    QString name = namefile.arg(1);
    m_somthread = new QSOMThread();

    m_somthread->setFileName(name.toStdString());
    m_somthread->setOutPrefix(QString("./som_train_class_%1").arg(1).toStdString());
    m_somthread->setNumEpoch(10);
    m_somthread->setSizeMap(12,12,"planar");
    m_somthread->setRadiusParam(12,1,"linear");
    m_somthread->setScaleParam(0.1,0.01,"linear");
    m_somthread->setKernelType(0,1);
    m_somthread->setSaveParam(0,"");
    m_somthread->start();
    m_somthread->moveToThread(m_somthread);
    m_somthread->startTrain();

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


}

MainWindow::~MainWindow()
{
    delete ui;
}
