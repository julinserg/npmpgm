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
    /*m_somthread = new QSOMThread();
    int x = 15;
    int y = 15;
    m_somthread->setFileName(name.toStdString());
    m_somthread->setOutPrefix(QString("./som_train_class_%1").arg(1).toStdString());
    m_somthread->setNumEpoch(500);
    m_somthread->setSizeMap(x,y,"planar");
    m_somthread->setRadiusParam(min(x,y),1,"linear");
    m_somthread->setScaleParam(0.1,0.01,"linear");
    m_somthread->setKernelType(0,4);
    m_somthread->setSaveParam(0,"");
    m_somthread->start();
    m_somthread->moveToThread(m_somthread);
    m_somthread->startTrain();*/

    double * codebook = NULL;
    unsigned int nDimensions = 0;
    unsigned int nVectors = 0;
    QString filecodebook = "./som_train_class_%1.wts";
    filecodebook = filecodebook.arg(1);
    codebook = readMatrixDoubleType(filecodebook.toStdString(), nVectors, nDimensions);
    mat CodeBook(codebook, nVectors, nDimensions);
    int g = 0;





}

MainWindow::~MainWindow()
{
    delete ui;
}
