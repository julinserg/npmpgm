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
    int nClass = dataForSOM.n_rows;
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
    field<mat> ProbabilityA;
    field<mat> ProbabilityAt;
    ProbabilityA.set_size(nClass,1);
    ProbabilityAt.set_size(nClass,1);
    for(int k=0; k < data.n_rows; ++k)
    {
        mat A = data(k,0);
        rowvec arrayWinUnit;
        arrayWinUnit.set_size(A.n_rows);
        for(int i=0; i < A.n_rows; ++i)
        {
            rowvec dataVect = A(i,span::all);
            int winUnit = 0;
            double minDist = 10000000000000000;
            for (int j=0; j < CodeBook.n_rows;++j)
            {
               rowvec wVect =  CodeBook(j,span::all);
               rowvec distVec = dataVect - wVect;
               double dist = norm(distVec,2);
               if (dist < minDist)
               {
                   winUnit = j+1;
                   minDist = dist;
               }
            }
            arrayWinUnit(i) = winUnit;
        }
        mat B;
        mat C;
        int classL = label(k,0);
        B = ProbabilityA(classL-1,0);
        C = ProbabilityAt(classL-1,0);
        if(B.n_cols == 0)
        {
           B.set_size(nVectors,nVectors);
           B.zeros();
        }
        if (C.n_cols == 0)
        {
           C.set_size(nVectors,1);
           C.zeros();
        }
        for(int i = 0; i < arrayWinUnit.n_cols - 1; ++i)
        {
            int pp = arrayWinUnit(i);
            int qq = arrayWinUnit(i+1);
            B(pp,qq) = B(pp,qq) + 1;
            C(pp,0) = C(pp,0) + 1;
        }
        ProbabilityA(classL-1,0) = B;
        ProbabilityAt(classL-1,0) = C;
    }
    field<mat> MatrixTransactA;
    MatrixTransactA.set_size(nClass,1);
    for(int i=0; i < ProbabilityA.n_rows; ++i)
    {
        mat B = ProbabilityA(i,0);
        mat C = ProbabilityAt(i,0);
        mat A;
        A.copy_size(B);
        A.zeros();
        for(int j=0; j < B.n_rows; ++j)
        {
           rowvec vect;
           vect.set_size(B.n_cols);
           vect.zeros();
           if (C(j,0) != 0)
           {
              rowvec vect2 = B(j,span::all);
              int cc = C(j,0);
              vect =  vect2 / cc;
           }
           A(j,span::all) = vect;
        }
        MatrixTransactA(i,0) = A;

    }

    int hhh = 0;



}

MainWindow::~MainWindow()
{
    delete ui;
}
