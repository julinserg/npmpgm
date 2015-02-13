#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <stdlib.h>
#include "cgetdata.h"
#include <QFile>
#include <QTextStream>
using namespace arma;
const QString C_SOMFILENAMEDATA = "./som_train_class_%1.%2";
const QString C_SOMFILENAMECODEBOOK = "./som_train_class_%1";
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);  

    readTrainData("data/characterTrainData.csv","data/characterTrainLabel.csv");
    writeSOMtrainfiles(C_SOMFILENAME);

    for(int i=1; i<=m_nClass; ++i)
    {
         QString FileName = C_SOMFILENAMECODEBOOK.arg(i).arg("csv");
         QString OutPrefix = C_SOMFILENAMECODEBOOK.arg(i);
         QSOMThread*  somthread = new QSOMThread();
         int x = 15;
         int y = 15;
         somthread->setFileName(FileName.toStdString());
         somthread->setOutPrefix(OutPrefix.toStdString());
         somthread->setNumEpoch(500);
         somthread->setSizeMap(x,y,"planar");
         somthread->setRadiusParam(min(x,y),1,"linear");
         somthread->setScaleParam(0.1,0.01,"linear");
         somthread->setKernelType(0,4);
         somthread->setSaveParam(0,"");
         somthread->start();
         somthread->startTrain();
    }

    int t =0;
    QString name = namefile.arg(1);


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

void MainWindow::readTrainData(const QString &datafile, const QString &lablefile)
{
    CGetData::getCellFromFile(datafile,lablefile,m_TrainData,m_TrainLabel);
    CGetData::formingTrainDataForSOM(m_TrainData,m_TrainLabel,m_TrainDataForSOM);
    m_nClass = m_TrainDataForSOM.n_rows;
}

void MainWindow::writeSOMtrainfiles(const QString &patternfilename)
{
    for(int k = 0; k < m_TrainDataForSOM.n_rows; ++k)
    {
        QString namefilecur = patternfilename.arg(k+1).append("csv");
        QFile file(namefilecur);
        mat MATRIX = m_TrainDataForSOM(k,0);
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

}
