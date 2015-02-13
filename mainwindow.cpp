#include "mainwindow.h"
#include "ui_mainwindow.h"

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
    writeSOMtrainfiles(C_SOMFILENAMEDATA);

    for(int i=1; i<=1/*m_nClass*/; ++i)
    {
         QString FileName = C_SOMFILENAMEDATA.arg(i).arg("csv");
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
         somthread->setObjectEvent(this,i);
         somthread->start();      
    }

}

MainWindow::~MainWindow()
{
    delete ui;
}

bool MainWindow::event(QEvent *ev)
{
    if(ev->type() == SOM_THREAD_STOP)
    {
        ESOMThreadStop* evend = (ESOMThreadStop*)ev;
        int label  =evend->labelClass;
        double * codebook = NULL;
        unsigned int nDimensions = 0;
        unsigned int nVectors = 0;
        QString filecodebook = C_SOMFILENAMEDATA.arg(label).arg("wts");
        codebook = readMatrixDoubleType(filecodebook.toStdString(), nVectors, nDimensions);
        mat CodeBook(codebook, nVectors, nDimensions);
        m_SOMCodeBookList(label-1,0) = CodeBook;
        mat TrasitionM = formMatrixTransaction(CodeBook,label);
        m_MatrixTransactA(label-1,0) = TrasitionM;
    }
}

void MainWindow::readTrainData(const QString &datafile, const QString &lablefile)
{
    CGetData::getCellFromFile(datafile,lablefile,m_TrainData,m_TrainLabel);
    CGetData::formingTrainDataForSOM(m_TrainData,m_TrainLabel,m_TrainDataForSOM);
    m_nClass = m_TrainDataForSOM.n_rows;
    m_SOMCodeBookList.set_size(m_nClass,1);
    m_MatrixTransactA.set_size(m_nClass,1);
}

void MainWindow::writeSOMtrainfiles(const QString &patternfilename)
{
    for(int k = 0; k < m_TrainDataForSOM.n_rows; ++k)
    {
        QString namefilecur = patternfilename.arg(k+1).arg("csv");
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

mat MainWindow::formMatrixTransaction(mat codebook, int label)
{
    mat ProbabilityA;
    mat ProbabilityAt;
    ProbabilityA.set_size(codebook.n_rows,codebook.n_rows);
    ProbabilityA.zeros();
    ProbabilityAt.set_size(codebook.n_rows,1);
    ProbabilityAt.zeros();
    for(int k=0; k < m_TrainData.n_rows; ++k)
    {
        if (m_TrainLabel(k,0) == label)
        {
            mat A = m_TrainData(k,0);
            rowvec arrayWinUnit;
            arrayWinUnit.set_size(A.n_rows);
            for(int i=0; i < A.n_rows; ++i)
            {
                rowvec dataVect = A(i,span::all);
                int winUnit = 0;
                double minDist = 10000000000000000;
                for (int j=0; j < codebook.n_rows;++j)
                {
                   rowvec wVect =  codebook(j,span::all);
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

            for(int i = 0; i < arrayWinUnit.n_cols - 1; ++i)
            {
                int pp = arrayWinUnit(i);
                int qq = arrayWinUnit(i+1);
                ProbabilityA(pp,qq) = ProbabilityA(pp,qq) + 1;
                ProbabilityAt(pp,0) = ProbabilityAt(pp,0) + 1;
            }
        }
    }

    mat A;
    A.copy_size(ProbabilityA);
    A.zeros();
    for(int j=0; j < ProbabilityA.n_rows; ++j)
    {
       rowvec vect;
       vect.set_size(ProbabilityA.n_cols);
       vect.zeros();
       if (ProbabilityAt(j,0) != 0)
       {
          rowvec vect2 = ProbabilityA(j,span::all);
          int cc = ProbabilityAt(j,0);
          vect =  vect2 / cc;
       }
       A(j,span::all) = vect;
    }

    return A;

}
