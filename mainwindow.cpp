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


    m_nClassComplete = 0;
    //train();
    m_timertrain = new QTimer();
    bool k = connect(m_timertrain,SIGNAL(timeout()),this,SLOT(timeoutAnalysisTrainComplete()));
    m_timertrain->start(1000);
    test();
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
        m_nClassComplete++;
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

void MainWindow::readTestData(const QString &datafile, const QString &lablefile)
{
   CGetData::getCellFromFile(datafile,lablefile,m_TestData,m_TestLabel);
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
                       winUnit = j;
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

void MainWindow::timeoutAnalysisTrainComplete()
{
   /* if (m_nClass == m_nClassComplete)
    {
        // если все классы обучены, то выполняем сохранение в файл моделей
        bool f = m_SOMCodeBookList.save("weight_som.model");
        bool f1 = m_MatrixTransactA.save("matrix_a.model");
    }*/
    bool g1 = m_SOMCodeBookList.load("weight_som.model");
    bool g2 = m_MatrixTransactA.load("matrix_a.model");
    int nClass = m_SOMCodeBookList.n_rows;
    for(int i=0; i < m_TestData.n_rows; ++i)
    {
        mat SEQ = m_TestData(i,0);
        int SizeSEQ = SEQ.n_rows;
        for (int cl =0; cl < nClass; ++cl)
        {
            mat WSOM = m_SOMCodeBookList(cl,0);
            int SizeSOM = WSOM.n_rows;
            mat Z;
            Z.set_size(SizeSOM,SizeSEQ);
            Z.zeros();
            for (uint t=0; t < SizeSOM; ++t)
            {
                rowvec R0 = WSOM(t,span::all);
                mat R1;
                R1.set_size(SizeSEQ,R0.n_elem);
                R1.zeros();
                R1.each_row() += R0;
                mat R2 = R1-SEQ;
                mat R3 = square(R2);
                mat R4 = sum(R3,1);
                rowvec R5 = R4.t();
                Z(t,span::all) = R5;
            }
            Z = -sqrt(Z);
            mat L = logsumexp(Z.t(),2);
            mat LMAT = repmat(L, 1, Z.n_rows);
            mat Znorm = Z - LMAT.t();
            mat B = exp(Znorm);
            mat G;
            G.set_size(1,1);
            G(0,0) = 0;
            mat PI = repmat(G,1,SizeSOM);
            PI(0,0) = 1;
            mat A = m_MatrixTransactA(cl,0);
            double logp = hmmFilter(PI,A,B);
        }
    }

    int g = 0;
}

void MainWindow::train()
{
    readTrainData("data/characterTrainData.csv","data/characterTrainLabel.csv");
    writeSOMtrainfiles(C_SOMFILENAMEDATA);
    for(int i=1; i<=m_nClass; ++i)
    {
         QString FileName = C_SOMFILENAMEDATA.arg(i).arg("csv");
         QString OutPrefix = C_SOMFILENAMECODEBOOK.arg(i);
         QSOMThread*  somthread = new QSOMThread();
         int x = 10;
         int y = 10;
         somthread->setFileName(FileName.toStdString());
         somthread->setOutPrefix(OutPrefix.toStdString());
         somthread->setNumEpoch(50);
         somthread->setSizeMap(x,y,"planar");
         somthread->setRadiusParam(min(x,y),1,"linear");
         somthread->setScaleParam(0.1,0.01,"linear");
         somthread->setKernelType(0,4);
         somthread->setSaveParam(0,"");
         somthread->setObjectEvent(this,i);
         somthread->start();
    }
}

void MainWindow::test()
{
     readTestData("data/characterTestData.csv","data/characterTestLabel.csv");
}

mat MainWindow::logsumexp(mat a, int dim)
{
    mat y = max(a,dim-1);
    mat dims = ones(1,2);
    if (dim == 1)
    {
        dims(0,dim-1) = a.n_rows;
    }
    if (dim == 2)
    {
        dims(0,dim-1) = a.n_cols;
    }
    mat b = repmat(y,dims(0,0),dims(0,1));
    a = a - b;
    mat s = y + log(sum(exp(a),dim-1));
    return s;
}

double MainWindow::hmmFilter(mat initDist, mat transmat, mat softev)
{
  int K =  softev.n_rows;
  int T = softev.n_cols;

  colvec scale;
  scale.set_size(T,1);
  scale.zeros();
  mat AT = transmat.t();
  mat alpha;
  mat R1 = vectorise(initDist) % softev(span::all,0);
  double z = sum( vectorise(R1));
  if (z == 0)
  {
      z = 1;
  }
 // Z.elem( find(Z ==0) ).ones();
  R1 = R1 / z;
  alpha = R1;
  scale(0) = z;
  for(int t =1; t < T; ++t)
  {
      mat R2 = (AT * alpha) %  softev(span::all,t);
      double z = sum( vectorise(R2));
      if (z == 0)
      {
          z = 1;
      }
      R2 = R2 / z;
      alpha = R2;
      scale(t) = z;
  }
  scale.print("scale:");
  alpha.print("alpha:");
  double loglike = sum(log(scale+datum::eps));
  return loglike;

}
