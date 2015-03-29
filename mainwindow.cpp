#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "util_fns.h"
#include "embeddings.h"
#include "gaussian_kernel.h"
#include "dmaps.h"
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
    train();
    m_timertrain = new QTimer();
    bool k = connect(m_timertrain,SIGNAL(timeout()),this,SLOT(timeoutAnalysisTrainComplete()));
    m_timertrain->start(1000);
   // test();
   // test2();

     // testEnsemble();

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
        mat CodeBook;
        CodeBook.set_size(nVectors,nDimensions);
        int index = 0;
        for(int i = 0; i<nVectors;++i)
        {
            for(int j =0; j < nDimensions;++j)
            {
                CodeBook(i,j) = codebook[index];
                ++index;
            }
        }
        mat UMATRIX =  calculateDistNodeMatrix(CodeBook);
        m_UmatrixGraphList(label-1,0) = UMATRIX;
        m_SOMCodeBookList(label-1,0) = CodeBook;
        mat TrasitionM;
        mat PiM;
        formMatrixTransaction(CodeBook,label,TrasitionM,PiM);
        m_MatrixTransactA(label-1,0) = TrasitionM;
        m_MatrixPI(label-1,0) = PiM;
        m_nClassComplete++;
    }
}

void MainWindow::readTrainData(const QString &datafile)
{
    CGetData::getCellFromFile(datafile,m_TrainData,m_TrainLabel);
    CGetData::formingTrainDataForSOM(m_TrainData,m_TrainLabel,m_TrainDataForSOM);
    m_nClass = m_TrainDataForSOM.n_rows;
    m_SOMCodeBookList.set_size(m_nClass,1);
    m_MatrixTransactA.set_size(m_nClass,1);
    m_UmatrixGraphList.set_size(m_nClass,1);
    m_MatrixPI.set_size(m_nClass,1);
}

void MainWindow::readTestData(const QString &datafile)
{
   CGetData::getCellFromFile(datafile,m_TestData,m_TestLabel);
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

void MainWindow::formMatrixTransaction(mat codebook, int label, mat& TR, mat& PI)
{
    mat ProbabilityA;
    mat ProbabilityAt;

    ProbabilityA.set_size(codebook.n_rows,codebook.n_rows);
    ProbabilityA.zeros();
    ProbabilityAt.set_size(codebook.n_rows,1);
    ProbabilityAt.zeros();
    PI.set_size(1,codebook.n_rows);
    PI.zeros();
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
                PI(0,winUnit) = PI(0,winUnit) + 1;
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


    TR.copy_size(ProbabilityA);
    TR.zeros();
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
       TR(j,span::all) = vect;
    }
   // PI.print("PI:");
    double sumPi = 0;
    sumPi = sum(sum(PI));
    if (sumPi == 0)
    {
        sumPi = 1;
    }
    PI = PI / sumPi;
   /* mat LPI = logsumexp(PI,2);
    mat LMATPI = repmat(LPI, 1, PI.n_cols);
  //  LMATPI.print("LMATPI:");
    mat PInorm = PI - LMATPI;
    PI = exp(PInorm);*/
   // PI.print("PI:");
   // int g = 0;

}

void MainWindow::timeoutAnalysisTrainComplete()
{
    if (m_nClass == m_nClassComplete)
    {
        // если все классы обучены, то выполняем сохранение в файл моделей
        bool f = m_SOMCodeBookList.save("weight_som.model");
        bool f1 = m_MatrixTransactA.save("matrix_a.model");
        bool f2 = m_UmatrixGraphList.save("matrix_u.model");
        bool f3 = m_MatrixPI.save("matrix_pi.model");
        m_timertrain->stop();
        test();
       // test2();
       // testEnsemble();
    }

}

void MainWindow::train()
{
    readTrainData("data/characterTrainData.csv");
    writeSOMtrainfiles(C_SOMFILENAMEDATA);
    for(int i=1; i<=m_nClass; ++i)
    {
         QString FileName = C_SOMFILENAMEDATA.arg(i).arg("csv");
         QString OutPrefix = C_SOMFILENAMECODEBOOK.arg(i);
         QSOMThread*  somthread = new QSOMThread();
         m_nSOM_X = 15;
         m_nSOM_Y = 15;
         m_mapType = "planar";
         somthread->setFileName(FileName.toStdString());
         somthread->setOutPrefix(OutPrefix.toStdString());
         somthread->setNumEpoch(1000);
         somthread->setSizeMap(m_nSOM_X,m_nSOM_Y,m_mapType);
         somthread->setRadiusParam(3,1,"linear");
         somthread->setScaleParam(0.1,0.01,"linear");
         somthread->setKernelType(0,4);
         somthread->setSaveParam(0,"");
         somthread->setObjectEvent(this,i);
         somthread->start();
    }
}

void MainWindow::test()
{
    readTestData("data/characterTestData.csv");
    bool g1 = m_SOMCodeBookList.load("weight_som.model");
    bool g2 = m_MatrixTransactA.load("matrix_a.model");
    bool g3 = m_MatrixPI.load("matrix_pi.model");
    int nClass = m_SOMCodeBookList.n_rows;
    mat arrayLL;
    arrayLL.set_size(m_TestData.n_rows,nClass);
    for(int i=0; i < m_TestData.n_rows; ++i)
    {
      mat SEQ = m_TestData(i,0);
     // SEQ.print("SEQ");
      int SizeSEQ = SEQ.n_rows;
      for (int cl =0; cl < nClass; ++cl)
      {
          mat WSOM = m_SOMCodeBookList(cl,0);
         // WSOM.print("SOM");
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
        //  L.print("L");
          mat LMAT = repmat(L, 1, Z.n_rows);
          mat Znorm = Z - LMAT.t();
          mat B = exp(Znorm);
         /* mat G;
          G.set_size(1,1);
          G(0,0) = 5;
          mat PI = repmat(G,1,SizeSOM);
          colvec vecmax = Z(span::all,0);
          uword maxindex;
          double pimax = vecmax.max(maxindex);
          PI(0,maxindex) = 10;

          mat LPI = logsumexp(PI,2);
          mat LMATPI = repmat(LPI, 1, PI.n_cols);

          mat PInorm = PI - LMATPI;
          PI = exp(PInorm);*/

          mat A = m_MatrixTransactA(cl,0);
          mat PI = m_MatrixPI(cl,0);
          double logp = hmmFilter(PI,A,B);
          logp = logp + sum(vectorise(L));
          arrayLL(i,cl) =logp;
      }
    }
    rowvec labelDetect;
    labelDetect.set_size(arrayLL.n_rows);
    for(int i = 0; i < arrayLL.n_rows; ++i)
    {
      rowvec vec = arrayLL(i,span::all);
      uword  index;
      double max = vec.max(index);
      labelDetect(i) = index+1;
    }
    rowvec labelTrue = m_TestLabel.t();
    //labelTrue.resize(200);
    double accuracy = 0;
    for(int i=0; i < labelTrue.n_elem; ++i)
    {
       if (labelTrue(i) == labelDetect(i))
       {
           accuracy++;
       }
    }
    accuracy = accuracy / labelTrue.n_elem;
    qDebug("accuracy: %f",accuracy);
    double fmeasure;
    double precision;
    double recall;
    quality(labelDetect,labelTrue,nClass,fmeasure,precision,recall);
    qDebug("fmeasure: %f",fmeasure);
    qDebug("precision: %f",precision);
    qDebug("recall: %f",recall);

    mat DP = mapminmax(arrayLL,0,1);
    DP.save("classif_1.result");  
    int gg = 0;
}

void MainWindow::test2()
{
    const int n_spec_params = 10;
    std::vector<double> spectral_params(n_spec_params);
    const double spec_param_max = 0.01;
    const double spec_param_min = 0.0001;
    const double dspec_param = (spec_param_max - spec_param_min)/(n_spec_params - 1);
    for(int i = 0; i < n_spec_params; i++) {
      spectral_params[i] = spec_param_min + i*dspec_param;
    }


    readTestData("data/characterTestData.csv");
    bool g1 = m_SOMCodeBookList.load("weight_som.model");
    bool g2 = m_UmatrixGraphList.load("matrix_u.model");
    int nClass = m_SOMCodeBookList.n_rows;
    std::vector< std::vector<double> > specmodelvec(nClass, std::vector<double>(n_spec_params));
    for(int i=0;i<nClass;++i)
    {
       mat UMAT = m_UmatrixGraphList(i,0);
       Eigen::MatrixXd g_model =  convertArmadilloToEngineMatrix(UMAT);
       std::vector<double> vec = spectral_embedding(g_model, spectral_params);
       specmodelvec[i] = vec;
    }
    mat arrayLL;
    arrayLL.set_size(m_TestData.n_rows,nClass);
    for(int i=0; i < m_TestData.n_rows; ++i)
    {
      mat SEQ = m_TestData(i,0);
     // SEQ.print("SEQ");
      int SizeSEQ = SEQ.n_rows;
      for (int cl =0; cl < nClass; ++cl)
      {
          mat WSOM = m_SOMCodeBookList(cl,0);
          mat UMAT = m_UmatrixGraphList(cl,0);
         // WSOM.print("SOM");
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
          uvec maxvec;
          maxvec.set_size(Z.n_cols);
          for(int i=0; i < Z.n_cols; ++i )
          {
              colvec vecmax = Z(span::all,i);
              uword maxindex;
              double max = vecmax.max(maxindex);
              maxvec(i) = maxindex;
          }
         // UMAT.print("UMAT:");
          mat GraphData;
          GraphData.set_size(UMAT.n_rows,UMAT.n_cols);
          GraphData.zeros();
          for (int i=0;i< UMAT.n_rows;++i)
          {
              for (int j=0;j< UMAT.n_cols;++j)
              {
                  uvec z1 = maxvec.elem(arma::find(maxvec == i));
                  uvec z2 = maxvec.elem(arma::find(maxvec == j));
                  if (z1.n_elem > 0 || z2.n_elem > 0)
                  {
                      GraphData(i,j) = UMAT(i,j);
                  }
              }
          }
         //GraphData.print("GraphData:");
          std::vector< std::vector<double> > graph_embedding(2);

          Eigen::MatrixXd g_data =  convertArmadilloToEngineMatrix(GraphData);
          std::vector<double> vecModel(n_spec_params);
          vecModel = specmodelvec[cl];

          graph_embedding[0] = vecModel;
          graph_embedding[1] = spectral_embedding(g_data, spectral_params);
          const int k = 2;
          std::vector<double> eigvals_EIG(k);
          std::vector< std::vector<double> > eigvects_EIG(k);
          std::vector< std::vector<double> > W_EIG(k);
          double epsilon = 1e-3;
          Gaussian_Kernel* gk = new Gaussian_Kernel(epsilon);
          const int dmaps_success = dmaps::map(graph_embedding, gk, eigvals_EIG, eigvects_EIG, W_EIG,k, 1e-12);
          //double l2n = get_l2norm(get_squared_distances(graph_embedding));
         // arrayLL(i,cl) = -l2n;
          arrayLL(i,cl) = W_EIG[0][1];
      }
    }
    rowvec labelDetect;
    labelDetect.set_size(arrayLL.n_rows);
    for(int i = 0; i < arrayLL.n_rows; ++i)
    {
      rowvec vec = arrayLL(i,span::all);
      uword  index;
      double max = vec.max(index);
      labelDetect(i) = index+1;
    }
    labelDetect.print("labelDetect:");
    rowvec labelTrue = m_TestLabel.t();
    //labelTrue.resize(200);
    double accuracy = 0;
    for(int i=0; i < labelTrue.n_elem; ++i)
    {
       if (labelTrue(i) == labelDetect(i))
       {
           accuracy++;
       }
    }
    accuracy = accuracy / labelTrue.n_elem;
    qDebug("accuracy: %f",accuracy);
    double fmeasure;
    double precision;
    double recall;
    quality(labelDetect,labelTrue,nClass,fmeasure,precision,recall);
    qDebug("fmeasure: %f",fmeasure);
    qDebug("precision: %f",precision);
    qDebug("recall: %f",recall);

    mat DP = mapminmax(arrayLL,0,1);
    DP.save("classif_2.result");

}

void MainWindow::testEnsemble()
{
    readTestData("data/characterTestData.csv");
    mat DP1;
    mat DP2;
    bool g1 = DP1.load("classif_1.result");
    bool g2 = DP2.load("classif_2.result");
    int nClass = DP1.n_cols;
    mat arrayLL_MAX;
    arrayLL_MAX.set_size(DP1.n_rows,DP1.n_cols);
    mat arrayLL_MIN;
    arrayLL_MIN.set_size(DP1.n_rows,DP1.n_cols);
    mat arrayLL_SUM;
    arrayLL_SUM.set_size(DP1.n_rows,DP1.n_cols);
    mat arrayLL_AVR;
    arrayLL_AVR.set_size(DP1.n_rows,DP1.n_cols);
    mat arrayLL_PRO;
    arrayLL_PRO.set_size(DP1.n_rows,DP1.n_cols);
    for(int i=0; i < DP1.n_rows; ++i)
    {
        for(int j=0; j < DP1.n_cols; ++j)
        {
            rowvec vec;
            vec.set_size(2);
            vec(0) = DP1(i,j);
            vec(1) = DP2(i,j);
            arrayLL_MAX(i,j) = max(vec);
            arrayLL_MIN(i,j) = min(vec);
            arrayLL_SUM(i,j) = sum(vec);
            arrayLL_AVR(i,j) = mean(vec);
            arrayLL_PRO(i,j) = prod(vec);
        }
    }

    rowvec labelTrue = m_TestLabel.t();


    rowvec labelDetect_MAX = calcLabelDetect(arrayLL_MAX);
    rowvec labelDetect_MIN = calcLabelDetect(arrayLL_MIN);
    rowvec labelDetect_SUM = calcLabelDetect(arrayLL_SUM);
    rowvec labelDetect_AVR = calcLabelDetect(arrayLL_AVR);
    rowvec labelDetect_PRO = calcLabelDetect(arrayLL_PRO);
    labelDetect_MAX.print("labelDetect_MAX:");
    //labelTrue.resize(200);
    double accuracy_max = 0;
    for(int i=0; i < labelTrue.n_elem; ++i)
    {
       if (labelTrue(i) == labelDetect_MAX(i))
       {
           accuracy_max++;
       }
    }
    accuracy_max = accuracy_max / labelTrue.n_elem;
    qDebug("accuracy_max: %f",accuracy_max);

    double accuracy_min = 0;
    for(int i=0; i < labelTrue.n_elem; ++i)
    {
       if (labelTrue(i) == labelDetect_MIN(i))
       {
           accuracy_min++;
       }
    }
    accuracy_min = accuracy_min / labelTrue.n_elem;
    qDebug("accuracy_min: %f",accuracy_min);

    double accuracy_sum = 0;
    for(int i=0; i < labelTrue.n_elem; ++i)
    {
       if (labelTrue(i) == labelDetect_SUM(i))
       {
           accuracy_sum++;
       }
    }
    accuracy_sum = accuracy_sum / labelTrue.n_elem;
    qDebug("accuracy_sum: %f",accuracy_sum);

    double accuracy_avr = 0;
    for(int i=0; i < labelTrue.n_elem; ++i)
    {
       if (labelTrue(i) == labelDetect_AVR(i))
       {
           accuracy_avr++;
       }
    }
    accuracy_avr = accuracy_avr / labelTrue.n_elem;
    qDebug("accuracy_avr: %f",accuracy_avr);

    double accuracy_pro = 0;
    for(int i=0; i < labelTrue.n_elem; ++i)
    {
       if (labelTrue(i) == labelDetect_PRO(i))
       {
           accuracy_pro++;
       }
    }
    accuracy_pro = accuracy_pro / labelTrue.n_elem;
    qDebug("accuracy_pro: %f",accuracy_pro);


    double fmeasure_max;
    double precision_max;
    double recall_max;
    quality(labelDetect_MAX,labelTrue,nClass,fmeasure_max,precision_max,recall_max);
    qDebug("fmeasure_max: %f",fmeasure_max);
    qDebug("precision_max: %f",precision_max);
    qDebug("recall_max: %f",recall_max);

    double fmeasure_min;
    double precision_min;
    double recall_min;
    quality(labelDetect_MIN,labelTrue,nClass,fmeasure_min,precision_min,recall_min);
    qDebug("fmeasure_min: %f",fmeasure_min);
    qDebug("precision_min: %f",precision_min);
    qDebug("recall_min: %f",recall_min);

    double fmeasure_sum;
    double precision_sum;
    double recall_sum;
    quality(labelDetect_SUM,labelTrue,nClass,fmeasure_sum,precision_sum,recall_sum);
    qDebug("fmeasure_sum: %f",fmeasure_sum);
    qDebug("precision_sum: %f",precision_sum);
    qDebug("recall_sum: %f",recall_sum);

    double fmeasure_avr;
    double precision_avr;
    double recall_avr;
    quality(labelDetect_AVR,labelTrue,nClass,fmeasure_avr,precision_avr,recall_avr);
    qDebug("fmeasure_avr: %f",fmeasure_avr);
    qDebug("precision_avr: %f",precision_avr);
    qDebug("recall_avr: %f",recall_avr);

    double fmeasure_prod;
    double precision_prod;
    double recall_prod;
    quality(labelDetect_PRO,labelTrue,nClass,fmeasure_prod,precision_prod,recall_prod);
    qDebug("fmeasure_prod: %f",fmeasure_prod);
    qDebug("precision_prod: %f",precision_prod);
    qDebug("recall_prod: %f",recall_prod);


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
  double loglike = sum(log(scale+datum::eps));
  return loglike;

}

void MainWindow::quality(rowvec labeldetect, rowvec labeltrue, int nClass, double& fmesure, double& precision, double& recall)
{
    mat Confusion_Matrix;
    Confusion_Matrix.set_size(nClass,nClass);
    Confusion_Matrix.zeros();
   // labeldetect.print("labeldetect");
   // labeltrue.print("labeltrue");
    if (labeldetect.n_elem == labeltrue.n_elem )
    {
        for(int i = 0; i < labeldetect.n_elem ; ++i)
        {
            int column =  labeltrue(i) - 1;
            int row = labeldetect(i) - 1;
            Confusion_Matrix(row,column) = Confusion_Matrix(row,column) + 1;
        }
    }
    rowvec PrecisionVec;
    PrecisionVec.set_size(nClass);
    rowvec RecallVec;
    RecallVec.set_size(nClass);
    for (int i= 0; i < nClass; ++i)
    {
        rowvec p = Confusion_Matrix(i,span::all);
        if (sum(p) == 0)
        {
           PrecisionVec(i) = 0;
        }
        else
        {
           PrecisionVec(i) =  Confusion_Matrix(i,i) / sum(p);
        }
        colvec r = Confusion_Matrix(span::all,i);
        if (sum(r) == 0)
        {
           RecallVec(i) = 0;
        }
        else
        {
           RecallVec(i) =  Confusion_Matrix(i,i) / sum(r);
        }
    }
    precision = sum(PrecisionVec) / nClass;
    recall = sum(RecallVec) / nClass;
    fmesure = 2*precision*recall / (precision + recall);
    return ;
}

Eigen::MatrixXd MainWindow::convertArmadilloToEngineMatrix(mat matrix)
{
    Eigen::MatrixXd G_eigen(matrix.n_rows,matrix.n_cols);
    for(int i = 0; i < matrix.n_rows; i++) {
      for(int j = 0; j < matrix.n_cols; j++) {
        G_eigen(i,j) = matrix(i,j);
      }
    }
    return G_eigen;
}

mat MainWindow::calculateDistNodeMatrix(mat codebook)
{
    mat Dist;
    int numNod = codebook.n_rows;
    Dist.set_size(numNod,numNod);
    Dist.zeros();
    for(int i=0; i < numNod; ++i)
    {
        int som_x1 = 0;
        int som_y1 = 0;
        twoFromOne((ulong)i,(ushort)m_nSOM_Y,(ushort&)som_x1,(ushort&)som_y1);
        for(int j=0; j < numNod; ++j)
        {
            int som_x2 = 0;
            int som_y2 = 0;
            twoFromOne((ulong)j,(ushort)m_nSOM_Y,(ushort&)som_x2,(ushort&)som_y2);
            float tmp = 0.0f;
            if (m_mapType == "planar")
            {
                tmp = euclideanDistanceOnPlanarMap(som_x1, som_y1, som_x2, som_y2);
            }
            if (m_mapType == "toroid")
            {
               tmp = euclideanDistanceOnToroidMap(som_x1, som_y1, som_x2, som_y2, m_nSOM_X, m_nSOM_Y);
            }
            if (tmp <= 1.001f)
            {
                Dist(i,j) = norm(codebook(i,span::all) - codebook(j,span::all),2);
            }
        }
    }
    colvec A = Dist.elem(arma::find(Dist != 0));
    float maxV = max(A);
    float minV = min(A);
    for(int i=0; i < numNod; ++i)
    {
        for(int j=0; j < numNod; ++j)
        {
           if(Dist(i,j) != 0)
           {
              Dist(i,j) = 1 - ((Dist(i,j) -  minV) / (maxV - minV));
           }
        }
    }
    //Dist.print("Dist:");
    return Dist;
}

double MainWindow::get_l2norm(std::vector<double> vec)
{
    double res = 0;
    for (int i= 0; i < vec.size(); ++i)
    {
        res = res + pow(vec[i],2);
    }
    res = sqrt(res);
    return res;
}

void MainWindow::twoFromOne (ulong z, ushort max_y, ushort& x, ushort& y)
{
    if ( z % max_y == 0)
    {
        x = z/max_y ;
    }
    else
    {
        x = z/max_y + 1;
    }
    y = z - max_y * (x-1);

}

mat MainWindow::mapminmax(mat matrix, double ymin, double ymax)
{
    mat resMatrix;
    resMatrix.set_size(matrix.n_rows,matrix.n_cols);
    for(int i=0; i < matrix.n_rows; ++i)
    {
        rowvec vec = matrix(i,span::all);
        double xmax = vec.max();
        double xmin = vec.min();
        for(int j=0;j < vec.n_elem; ++j)
        {
            resMatrix(i,j) = (ymax-ymin)*(vec(j)-xmin)/(xmax-xmin) + ymin;
        }
    }
    return resMatrix;
}

rowvec MainWindow::calcLabelDetect(mat arrayLL)
{
    rowvec labelDetect;
    labelDetect.set_size(arrayLL.n_rows);
    for(int i = 0; i < arrayLL.n_rows; ++i)
    {
      rowvec vec = arrayLL(i,span::all);
      uword  index;
      double max = vec.max(index);
      labelDetect(i) = index+1;
    }
    return labelDetect;
}







