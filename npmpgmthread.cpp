#include "npmpgmthread.h"
#include "util_fns.h"
#include "embeddings.h"
#include "gaussian_kernel.h"
#include "dmaps.h"
#include <QFile>
#include <QTextStream>
#include<QDebug>
using namespace arma;
const QString C_SOMFILENAMEDATA = "som_train_class_%1.%2";
const QString C_SOMFILENAMECODEBOOK = "som_train_class_%1";
class_npmpgm_model::class_npmpgm_model()
{
    m_nclassComplete = 0;   
    m_pathtomodel = "";
}

class_npmpgm_model::~class_npmpgm_model()
{
    delete m_timertrain;
}

void class_npmpgm_model::run()
{
    m_timertrain = new QTimer();
    bool k = connect(m_timertrain,SIGNAL(timeout()),this,SLOT(timeout_analysis_train_complete()));
    exec();
}

bool class_npmpgm_model::event(QEvent *ev)
{
    if(ev->type() == SOM_THREAD_STOP)
    {
        ESOMThreadStop* evend = (ESOMThreadStop*)ev;
        int label  =evend->labelClass;
        double * codebook = NULL;
        unsigned int nDimensions = 0;
        unsigned int nVectors = 0;
        QString filecodebook = m_pathtomodel + C_SOMFILENAMEDATA.arg(label).arg("wts");
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
        mat UMATRIX =  calculate_dist_node_matrix(CodeBook);
        m_umatrix_graph_list(label-1,0) = UMATRIX;
        m_som_codebook_list(label-1,0) = CodeBook;
        mat TrasitionM;
        mat PiM;
        form_matrix_transaction(CodeBook,label,TrasitionM,PiM);
        m_matrix_transact_a(label-1,0) = TrasitionM;
        m_matrix_pi(label-1,0) = PiM;
        m_nclassComplete++;
        emit number_class_complet(m_nclassComplete);
        return true;
    }
    return QThread::event(ev);
}

void class_npmpgm_model::read_train_data(QString datafile)
{
    if (CGetData::get_cell_from_file(datafile,m_train_data,m_train_label))
    {
        CGetData::forming_train_data_for_som(m_train_data,m_train_label,m_train_data_for_som);
        m_nclass = m_train_data_for_som.n_rows;
        m_som_codebook_list.set_size(m_nclass,1);
        m_matrix_transact_a.set_size(m_nclass,1);
        m_umatrix_graph_list.set_size(m_nclass,1);
        m_matrix_pi.set_size(m_nclass,1);

        emit end_load_traindata(true);
    }
    else
    {
        qDebug() << "No train data!";
        emit end_load_traindata(false);
    }
}

void class_npmpgm_model::read_test_data(QString datafile)
{
    if (CGetData::get_cell_from_file(datafile,m_test_data,m_test_label))
    {
        emit end_load_testdata(true);
    }
    else
    {
        qDebug() << "No test data!";
        emit end_load_testdata(false);
    }
}

void class_npmpgm_model::set_path_to_model(QString path)
{
    m_pathtomodel = path;
}

void class_npmpgm_model::read_test_model(QString path)
{
    bool g1 = m_som_codebook_list.load(path.toStdString()+"weight_som.model");
    bool g2 = m_matrix_transact_a.load(path.toStdString()+"matrix_a.model");
    bool g3 = m_matrix_pi.load(path.toStdString()+"matrix_pi.model");
    if (g1 && g2 && g3)
    {
        emit end_load_testmodel(true);
    }
    else
    {
        qDebug() << "No model!";
        emit end_load_testmodel(false);
    }
}

void class_npmpgm_model::write_som_trainfiles(const QString &patternfilename)
{
    for(int k = 0; k < m_train_data_for_som.n_rows; ++k)
    {
        QString namefilecur = patternfilename.arg(k+1).arg("csv");
        QFile file(namefilecur);
        mat MATRIX = m_train_data_for_som(k,0);
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

void class_npmpgm_model::form_matrix_transaction(mat codebook, int label, mat& TR, mat& PI)
{
    mat ProbabilityA;
    mat ProbabilityAt;

    ProbabilityA.set_size(codebook.n_rows,codebook.n_rows);
    ProbabilityA.zeros();
    ProbabilityAt.set_size(codebook.n_rows,1);
    ProbabilityAt.zeros();
    PI.set_size(1,codebook.n_rows);
    PI.zeros();
    for(int k=0; k < m_train_data.n_rows; ++k)
    {
        if (m_train_label(k,0) == label)
        {
            mat A = m_train_data(k,0);
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

void class_npmpgm_model::timeout_analysis_train_complete()
{
    if (m_nclass == m_nclassComplete)
    {
        // если все классы обучены, то выполняем сохранение в файл моделей
        bool f = m_som_codebook_list.save(m_pathtomodel.toStdString()+"weight_som.model");
        bool f1 = m_matrix_transact_a.save(m_pathtomodel.toStdString()+"matrix_a.model");
        bool f2 = m_umatrix_graph_list.save(m_pathtomodel.toStdString()+"matrix_u.model");
        bool f3 = m_matrix_pi.save(m_pathtomodel.toStdString()+"matrix_pi.model");
        m_timertrain->stop();
    }
}

void class_npmpgm_model::train(int numEpoch,int nSOM_X, int nSOM_Y,
                         QString mapType, int bRadius, int eRadiusm,
                         QString typeRadius, int bScale, int eScale, QString typeScale)
{
    m_nclassComplete = 0;
    m_timertrain->start(1000);
    write_som_trainfiles(m_pathtomodel + C_SOMFILENAMEDATA);
    for(int i=1; i<=m_nclass; ++i)
    {
         QString FileName = m_pathtomodel + C_SOMFILENAMEDATA.arg(i).arg("csv");
         QString OutPrefix = m_pathtomodel + C_SOMFILENAMECODEBOOK.arg(i);
         class_som_thread*  somthread = new class_som_thread();
         m_nsom_x = nSOM_X;
         m_nsom_y = nSOM_Y;
         m_map_type = mapType.toStdString();
         somthread->setFileName(FileName.toStdString());
         somthread->setOutPrefix(OutPrefix.toStdString());
         somthread->setNumEpoch(numEpoch);
         somthread->setSizeMap(m_nsom_x,m_nsom_y,m_map_type);
         somthread->setRadiusParam(bRadius,eRadiusm,typeRadius.toStdString());
         somthread->setScaleParam(bScale,eScale,typeScale.toStdString());
         somthread->setKernelType(0,4);
         somthread->setSaveParam(0,"");
         somthread->setObjectEvent(this,i);
         somthread->start();
    }
}

void class_npmpgm_model::test()
{
    int nClass = m_som_codebook_list.n_rows;
    if (nClass <= 0)
    {
        return ;
    }
    mat arrayLL;
    arrayLL.set_size(m_test_data.n_rows,nClass);
    for(int i=0; i < m_test_data.n_rows; ++i)
    {
      mat SEQ = m_test_data(i,0);      
      int SizeSEQ = SEQ.n_rows;
      for (int cl =0; cl < nClass; ++cl)
      {
          mat WSOM = m_som_codebook_list(cl,0);
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

          mat A = m_matrix_transact_a(cl,0);
          mat PI = m_matrix_pi(cl,0);
          double logp = hmm_filter(PI,A,B);
          logp = logp + sum(vectorise(L));
          arrayLL(i,cl) =logp;
      }
      emit number_testdata_complet(i);
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
    rowvec labelTrue = m_test_label.t();
    //labelTrue.resize(200);
    double accuracy = 0;
    for(int i=0; i < labelTrue.n_elem; ++i)
    {
       if (labelTrue(i) == labelDetect(i))
       {
           accuracy++;
       }
    }
    labelTrue.print("labelTrue:");
    labelDetect.print("labelDetect:");
    accuracy = accuracy / labelTrue.n_elem;
    qDebug("accuracy: %f",accuracy);
    double fmeasure;
    double precision;
    double recall;
    quality(labelDetect,labelTrue,nClass,fmeasure,precision,recall);
    qDebug("fmeasure: %f",fmeasure);
    qDebug("precision: %f",precision);
    qDebug("recall: %f",recall);
    emit result_testing(fmeasure,precision,recall);

    ///!!!!!!!!! РАСКОМЕНТИРОВАТЬ ДЛЯ АНСАМБЛЯ
   /* mat DP = mapminmax(arrayLL,0,1);
    DP.save("classif_1.result"); */
}

void class_npmpgm_model::test2()
{
    const int n_spec_params = 10;
    std::vector<double> spectral_params(n_spec_params);
    const double spec_param_max = 0.01;
    const double spec_param_min = 0.0001;
    const double dspec_param = (spec_param_max - spec_param_min)/(n_spec_params - 1);
    for(int i = 0; i < n_spec_params; i++) {
      spectral_params[i] = spec_param_min + i*dspec_param;
    }


    read_test_data("data/characterTestData.csv");
    bool g1 = m_som_codebook_list.load("weight_som.model");
    bool g2 = m_umatrix_graph_list.load("matrix_u.model");
    int nClass = m_som_codebook_list.n_rows;
    std::vector< std::vector<double> > specmodelvec(nClass, std::vector<double>(n_spec_params));
    for(int i=0;i<nClass;++i)
    {
       mat UMAT = m_umatrix_graph_list(i,0);
       Eigen::MatrixXd g_model =  convert_armadillo_to_engine_matrix(UMAT);
       std::vector<double> vec = spectral_embedding(g_model, spectral_params);
       specmodelvec[i] = vec;
    }
    mat arrayLL;
    arrayLL.set_size(m_test_data.n_rows,nClass);
    for(int i=0; i < m_test_data.n_rows; ++i)
    {
      mat SEQ = m_test_data(i,0);
     // SEQ.print("SEQ");
      int SizeSEQ = SEQ.n_rows;
      for (int cl =0; cl < nClass; ++cl)
      {
          mat WSOM = m_som_codebook_list(cl,0);
          mat UMAT = m_umatrix_graph_list(cl,0);
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

          Eigen::MatrixXd g_data =  convert_armadillo_to_engine_matrix(GraphData);
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
    rowvec labelTrue = m_test_label.t();
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

    ///!!!!!!!!! РАСКОМЕНТИРОВАТЬ ДЛЯ АНСАМБЛЯ
   /* mat DP = mapminmax(arrayLL,0,1);
    DP.save("classif_2.result");*/

}

void class_npmpgm_model::test_ensemble()
{
    read_test_data("data/characterTestData.csv");
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

    rowvec labelTrue = m_test_label.t();


    rowvec labelDetect_MAX = calc_label_detect(arrayLL_MAX);
    rowvec labelDetect_MIN = calc_label_detect(arrayLL_MIN);
    rowvec labelDetect_SUM = calc_label_detect(arrayLL_SUM);
    rowvec labelDetect_AVR = calc_label_detect(arrayLL_AVR);
    rowvec labelDetect_PRO = calc_label_detect(arrayLL_PRO);
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



mat class_npmpgm_model::logsumexp(mat a, int dim)
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

double class_npmpgm_model::hmm_filter(mat initDist, mat transmat, mat softev)
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

void class_npmpgm_model::quality(rowvec labeldetect, rowvec labeltrue, int nClass, double& fmesure, double& precision, double& recall)
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
    if ((precision + recall) != 0)
    {
        fmesure = 2*precision*recall / (precision + recall);
    }
    else
    {
        fmesure = 0;
    }

    return ;
}

Eigen::MatrixXd class_npmpgm_model::convert_armadillo_to_engine_matrix(mat matrix)
{
    Eigen::MatrixXd G_eigen(matrix.n_rows,matrix.n_cols);
    for(int i = 0; i < matrix.n_rows; i++) {
      for(int j = 0; j < matrix.n_cols; j++) {
        G_eigen(i,j) = matrix(i,j);
      }
    }
    return G_eigen;
}

mat class_npmpgm_model::calculate_dist_node_matrix(mat codebook)
{
    mat Dist;
    int numNod = codebook.n_rows;
    Dist.set_size(numNod,numNod);
    Dist.zeros();
    for(int i=0; i < numNod; ++i)
    {
        int som_x1 = 0;
        int som_y1 = 0;
        two_from_one((ulong)i,(ushort)m_nsom_y,(ushort&)som_x1,(ushort&)som_y1);
        for(int j=0; j < numNod; ++j)
        {
            int som_x2 = 0;
            int som_y2 = 0;
            two_from_one((ulong)j,(ushort)m_nsom_y,(ushort&)som_x2,(ushort&)som_y2);
            float tmp = 0.0f;
            if (m_map_type == "planar")
            {
                tmp = euclideanDistanceOnPlanarMap(som_x1, som_y1, som_x2, som_y2);
            }
            if (m_map_type == "toroid")
            {
               tmp = euclideanDistanceOnToroidMap(som_x1, som_y1, som_x2, som_y2, m_nsom_x, m_nsom_y);
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

double class_npmpgm_model::get_l2norm(std::vector<double> vec)
{
    double res = 0;
    for (int i= 0; i < vec.size(); ++i)
    {
        res = res + pow(vec[i],2);
    }
    res = sqrt(res);
    return res;
}

void class_npmpgm_model::two_from_one (ulong z, ushort max_y, ushort& x, ushort& y)
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

mat class_npmpgm_model::mapminmax(mat matrix, double ymin, double ymax)
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

rowvec class_npmpgm_model::calc_label_detect(mat arrayLL)
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



