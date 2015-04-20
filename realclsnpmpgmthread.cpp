#include "realclsnpmpgmthread.h"
#include <QDebug>
using namespace arma;
const QString C_SOMFILENAMEDATA = "som_train_class_%1.%2";
const QString C_SOMFILENAMECODEBOOK = "som_train_class_%1";
thread_npmpgm_classify::thread_npmpgm_classify()
{

}

thread_npmpgm_classify::~thread_npmpgm_classify()
{

}

void thread_npmpgm_classify::run()
{   
    exec();
}

void thread_npmpgm_classify::set_model(QString path)
{
    bool g1 = m_som_codebook_list.load(path.toStdString()+"weight_som.model");
    bool g2 = m_matrix_transact_a.load(path.toStdString()+"matrix_a.model");
    bool g3 = m_matrix_pi.load(path.toStdString()+"matrix_pi.model");
    if (g1 && g2 && g3)
    {
        emit set_model_result(true);
    }
    else
    {
        qDebug() << "No model!";
        emit set_model_result(false);
    }
}


void thread_npmpgm_classify::classifyData(std::vector<std::vector<double> > data, int num)
{
    int nClass = m_som_codebook_list.n_rows;
    if (nClass <= 0)
    {
        return ;
    }
    mat arrayLL;
    arrayLL.set_size(1,nClass);

    int row = data.size();
    std::vector<double> v = data[0];
    int col = v.size();
    mat SEQ(row,col);
    for (int i=0;i< row; ++i)
    {
        for (int j=0; j< col; ++j)
        {
           SEQ(i,j) = data[i][j];
        }
    }
    int SizeSEQ = SEQ.n_rows;

    for (int cl =0; cl < nClass; ++cl)
    {
      mat WSOM = m_som_codebook_list(cl,0);
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


      mat A = m_matrix_transact_a(cl,0);
      mat PI = m_matrix_pi(cl,0);
      double logp = hmm_filter(PI,A,B);
      logp = logp + sum(vectorise(L));
      arrayLL(0,cl) =logp;
    }

    int  labelDetect;
    rowvec vec = arrayLL(0,span::all);
    uword  index;
    double max = vec.max(index);
    labelDetect = index+1;

    emit result_classify(labelDetect,num);
}

mat thread_npmpgm_classify::logsumexp(mat a, int dim)
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

double thread_npmpgm_classify::hmm_filter(mat initDist, mat transmat, mat softev)
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
