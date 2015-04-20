#ifndef REALCLSNPMPGMTHREAD_H
#define REALCLSNPMPGMTHREAD_H


#include <stdlib.h>
#include <QThread>
#include "armadillo"
using namespace arma;
class thread_npmpgm_classify : public QThread
{
    Q_OBJECT

public:
     thread_npmpgm_classify();
    ~thread_npmpgm_classify();
    void run();
private:
    /// количество классов
    int m_nclass;
    /// веса карты кохонена для каждого класса
    field<mat> m_som_codebook_list;
    /// матрица переходов между состояними для каждого класса
    field<mat> m_matrix_transact_a;
    /// вектор инициализации переходов между состояними для каждого класса
    field<mat> m_matrix_pi;
private:
    mat logsumexp(mat a, int dim);

    double hmm_filter(mat initDist, mat transmat, mat softev);

public slots:
    /// загрузка обученной модели
    void set_model(QString path);
    /// классифицировать массив данных (временную последовательность)
    void classifyData( std::vector< std::vector<double> > data, int num);

signals:
    void set_model_result(bool);
    void result_classify(int,int);

};

#endif // REALCLSNPMPGMTHREAD_H
