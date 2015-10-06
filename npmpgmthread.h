#ifndef NPMPGMTHREAD_H
#define NPMPGMTHREAD_H

#include <QTimer>
#include "qsomthread.h"
#include "armadillo"
#include <stdlib.h>
#include "cgetdata.h"

class class_npmpgm_model : public QThread
{
    Q_OBJECT

public:
     class_npmpgm_model();
    ~class_npmpgm_model();
    ///Обработчик событий
    virtual bool event(QEvent* ev);
    void run();
    int getNumberClass()
    {
        return m_nclass;
    }
    int getNumberDataTest()
    {
        return m_test_data.n_rows;
    }

private:

    class_som_thread* m_somthread;

    /// обучающие данные (мапятся с метками по индексу строки)
    field<mat> m_train_data;
    /// метки классов (мапятся с данными по индексу строки)
    mat m_train_label;
    /// тестовые данные (мапятся с метками по индексу строки)
    field<mat> m_test_data;
    /// метки классов если таковые имеются (мапятся с данными по индексу строки)
    mat m_test_label;
    /// все оубчающие данные дял каждого класса представлены как одна матрциа
    ///(без разделений на последовательности)
    field<mat> m_train_data_for_som;
    /// количество классов
    int m_nclass;
    /// веса карты кохонена для каждого класса
    field<mat> m_som_codebook_list;
    /// графы карты кохонена для каждого класса
    field<mat> m_umatrix_graph_list;
    /// матрциа переходов между состояними для каждого класса
    field<mat> m_matrix_transact_a;
    /// вектор инициализации переходов между состояними для каждого класса
    field<mat> m_matrix_pi;
    ///количество классов для которых завершено обучение
    int m_nclassComplete;
    /// таймер времени обучения
    QTimer* m_timertrain;
    /// путь к файлам модели
    QString m_pathtomodel;

    string m_map_type;
    int m_nsom_x;
    int m_nsom_y;
private:   
    /// формирвоание набора файлов для оубчения SOM
    void write_som_trainfiles(const QString& patternfilename);
    /// формирвоание матрицы переходов между состояними
    void form_matrix_transaction(mat codebook,int label, mat& TR, mat& PI);

    mat logsumexp(mat a, int dim);

    double hmm_filter(mat initDist, mat transmat, mat softev);

    /*double FilterFwd(const Eigen::Map<MatrixType>& transmat, const Eigen::Map<MatrixType>& softev,
                     const Eigen::Map<VectorType>& init);*/

    void quality(rowvec labeldetect, rowvec labeltrue, int nClass, double& fmesure, double& precision, double& recall);

    /*Eigen::MatrixXd convert_armadillo_to_engine_matrix(mat matrix);*/

    mat calculate_dist_node_matrix(mat codebook);

    double get_l2norm(std::vector<double> vec);

    void two_from_one(ulong z, ushort max_y, ushort& x, ushort& y);

    mat mapminmax(mat matrix, double ymin, double ymax);

    rowvec calc_label_detect(mat arrayLL);
public slots:
    /// читаем данные из csv файла и заполняем m_train_data m_train_label m_train_data_for_som m_nclass
    void read_train_data(QString datafile);
    /// читаем данные из csv файла и заполняем m_test_data m_test_label
    void read_test_data(QString datafile);
    /// путь к файлам модели
    void set_path_to_model(QString path);
    /// загрузка обученной модели
    void read_test_model(QString path);
    /// обучение
    void train(int numEpoch,int nSOM_X, int nSOM_Y,
               QString mapType, int bRadius, int eRadiusm,
               QString typeRadius, int bScale, int eScale, QString typeScale);
    /// тестирование
    void test();
    /// тестирование сравнения графов
    //void test2();
    /// тестирование объединение классификаторов
    //void test_ensemble();
private slots:
    /// анализ времени обучения
    void timeout_analysis_train_complete();
signals:
    void end_load_traindata(bool);
    void end_load_testdata(bool);
    void end_load_testmodel(bool);
    void number_class_complet(int);
    void result_testing(float,float,float);
    void number_testdata_complet(int);

};

#endif // NPMPGMTHREAD_H
