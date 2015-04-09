#ifndef NPMPGMTHREAD_H
#define NPMPGMTHREAD_H

#include <QTimer>
#include "qsomthread.h"
#include "armadillo"
#include <stdlib.h>
#include "cgetdata.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>


class NPMPGMThread : public QThread
{
    Q_OBJECT

public:
     NPMPGMThread();
    ~NPMPGMThread();
    ///Обработчик событий
    virtual bool event(QEvent* ev);
    void run();
private:

    QSOMThread* m_somthread;

    /// обучающие данные (мапятся с метками по индексу строки)
    field<mat> m_TrainData;
    /// метки классов (мапятся с данными по индексу строки)
    mat m_TrainLabel;
    /// тестовые данные (мапятся с метками по индексу строки)
    field<mat> m_TestData;
    /// метки классов если таковые имеются (мапятся с данными по индексу строки)
    mat m_TestLabel;
    /// все оубчающие данные дял каждого класса представлены как одна матрциа
    ///(без разделений на последовательности)
    field<mat> m_TrainDataForSOM;
    /// количество классов
    int m_nClass;
    /// веса карты кохонена для каждого класса
    field<mat> m_SOMCodeBookList;
    /// графы карты кохонена для каждого класса
    field<mat> m_UmatrixGraphList;
    /// матрциа переходов между состояними для каждого класса
    field<mat> m_MatrixTransactA;
    /// вектор инициализации переходов между состояними для каждого класса
    field<mat> m_MatrixPI;
    ///количество классов для которых завершено обучение
    int m_nClassComplete;
    /// таймер времени обучения
    QTimer* m_timertrain;
    /// путь к файлам модели
    QString m_pathtomodel;

    string m_mapType;
    int m_nSOM_X;
    int m_nSOM_Y;
private:   
    /// формирвоание набора файлов для оубчения SOM
    void writeSOMtrainfiles(const QString& patternfilename);
    /// формирвоание матрицы переходов между состояними
    void formMatrixTransaction(mat codebook,int label, mat& TR, mat& PI);

    mat logsumexp(mat a, int dim);

    double hmmFilter(mat initDist, mat transmat, mat softev);

    void quality(rowvec labeldetect, rowvec labeltrue, int nClass, double& fmesure, double& precision, double& recall);

    Eigen::MatrixXd convertArmadilloToEngineMatrix(mat matrix);

    mat calculateDistNodeMatrix(mat codebook);

    double get_l2norm(std::vector<double> vec);

    void twoFromOne(ulong z, ushort max_y, ushort& x, ushort& y);

    mat mapminmax(mat matrix, double ymin, double ymax);

    rowvec calcLabelDetect(mat arrayLL);
public slots:
    /// читаем данные из csv файла и заполняем m_TrainData m_TrainLabel m_TrainDataForSOM m_nClass
    void readTrainData(QString datafile);
    /// читаем данные из csv файла и заполняем m_TestData m_TestLabel
    void readTestData(QString datafile);
    /// путь к файлам модели
    void setPathToModel(QString path);
    /// обучение
    void train(int numEpoch,int nSOM_X, int nSOM_Y,
               QString mapType, int bRadius, int eRadiusm,
               QString typeRadius, int bScale, int eScale, QString typeScale);
    /// тестирование
    void test();
    /// тестирование сравнения графов
    void test2();
    /// тестирование объединение классификаторов
    void testEnsemble();
private slots:
    /// анализ времени обучения
    void timeoutAnalysisTrainComplete();

};

#endif // NPMPGMTHREAD_H
