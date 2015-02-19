#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QTimer>
#include "qsomthread.h"
#include "armadillo"
#include <stdlib.h>
#include "cgetdata.h"
namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    ///Обработчик событий
    virtual bool event(QEvent* ev);

private:
    Ui::MainWindow *ui;

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
    /// матрциа переходов между состояними для каждого класса
    field<mat> m_MatrixTransactA;    
    ///количество классов для которых завершено обучение
    int m_nClassComplete;
    /// таймер времени обучения
    QTimer* m_timertrain;

private:
    /// читаем данные из двух csv файлов и заполняем m_TrainData m_TrainLabel m_TrainDataForSOM m_nClass
    void readTrainData(const QString &datafile, const QString &lablefile);
    /// читаем данные из двух csv файлов и заполняем m_TestData m_TestLabel
    void readTestData(const QString &datafile, const QString &lablefile);
    /// формирвоание набора файлов для оубчения SOM
    void writeSOMtrainfiles(const QString& patternfilename);
    /// формирвоание матрицы переходов между состояними
    mat formMatrixTransaction(mat codebook,int label);

    mat logsumexp(mat a, int dim);

    double hmmFilter(mat initDist, mat transmat, mat softev);
public slots:
    /// анализ времени обучения
    void timeoutAnalysisTrainComplete();
    /// обучение
    void train();
    /// тестирование
    void test();
};

#endif // MAINWINDOW_H
