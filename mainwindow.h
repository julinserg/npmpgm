#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "npmpgmthread.h"



namespace Ui {
class MainWindow;
}

class NPMPGM_GUI : public QMainWindow
{
    Q_OBJECT

public:
    explicit NPMPGM_GUI(QWidget *parent = 0);
    ~NPMPGM_GUI();


private:
    Ui::MainWindow *ui;
    NPMPGMThread* m_ModelThread;
    QString m_fileNameTrainData;
    QString m_dirNameModel;
private:


private slots:
    void on_pbSetTrainData_clicked();
    void on_pbLoadTrainData_clicked();
    void on_pbSetModelPath_clicked();
    void on_pbTrain_clicked();

signals:
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
};

#endif // MAINWINDOW_H
