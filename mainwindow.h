#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "npmpgmthread.h"



namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();


private:
    Ui::MainWindow *ui;
    NPMPGMThread* m_ModelThread;

private:


private slots:
    void on_pbLoadTrainData_clicked();
    void on_pbSetModelPath_clicked();
    void on_pbTrain_clicked();

signals:
    /// читаем данные из csv файла и заполняем m_TrainData m_TrainLabel m_TrainDataForSOM m_nClass
    void readTrainData(const QString &datafile);
    /// читаем данные из csv файла и заполняем m_TestData m_TestLabel
    void readTestData(const QString &datafile);
    /// путь к файлам модели
    void setPathToModel(const QString &path);
    /// обучение
    void train(int numEpoch,int nSOM_X, int nSOM_Y,
               string mapType, int bRadius, int eRadiusm,
               string typeRadius, int bScale, int eScale, string typeScale);
    /// тестирование
    void test();
};

#endif // MAINWINDOW_H
