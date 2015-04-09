#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "npmpgmthread.h"



namespace Ui {
class MainWindow;
}

class class_npmpgm_gui : public QMainWindow
{
    Q_OBJECT

public:
    explicit class_npmpgm_gui(QWidget *parent = 0);
    ~class_npmpgm_gui();


private:
    Ui::MainWindow *ui;
    class_npmpgm_model* m_model;
    QString m_filename_traindata;
    QString m_dirname_model;
private:


private slots:
    void on_pbSetTrainData_clicked();
    void on_pbLoadTrainData_clicked();
    void on_pbSetModelPath_clicked();
    void on_pbTrain_clicked();
    void end_load_traindata(bool res);
signals:
    /// читаем данные из csv файла и заполняем m_train_data m_train_label m_train_data_for_som m_nclass
    void read_train_data(QString datafile);
    /// читаем данные из csv файла и заполняем m_test_data m_test_label
    void read_test_data(QString datafile);
    /// путь к файлам модели
    void set_path_to_model(QString path);
    /// обучение
    void train(int numEpoch,int nSOM_X, int nSOM_Y,
               QString mapType, int bRadius, int eRadiusm,
               QString typeRadius, int bScale, int eScale, QString typeScale);
    /// тестирование
    void test();
};

#endif // MAINWINDOW_H
