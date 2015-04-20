#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "npmpgmthread.h"
#include "realclsnpmpgmthread.h"


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
    thread_npmpgm_classify* m_classify;
    QString m_filename_traindata;
    QString m_filename_testdata;
    QString m_dirname_model;
    QString m_dirname_testmodel;
    QString m_dirname_classifymodel;
    int m_nClass;
    int m_nTestData;
    std::vector< std::vector<double> > m_matrixdata;
private:


private slots:
    void on_pbSetTrainData_clicked();
    void on_pbLoadTrainData_clicked();
    void on_pbSetTestData_clicked();
    void on_pbLoadTestData_clicked();
    void on_pbSetTestModel_clicked();
    void on_pbLoadTestModel_clicked();
    void on_pbSetModelPath_clicked();
    void on_pbSetClassifyData_clicked();
    void on_pbSetClassifyModel_clicked();
    void on_pbLoadClassifyModel_clicked();

    void on_pbTrain_clicked();
    void on_pbTest_clicked();
    void on_pbClassify_clicked();
    void end_load_traindata(bool res);
    void end_load_testdata(bool res);
    void end_load_testmodel(bool res);
    void number_class_complet(int n);
    void number_testdata_complet(int n);
    void result_testing(float fmeasure,float precision,float recall);
    void set_model_result(bool res);
    void result_classify(int nclass,int num);
signals:
    /// читаем данные из csv файла и заполняем m_train_data m_train_label m_train_data_for_som m_nclass
    void read_train_data(QString datafile);
    /// читаем данные из csv файла и заполняем m_test_data m_test_label
    void read_test_data(QString datafile);
    /// путь к файлам модели
    void set_path_to_model(QString path);
    /// путь к файлам модели
    void set_model(QString path);
    /// обучение
    void train(int numEpoch,int nSOM_X, int nSOM_Y,
               QString mapType, int bRadius, int eRadiusm,
               QString typeRadius, int bScale, int eScale, QString typeScale);
    /// тестирование
    void test();
    /// загрузка обученной модели
    void read_test_model(QString path);
    /// классифицировать массив данных (временную последовательность)
    void classifyData( std::vector< std::vector<double> > data, int num);

};

#endif // MAINWINDOW_H
