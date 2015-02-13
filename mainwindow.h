#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "qsomthread.h"
#include "armadillo"
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

    QSOMThread* m_somthread;

    /// обучающие данные (мапятся с метками по индексу строки)
    field<mat> m_TrainData;
    /// метки классов (мапятся с данными по индексу строки)
    mat m_TrainLabel;
    /// все оубчающие данные дял каждого класса представлены как одна матрциа
    ///(без разделений на последовательности)
    field<mat> m_TrainDataForSOM;
    /// количество классов
    int m_nClass;


private:
    /// читаем данные из двух csv файлов и заполняем m_TrainData m_TrainLabel m_TrainDataForSOM m_nClass
    void readTrainData(const QString &datafile, const QString &lablefile);
    /// формирвоание набора файлов для оубчения SOM
    void writeSOMtrainfiles(const QString& patternfilename);
};

#endif // MAINWINDOW_H
