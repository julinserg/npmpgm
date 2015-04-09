#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
//QString FileName = C_SOMFILENAMEDATA.arg(i).arg("csv");
//QString OutPrefix = C_SOMFILENAMECODEBOOK.arg(i);
//QSOMThread*  somthread = new QSOMThread();
//m_nSOM_X = 15;
//m_nSOM_Y = 15;
//m_mapType = "planar";
//somthread->setFileName(FileName.toStdString());
//somthread->setOutPrefix(OutPrefix.toStdString());
//somthread->setNumEpoch(1000);
//somthread->setSizeMap(m_nSOM_X,m_nSOM_Y,m_mapType);
//somthread->setRadiusParam(3,1,"linear");
//somthread->setScaleParam(0.1,0.01,"linear");
//somthread->setKernelType(0,4);
//somthread->setSaveParam(0,"");
//somthread->setObjectEvent(this,i);
//somthread->start();

const QString C_DEF_PATH_TO_TRAIN_DATA = "./data/characterTrainData.csv";
const QString C_DEF_PATH_TO_MODEL = "./model";

class_npmpgm_gui::class_npmpgm_gui(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);



    m_model = new class_npmpgm_model();
    m_model->moveToThread(m_model);
    m_model->start();
    bool isconnected1 = false;
    bool isconnected2 = false;
    bool isconnected3 = false;
    bool isconnected4 = false;
    isconnected1 = QObject::connect(this,SIGNAL(read_train_data(QString)),m_model,SLOT(read_train_data(QString)),Qt::QueuedConnection);
    isconnected2 = QObject::connect(this,SIGNAL(set_path_to_model(QString)),m_model,SLOT(set_path_to_model(QString)),Qt::QueuedConnection);
    isconnected3 = QObject::connect(this,SIGNAL(train(int ,int , int , QString , int , int , QString , int , int , QString )),
            m_model,SLOT(train(int ,int , int , QString , int , int , QString , int , int , QString )),Qt::QueuedConnection);
    isconnected4 = QObject::connect(m_model,SIGNAL(end_load_traindata(bool)),this,SLOT(end_load_traindata(bool)),Qt::QueuedConnection);

    m_filename_traindata = C_DEF_PATH_TO_TRAIN_DATA;
    ui->leLoadTrainData->setText(C_DEF_PATH_TO_TRAIN_DATA);
    m_dirname_model = C_DEF_PATH_TO_MODEL;
    ui->leSetModelPath->setText(C_DEF_PATH_TO_MODEL);
    ui->leStatusTrainData->setText("Данные не загружены.");

}

class_npmpgm_gui::~class_npmpgm_gui()
{
    delete ui;
}

void class_npmpgm_gui::on_pbSetTrainData_clicked()
{
   QString fileName = QFileDialog::getOpenFileName(this,
        tr("Load Train Data File"), "./", tr("Train Data File (*.csv)"));
   ui->leLoadTrainData->setText(fileName);
   m_filename_traindata = fileName;

}

void class_npmpgm_gui::on_pbLoadTrainData_clicked()
{
   ui->leStatusTrainData->setText("Идет загрузка данных...");
   emit read_train_data(m_filename_traindata);
}

void class_npmpgm_gui::on_pbSetModelPath_clicked()
{
    QString dirName = QFileDialog::getExistingDirectory(this,
         tr("Select Models Path"), "./");
    ui->leSetModelPath->setText(dirName);
    m_dirname_model = dirName;
}

void class_npmpgm_gui::on_pbTrain_clicked()
{
    emit set_path_to_model(m_dirname_model + "/");
    emit train(500,10, 10,"planar", 3, 1,"linear", 0.1, 0.01, "linear");
}

void class_npmpgm_gui::end_load_traindata(bool res)
{
   if (res)
   {
       ui->leStatusTrainData->setText("Данные успешно загружены.");
   }
   else
   {
       ui->leStatusTrainData->setText("Ошибка загрузки данных.");
   }


}








