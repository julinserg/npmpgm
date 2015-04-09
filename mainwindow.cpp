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

NPMPGM_GUI::NPMPGM_GUI(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);



    m_ModelThread = new NPMPGMThread();
    m_ModelThread->moveToThread(m_ModelThread);
    m_ModelThread->start();
    bool isconnected1 = false;
    bool isconnected2 = false;
    bool isconnected3 = false;
    isconnected1 = QObject::connect(this,SIGNAL(readTrainData(QString)),m_ModelThread,SLOT(readTrainData(QString)),Qt::QueuedConnection);
    isconnected2 = QObject::connect(this,SIGNAL(setPathToModel(QString)),m_ModelThread,SLOT(setPathToModel(QString)),Qt::QueuedConnection);
    isconnected3 = QObject::connect(this,SIGNAL(train(int ,int , int , QString , int , int , QString , int , int , QString )),
            m_ModelThread,SLOT(train(int ,int , int , QString , int , int , QString , int , int , QString )),Qt::QueuedConnection);

    m_fileNameTrainData = C_DEF_PATH_TO_TRAIN_DATA;
    ui->leLoadTrainData->setText(C_DEF_PATH_TO_TRAIN_DATA);
    m_dirNameModel = C_DEF_PATH_TO_MODEL;
    ui->leSetModelPath->setText(C_DEF_PATH_TO_MODEL);

}

NPMPGM_GUI::~NPMPGM_GUI()
{
    delete ui;
}

void NPMPGM_GUI::on_pbSetTrainData_clicked()
{
   QString fileName = QFileDialog::getOpenFileName(this,
        tr("Load Train Data File"), "./", tr("Train Data File (*.csv)"));
   ui->leLoadTrainData->setText(fileName);
   m_fileNameTrainData = fileName;

}

void NPMPGM_GUI::on_pbLoadTrainData_clicked()
{
   emit readTrainData(m_fileNameTrainData);
}

void NPMPGM_GUI::on_pbSetModelPath_clicked()
{
    QString dirName = QFileDialog::getExistingDirectory(this,
         tr("Select Models Path"), "./");
    ui->leSetModelPath->setText(dirName);
    m_dirNameModel = dirName;
}

void NPMPGM_GUI::on_pbTrain_clicked()
{
    emit setPathToModel(m_dirNameModel);
    emit train(500,10, 10,"planar", 3, 1,"linear", 0.1, 0.01, "linear");
}








