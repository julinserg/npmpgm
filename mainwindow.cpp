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

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);  

    m_ModelThread = new NPMPGMThread();
    m_ModelThread->moveToThread(m_ModelThread);
    bool isconnected1 = false;
    bool isconnected2 = false;
    bool isconnected3 = false;
    isconnected1 = connect(this,SIGNAL(readTrainData(const QString &)),m_ModelThread,SLOT(readTrainData(const QString &)),Qt::QueuedConnection);
    isconnected2 = connect(this,SIGNAL(setPathToModel(const QString &)),m_ModelThread,SLOT(setPathToModel(const QString &)),Qt::QueuedConnection);
    isconnected3 = connect(this,SIGNAL(train(int ,int , int , string , int , int , string , int , int , string )),
            m_ModelThread,SLOT(train(int ,int , int , string , int , int , string , int , int , string )),Qt::QueuedConnection);

    m_ModelThread->start();



}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pbLoadTrainData_clicked()
{
   QString fileName = QFileDialog::getOpenFileName(this,
        tr("Load Train Data File"), "./", tr("Train Data File (*.csv)"));
   emit readTrainData(fileName);
}

void MainWindow::on_pbSetModelPath_clicked()
{
  /*  QString dirName = QFileDialog::getExistingDirectory(this,
         tr("Select Models Path"), "./");*/
    emit setPathToModel(/*dirName*/"123123");
}

void MainWindow::on_pbTrain_clicked()
{
    emit train(500,10, 10,"planar", 3, 1,"linear", 0.1, 0.01, "linear");
}








