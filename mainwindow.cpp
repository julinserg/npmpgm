#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <QTextStream>
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
const QString C_DEF_PATH_TO_TEST_DATA = "./data/characterTestData.csv";
const QString C_DEF_PATH_TO_MODEL = "./model";

class_npmpgm_gui::class_npmpgm_gui(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);



    m_model = new class_npmpgm_model();
    m_model->moveToThread(m_model);
    m_model->start();

    m_classify = new thread_npmpgm_classify();
    m_classify->moveToThread(m_classify);
    m_classify->start();

    bool isconnected = false;
    isconnected = QObject::connect(this,SIGNAL(read_train_data(QString)),m_model,SLOT(read_train_data(QString)),Qt::QueuedConnection);
    isconnected = QObject::connect(this,SIGNAL(set_path_to_model(QString)),m_model,SLOT(set_path_to_model(QString)),Qt::QueuedConnection);
    isconnected = QObject::connect(this,SIGNAL(train(int ,int , int , QString , int , int , QString , int , int , QString )),
            m_model,SLOT(train(int ,int , int , QString , int , int , QString , int , int , QString )),Qt::QueuedConnection);
    isconnected = QObject::connect(m_model,SIGNAL(end_load_traindata(bool)),this,SLOT(end_load_traindata(bool)),Qt::QueuedConnection);
    isconnected = QObject::connect(m_model,SIGNAL(number_class_complet(int)),this,SLOT(number_class_complet(int)),Qt::QueuedConnection);
    isconnected = QObject::connect(m_model,SIGNAL(end_load_testdata(bool)),this,SLOT(end_load_testdata(bool)),Qt::QueuedConnection);
    isconnected = QObject::connect(this,SIGNAL(read_test_data(QString)),m_model,SLOT(read_test_data(QString)),Qt::QueuedConnection);
    isconnected = QObject::connect(this,SIGNAL(read_test_model(QString)),m_model,SLOT(read_test_model(QString)),Qt::QueuedConnection);
    isconnected = QObject::connect(this,SIGNAL(test()),m_model,SLOT(test()),Qt::QueuedConnection);
    isconnected = QObject::connect(m_model,SIGNAL(result_testing(float,float,float)),this,SLOT(result_testing(float,float,float)),Qt::QueuedConnection);
    isconnected = QObject::connect(m_model,SIGNAL(number_testdata_complet(int)),this,SLOT(number_testdata_complet(int)),Qt::QueuedConnection);
    isconnected = QObject::connect(m_model,SIGNAL(end_load_testmodel(bool)),this,SLOT(end_load_testmodel(bool)),Qt::QueuedConnection);


    qRegisterMetaType<std::vector< std::vector<double> > >("std::vector< std::vector<double> >");
    isconnected = QObject::connect(this,SIGNAL(classifyData(std::vector< std::vector<double> >, int)),
                                   m_classify,SLOT(classifyData(std::vector< std::vector<double> >, int)),Qt::QueuedConnection);
    isconnected = QObject::connect(this,SIGNAL(set_model(QString)),m_classify,SLOT(set_model(QString)),Qt::QueuedConnection);
    isconnected = QObject::connect(m_classify,SIGNAL(set_model_result(bool)),this,SLOT(set_model_result(bool)),Qt::QueuedConnection);
    isconnected = QObject::connect(m_classify,SIGNAL(result_classify(int,int)),this,SLOT(result_classify(int,int)),Qt::QueuedConnection);
    m_filename_traindata = C_DEF_PATH_TO_TRAIN_DATA;
    ui->leLoadTrainData->setText(C_DEF_PATH_TO_TRAIN_DATA);
    m_filename_testdata = C_DEF_PATH_TO_TEST_DATA;
    ui->leLoadTestData->setText(C_DEF_PATH_TO_TEST_DATA);
    m_dirname_model = C_DEF_PATH_TO_MODEL;
    ui->leSetModelPath->setText(C_DEF_PATH_TO_MODEL);
    m_dirname_testmodel = C_DEF_PATH_TO_MODEL;
    ui->leLoadTestModel->setText(C_DEF_PATH_TO_MODEL);
    m_dirname_classifymodel = C_DEF_PATH_TO_MODEL;
    ui->leLoadClassifyModel->setText(C_DEF_PATH_TO_MODEL);
    ui->leStatusTrainData->setText("Данные не загружены.");
    ui->leStatusTrainData->setStyleSheet("background-color: #FFD700");
    ui->leStatusTestData->setText("Данные не загружены.");
    ui->leStatusTestData->setStyleSheet("background-color: #FFD700");
    ui->leStatusTestModel->setText("Модель не загружена.");
    ui->leStatusTestModel->setStyleSheet("background-color: #FFD700");
    ui->leStatusClassifyModel->setText("Модель не загружена.");
    ui->leStatusClassifyModel->setStyleSheet("background-color: #FFD700");
    ui->le_Result_Class->setText("---");

}

class_npmpgm_gui::~class_npmpgm_gui()
{
    delete ui;
}

void class_npmpgm_gui::on_pbSetTrainData_clicked()
{
   ui->leStatusTrainData->setText("Данные не загружены.");
   ui->leStatusTrainData->setStyleSheet("background-color: #FFD700");
   QString fileName = QFileDialog::getOpenFileName(this,
        tr("Выбор файла с обучающими данными"), "./", tr("Обучающие данные (*.csv)"));
   ui->leLoadTrainData->setText(fileName);
   m_filename_traindata = fileName;

}

void class_npmpgm_gui::on_pbLoadTrainData_clicked()
{
   ui->leStatusTrainData->setText("Идет загрузка данных...");
   ui->leStatusTrainData->setStyleSheet("background-color: #FFFFFF");
   emit read_train_data(m_filename_traindata);
}

void class_npmpgm_gui::on_pbSetTestData_clicked()
{
    ui->leStatusTestData->setText("Данные не загружены.");
    ui->leStatusTestData->setStyleSheet("background-color: #FFD700");
    QString fileName = QFileDialog::getOpenFileName(this,
         tr("Выбор файла с тестовыми данными"), "./", tr("Тестовые данные (*.csv)"));
    ui->leLoadTestData->setText(fileName);
    m_filename_testdata = fileName;
}

void class_npmpgm_gui::on_pbLoadTestData_clicked()
{
    ui->leStatusTestData->setText("Идет загрузка данных...");
    ui->leStatusTestData->setStyleSheet("background-color: #FFFFFF");
    emit read_test_data(m_filename_testdata);
}

void class_npmpgm_gui::on_pbSetTestModel_clicked()
{
    ui->leStatusTestModel->setText("Данные не загружены.");
    ui->leStatusTestModel->setStyleSheet("background-color: #FFD700");
    QString dirName = QFileDialog::getExistingDirectory(this,
         tr("Выбрать папку с файлами модели"), "./");
    ui->leLoadTestModel->setText(dirName);
    m_dirname_testmodel = dirName;
}

void class_npmpgm_gui::on_pbLoadTestModel_clicked()
{
    ui->leStatusTestModel->setText("Идет загрузка модели...");
    ui->leStatusTestModel->setStyleSheet("background-color: #FFFFFF");
    emit read_test_model(m_dirname_testmodel + "/");
}

void class_npmpgm_gui::on_pbSetModelPath_clicked()
{
    QString dirName = QFileDialog::getExistingDirectory(this,
         tr("Указать папку сохранения файлов модели"), "./");
    ui->leSetModelPath->setText(dirName);
    m_dirname_model = dirName;
}

void class_npmpgm_gui::on_pbSetClassifyData_clicked()
{
    QString namefile = QFileDialog::getOpenFileName(this,
         tr("Выбрать файл с данными для классификации"), "./", tr("Файл с данными (*.csv)"));
    QFile file(namefile);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
    {
         ui->leLoadClassifyData->setStyleSheet("background-color: #B22222");
         return ;
    }

    int lineNum = 0;
    QTextStream in(&file);
    QString line = in.readLine();
    QStringList strListVal = line.split(",");
    bool ok;
    int row = strListVal.at(0).toInt(&ok);
    int col = strListVal.at(1).toInt(&ok);
    if (row <= 1 || col <= 1)
    {
        ui->leLoadClassifyData->setStyleSheet("background-color: #B22222");
        return ;
    }
    std::vector< std::vector<double> > matrixdata(row, std::vector<double>(col));
    while (!in.atEnd()) {
        QString line = in.readLine();
        QStringList strListVal = line.split(",");

        for(int i=0;i< strListVal.size();++i)
        {
            bool ok;
            matrixdata[lineNum][i] = strListVal.at(i).toDouble(&ok);
        }
        lineNum++;
    }
    int g = 0;
    ui->leLoadClassifyData->setText(namefile);
    ui->leLoadClassifyData->setStyleSheet("background-color: #32CD32");
    m_matrixdata = matrixdata;
}

void class_npmpgm_gui::on_pbSetClassifyModel_clicked()
{
    ui->leStatusClassifyModel->setText("Данные не загружены.");
    ui->leStatusClassifyModel->setStyleSheet("background-color: #FFD700");
    QString dirName = QFileDialog::getExistingDirectory(this,
         tr("Выбрать папку с файлами модели"), "./");
    ui->leLoadClassifyModel->setText(dirName);
    m_dirname_classifymodel = dirName;
}

void class_npmpgm_gui::on_pbLoadClassifyModel_clicked()
{
    ui->leStatusClassifyModel->setText("Идет загрузка модели...");
    ui->leStatusClassifyModel->setStyleSheet("background-color: #FFFFFF");
    emit set_model(m_dirname_classifymodel + "/");
}

void class_npmpgm_gui::on_pbTrain_clicked()
{
    emit set_path_to_model(m_dirname_model + "/");
    QString typeScale = "linear";
    QString typeRadius = "linear";
    QString typeMap = "planar";
    switch (ui->cb_mapType->currentIndex()) {
    case 0:
        typeMap = "planar";
        break;
    case 1:
        typeMap = "toroid";
        break;
    }
    switch (ui->cb_typeRadius->currentIndex()) {
    case 0:
        typeRadius = "linear";
        break;
    case 1:
        typeRadius = "exponential";
        break;
    }
    switch (ui->cb_typeScale->currentIndex()) {
    case 0:
        typeScale = "linear";
        break;
    case 1:
        typeScale = "exponential";
        break;
    }
    m_nClass = m_model->getNumberClass();
    ui->prbar_Learn->setMaximum(0);
    ui->prbar_Learn->setMinimum(0);
    ui->prbar_Learn->reset();
    ui->leStatusTestModel->setText("Модель не загружена.");
    ui->leStatusTestModel->setStyleSheet("background-color: #FFD700");
    emit train(ui->spb_numEpoch->value(),ui->spb_nSOM_X->value(), ui->spb_nSOM_Y->value(),typeMap
               ,ui->spb_bRadius->value(), ui->spb_eRadius->value(),typeRadius
               ,ui->spb_bScale->value(),ui->spb_eScale->value(), typeScale);
}

void class_npmpgm_gui::on_pbTest_clicked()
{
    m_nTestData = m_model->getNumberDataTest();
    ui->prbar_Test->reset();
    emit test();
}

void class_npmpgm_gui::on_pbClassify_clicked()
{
    ui->le_Result_Class->setText("---");
    emit classifyData(m_matrixdata,0);
}

void class_npmpgm_gui::end_load_traindata(bool res)
{
   if (res)
   {
       ui->leStatusTrainData->setText("Данные успешно загружены.");
       ui->leStatusTrainData->setStyleSheet("background-color: #32CD32");
   }
   else
   {
       ui->leStatusTrainData->setText("Ошибка загрузки данных.");
       ui->leStatusTrainData->setStyleSheet("background-color: #B22222");
   }


}

void class_npmpgm_gui::end_load_testdata(bool res)
{
    if (res)
    {
        ui->leStatusTestData->setText("Данные успешно загружены.");
        ui->leStatusTestData->setStyleSheet("background-color: #32CD32");
    }
    else
    {
        ui->leStatusTestData->setText("Ошибка загрузки данных.");
        ui->leStatusTestData->setStyleSheet("background-color: #B22222");
    }
}

void class_npmpgm_gui::end_load_testmodel(bool res)
{
    if (res)
    {
        ui->leStatusTestModel->setText("Модель успешно загружена.");
        ui->leStatusTestModel->setStyleSheet("background-color: #32CD32");
    }
    else
    {
        ui->leStatusTestModel->setText("Ошибка загрузки модели.");
        ui->leStatusTestModel->setStyleSheet("background-color: #B22222");
    }

}

void class_npmpgm_gui::number_class_complet(int n)
{
    ui->prbar_Learn->setMaximum(m_nClass);
    ui->prbar_Learn->setMinimum(0);
    ui->prbar_Learn->setValue(n);
}

void class_npmpgm_gui::number_testdata_complet(int n)
{
    ui->prbar_Test->setMaximum(m_nTestData);
    ui->prbar_Test->setMinimum(0);
    ui->prbar_Test->setValue(n+1);
}

void class_npmpgm_gui::result_testing(float fmeasure, float precision, float recall)
{
   ui->le_Result_Precision->setText(QString::number(precision,'f',5));
   ui->le_Result_Recall->setText(QString::number(recall,'f',5));
   ui->le_Result_Fmeasure->setText(QString::number(fmeasure,'f',5));
}

void class_npmpgm_gui::set_model_result(bool res)
{
    if (res)
    {
        ui->leStatusClassifyModel->setText("Модель успешно загружена.");
        ui->leStatusClassifyModel->setStyleSheet("background-color: #32CD32");
    }
    else
    {
        ui->leStatusClassifyModel->setText("Ошибка загрузки модели.");
        ui->leStatusClassifyModel->setStyleSheet("background-color: #B22222");
    }
}

void class_npmpgm_gui::result_classify(int nclass, int num)
{
    ui->le_Result_Class->setText(QString::number(nclass));
}








