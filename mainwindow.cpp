#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    m_somthread = new QSOMThread();
    m_somthread->start();
    m_somthread->moveToThread(m_somthread);

}

MainWindow::~MainWindow()
{
    delete ui;
}
