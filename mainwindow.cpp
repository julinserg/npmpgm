#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "csv_parser.hpp"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    m_somthread = new QSOMThread();
    m_somthread->start();
    m_somthread->moveToThread(m_somthread);

    csv_parser csv("characterTestLabel.csv");

    //For getting values according to row number and column number. Remember it starts from (1,1) and not (0,0)
    string value = csv.get_value(3,4);

    //For getting a particular row in the CSV file
    string line = csv.get_line(3);

    //For getting number of fields in a particular row.
    int total_fields = csv.fields(csv.get_line(3));

    cout<<"Value in (3,4) :"<<value<<endl;
    cout<<"Row 3: "<<line<<endl;
    cout<<"Total fields in row 3:"<<total_fields<<endl;


}

MainWindow::~MainWindow()
{
    delete ui;
}
