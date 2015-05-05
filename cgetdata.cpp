#include "cgetdata.h"
#include <QTime>
#include <QFile>
#include <QTextStream>
#include <QStringList>
//const int C_MAX_ROWS = 999999;
const int C_MAX_ROWS = 1100000;
const int C_MAX_COLUMNS = 20;
bool CGetData::get_mat_from_file(const QString& namefile, mat &data, mat &label)
{
    QFile file(namefile);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        return false;
    label.set_size(C_MAX_ROWS,2);
    data.set_size(C_MAX_ROWS,C_MAX_COLUMNS);
    QTextStream in(&file);
    int lineNum = 0;
    bool ok;
    int desi = 0;  
    while (!in.atEnd()) {
        QString line = in.readLine();
        QStringList strListVal = line.split(",");
        if (strListVal.size() > 3)
        {
           // label.resize(lineNum+1,2);
            label(lineNum,0) = strListVal.at(1).toInt(&ok);
            label(lineNum,1) = strListVal.at(2).toInt(&ok);
            desi = strListVal.size() - 3;
            for(int i=0;i< desi;++i)
            {
               // data.resize(lineNum+1,desi);
                data(lineNum,i) = strListVal.at(3+i).toDouble(&ok);

            }
        }       
        lineNum++;
    }
    label.resize(lineNum,2);
    data.resize(lineNum,desi);
    int g = 0;
    return true;

}

bool CGetData::get_cell_from_file(const QString& namefileData,field<mat>& data, mat& label)
{
    mat Data;
    mat Label;
    if (!CGetData::get_mat_from_file(namefileData,Data,Label))
    {
      return false;
    }


    field<mat> dataTemp;
    dataTemp.set_size(Data.n_rows,1);
    int predSeq = Label(0,0);
    int predLab = Label(0,1);
    mat A;
    int index = 0;
    int indexField = 0;
    label.set_size(C_MAX_ROWS,1);
    A.set_size(C_MAX_ROWS,Data.n_cols);
    for (int i = 0; i < Data.n_rows; ++i)
    {
       if (predSeq != Label(i,0))
       {
           A.resize(index,Data.n_cols);
           index = 0;
           dataTemp(indexField,0) = A;
           label(indexField,0) = predLab;
           int predLab0 = label(0,0);
           A.clear();
           A.set_size(C_MAX_ROWS,Data.n_cols);
           indexField++;
       }

       A(index,span::all) =  Data(i,span::all);
       index++;
       predSeq = Label(i,0);
       predLab = Label(i,1);
    }
    A.resize(index,Data.n_cols);
    dataTemp(indexField,0) = A;
    label(indexField,0) = predLab;
    data.set_size(indexField+1,1);
    label.resize(indexField+1,1);
    for (int i = 0; i < data.n_rows; ++i)
    {
       data(i,0) =  dataTemp(i,0);
    }
    int t = 0;

    return true;
}

void CGetData::forming_train_data_for_som(field<mat> data, mat label, field<mat> &dataForSOM)
{
    int predLab = 0;
    predLab = label(0,0);
    mat A;
    int index = 0;
    int indexField = 0;
    int first_row = 0;
    int last_row = 0;
    int sizeRow =label(label.n_rows-1,0) -  label(0,0) + 1;
    A.set_size(C_MAX_ROWS,data(0,0).n_cols);
    dataForSOM.set_size(sizeRow,1);
    int des = 0;
    for (int i = 0; i < data.n_rows; ++i)
    {
       if (predLab != label(i,0))
       {
           A.resize(first_row,data(0,0).n_cols);
           first_row = 0;
           last_row = 0;

           dataForSOM(indexField,0) = A;

           A.clear();
           A.set_size(C_MAX_ROWS,data(0,0).n_cols);
           indexField++;
       }
       mat D = data(i,0);     
       last_row = first_row + D.n_rows - 1;       
       A(span(first_row, last_row),span::all) = D;
       first_row = last_row + 1;
      /* for (int k = 0; k < D.n_rows; ++k)
       {
           A.resize(index+1,D.n_cols);
           for (int j = 0; j < D.n_cols; ++j)
           {
             A(index,j) =  D(k,j);
           }
           index++;
       }*/

       predLab = label(i,0);
    }
    A.resize(first_row,data(0,0).n_cols);
    dataForSOM(indexField,0) = A;

    return ;

}



















