#include "cgetdata.h"
#include "csv.h"
mat CGetData::getMatFromFile(const QString& namefile)
{
    mat B;
    QList<QStringList> listfile = CSV::parseFromFile(namefile);
    B.set_size(listfile.size(),listfile.first().size());
    int i = 0;
    int j = 0;
    foreach (QStringList rowList, listfile) {
        j = 0;
        foreach (QString valStr, rowList) {
             bool ok;
             B(i,j) = valStr.toDouble(&ok);
             j++;
        }
        i++;
    }
    return B;
}

void CGetData::getCellFromFile(const QString& namefileData, const QString& namefileLable, field<mat>& data, mat& label)
{
    mat Data = CGetData::getMatFromFile(namefileData);
    mat Label = CGetData::getMatFromFile(namefileLable);
    field<mat> dataTemp;
    dataTemp.set_size(Data.n_rows,1);
    int predSeq = Label(0,0);
    int predLab = Label(0,1);
    mat A;
    int index = 0;
    int indexField = 0;
    for (int i = 0; i < Data.n_rows; ++i)
    {
       if (predSeq != Label(i,0))
       {
           index = 0;           
           label.resize(indexField+1,1);
           dataTemp(indexField,0) = A;
           label(indexField,0) = predLab;
           int predLab0 = label(0,0);
           A.clear();
           indexField++;
       }
       A.resize(index+1,Data.n_cols);
       for (int j = 0; j < Data.n_cols; ++j)
       {         
         A(index,j) =  Data(i,j);         
       }
       index++;
       predSeq = Label(i,0);
       predLab = Label(i,1);
    }   
    label.resize(indexField+1,1);
    dataTemp(indexField,0) = A;
    label(indexField,0) = predLab;
    data.set_size(indexField+1,1);
    for (int i = 0; i < data.n_rows; ++i)
    {
       data(i,0) =  dataTemp(i,0);
    }

    return ;
}

void CGetData::formingTrainDataForSOM(field<mat> data, mat label, field<mat> &dataForSOM)
{
    int predLab = 0;
    predLab = label(0,0);
    mat A;
    int index = 0;
    int indexField = 0;
    int sizeRow =label(label.n_rows-1,0) -  label(0,0) + 1;
    dataForSOM.set_size(sizeRow,1);
    for (int i = 0; i < data.n_rows; ++i)
    {
       if (predLab != label(i,0))
       {
           index = 0;           

           dataForSOM(indexField,0) = A;

           A.clear();
           indexField++;
       }
       mat D = data(i,0);
       for (int k = 0; k < D.n_rows; ++k)
       {
           A.resize(index+1,D.n_cols);
           for (int j = 0; j < D.n_cols; ++j)
           {
             A(index,j) =  D(k,j);
           }
           index++;
       }

       predLab = label(i,0);
    }  
    dataForSOM(indexField,0) = A;

    return ;

}



















