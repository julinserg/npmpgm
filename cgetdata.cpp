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

field<mat> CGetData::getCellFromFile(const QString& namefileData, const QString& namefileLable)
{
    mat Data = CGetData::getMatFromFile(namefileData);
    mat Label = CGetData::getMatFromFile(namefileLable);

    int predSeq = Label(1,0);
    mat A;
    int index = 0;
    field<mat> B;
    int indexField = 0;
    for (int i = 0; i < Data.n_rows; ++i)
    {
       if (predSeq != Label(i,0))
       {
           index = 0;
           B.set_size(indexField+1,1);
           B(indexField,0) = A;
           A.clear();
           indexField++;
       }
       A.set_size(index+1,Data.n_cols);
       for (int j = 0; j < Data.n_cols; ++j)
       {         
         A(index,j) =  Data(i,j);         
       }
       index++;
       predSeq = Label(i,0);
    }
    B.set_size(indexField+1,1);
    B(indexField,0) = A;
    return B;
}

