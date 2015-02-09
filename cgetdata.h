#ifndef CGETDATA_H
#define CGETDATA_H

#include <QObject>
#include "armadillo"
using namespace arma;
class CGetData
{   
public:    
    static mat getMatFromFile(const QString& namefile);
    static void getCellFromFile(const QString& namefileData, const QString& namefileLable, field<mat>& data, mat& label);
    static void formingTrainDataForSOM(field<mat> data, mat label, field<mat>& dataForSOM);
};

#endif // CGETDATA_H
