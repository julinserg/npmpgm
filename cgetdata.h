#ifndef CGETDATA_H
#define CGETDATA_H

#include <QObject>
#include "armadillo"
using namespace arma;
class CGetData
{   
public:    
    static bool getMatFromFile(const QString& namefile, mat& data, mat& label);
    static bool getCellFromFile(const QString& namefileData, field<mat>& data, mat& label);
    static void formingTrainDataForSOM(field<mat> data, mat label, field<mat>& dataForSOM);
};

#endif // CGETDATA_H
