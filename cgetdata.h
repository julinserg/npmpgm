#ifndef CGETDATA_H
#define CGETDATA_H

#include <QObject>
#include "armadillo"
using namespace arma;
class CGetData
{   
public:    
    static mat getMatFromFile(const QString& namefile);
    static field<mat> getCellFromFile(const QString& namefileData, const QString& namefileLable);
};

#endif // CGETDATA_H
