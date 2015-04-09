#ifndef CGETDATA_H
#define CGETDATA_H

#include <QObject>
#include "armadillo"
using namespace arma;
class CGetData
{   
public:    
    static bool get_mat_from_file(const QString& namefile, mat& data, mat& label);
    static bool get_cell_from_file(const QString& namefileData, field<mat>& data, mat& label);
    static void forming_train_data_for_som(field<mat> data, mat label, field<mat>& dataForSOM);
};

#endif // CGETDATA_H
