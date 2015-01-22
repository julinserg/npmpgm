#ifndef QSOMTHREAD_H
#define QSOMTHREAD_H

#include "QThread"

#include <string>
#include <algorithm>

using namespace std;

#if HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_MPI
#include <mpi.h>
#endif

/// The neighbor_fuct value below which we consider
/// the impact zero for a given node in the map
#define NEIGHBOR_THRESHOLD 0.05

#include <cmath>
#include <cstdlib>
#include <iostream>

#ifdef _WIN32
#include "Windows/unistd.h"
#include "Windows/getopt.h"
#else
#include <unistd.h>
#include <getopt.h>
#endif

/// For synchronized timing
#ifndef MPI_WTIME_IS_GLOBAL
#define MPI_WTIME_IS_GLOBAL 1
#endif



/// Sparse structures and routines
struct svm_node
{
    int index;
    float value;
};

/// Core data structures
struct core_data
{
    float *codebook;
    int *globalBmus;
    float *uMatrix;
    int codebook_size;
    int globalBmus_size;
    int uMatrix_size;
};

class QSOMThread : public QThread
{
    Q_OBJECT
public:

    QSOMThread();

    ~QSOMThread();

    void run();

    // путь к файлу с данными
    void setFileName(string inFilename)
    {
        m_inFilename  = inFilename;
    }
    // какой то аут префикс
    void setOutPrefix(string outPrefix)
    {
        m_outPrefix = outPrefix;
    }
    // количество эпох обучения
    void setNumEpoch(uint number)
    {
        m_nEpoch = number;
    }
    // размер карты и тип карты(planar or toroid)
    void setNumEpoch(uint nSomX, uint nSomY, string mapType)
    {
        m_nSomX = nSomX;
        m_nSomY = nSomY;
        m_mapType = mapType;
    }
    // начальный радиус, конечный радиус, тип изменения радиуса (linear or exponential)
    void setRadiusParam(uint radius0, uint radiusN, string radiusCooling )
    {
        m_radius0 = radius0;
        m_radiusN = radiusN;
        m_radiusCooling = radiusCooling;
    }
    // начальный коэффициент обучения, конечный коэффициент обучения,
    // тип изменения коэффициента (linear or exponential)
    void setScaleParam(uint scale0, uint scaleN, string scaleCooling )
    {
        m_scale0 = scale0;
        m_scaleN = scaleN;
        m_scaleCooling = scaleCooling;
    }
    // тип вычислительного устрйоства и колчиество процессов/потоков
    /*#define DENSE_CPU 0
    #define DENSE_GPU 1
    #define SPARSE_CPU 2*/
    void setKernelType(uint kernelType, int numberProc)
    {
        m_kernelType = kernelType;
        m_nProcs = numberProc;
    }
    // тип сохраняемых в файл данных модели и имя файла
    void setSaveParam(uint snapshots, string initialCodebookFilename)
    {
       m_snapshots =  snapshots;
       m_initialCodebookFilename = initialCodebookFilename;
    }

 private:

    uint m_nEpoch;
    uint m_nSomX;
    uint m_nSomY;
    uint m_kernelType;
    string m_mapType;
    uint m_radius0;
    uint m_radiusN;
    string m_radiusCooling;
    float m_scale0;
    float m_scaleN;
    string m_scaleCooling;
    uint m_snapshots;
    string m_inFilename;
    string m_initialCodebookFilename;
    string m_outPrefix;
    int m_rank;
    int m_nProcs;
    int m_nThread;

 public slots:

    void startTrain();

};

#endif // QSOMTHREAD_H
