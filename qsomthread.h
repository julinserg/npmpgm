#ifndef QSOMTHREAD_H
#define QSOMTHREAD_H

#include "QThread"
#include <QEvent>
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

#ifdef _MSC_VER
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

float euclideanDistanceOnToroidMap(const unsigned int som_x, const unsigned int som_y, const unsigned int x, const unsigned int y, const unsigned int nSomX, const unsigned int nSomY);
float euclideanDistanceOnPlanarMap(const unsigned int som_x, const unsigned int som_y, const unsigned int x, const unsigned int y);
float getWeight(float distance, float radius, float scaling);
int saveCodebook(string cbFileName, float *codebook,
                unsigned int nSomX, unsigned int nSomY, unsigned int nDimensions);
float *calculateUMatrix(float *codebook, unsigned int nSomX,
             unsigned int nSomY, unsigned int nDimensions, string mapType);
int saveUMatrix(string fname, float *uMatrix, unsigned int nSomX,
              unsigned int nSomY);
int saveBmus(string filename, int *bmus, unsigned int nSomX,
             unsigned int nSomY, unsigned int nVectors);
//void printMatrix(float *A, int nRows, int nCols);
float *readMatrix(const string inFilename,
                  unsigned int &nRows, unsigned int &nCols);
double *readMatrixDoubleType(const string inFilename,
                  unsigned int &nRows, unsigned int &nCols);
void readSparseMatrixDimensions(const string filename, unsigned int &nRows,
                            unsigned int &nColumns);
svm_node** readSparseMatrixChunk(const string filename, unsigned int nRows,
                                 unsigned int nRowsToRead,
                                 unsigned int rowOffset);
core_data trainOneEpoch(int itask, float *data, svm_node **sparseData,
           core_data coreData, unsigned int nEpoch, unsigned int currentEpoch,
           bool enableCalculatingUMatrix,
           unsigned int nSomX, unsigned int nSomY,
           unsigned int nDimensions, unsigned int nVectors,
           unsigned int nVectorsPerRank,
           unsigned int radius0, unsigned int radiusN,
           string radiusCooling,
           float scale0, float scaleN,
           string scaleCooling,
           unsigned int kernelType, string mapType);
void train(int itask, float *data, svm_node **sparseData,
           unsigned int nSomX, unsigned int nSomY,
           unsigned int nDimensions, unsigned int nVectors,
           unsigned int nVectorsPerRank, unsigned int nEpoch,
           unsigned int radius0, unsigned int radiusN,
           string radiusCooling,
           float scale0, float scaleN,
           string scaleCooling,
           string outPrefix, unsigned int snapshots,
           unsigned int kernelType, string mapType,
           string initialCodebookFilename);
void trainOneEpochDenseCPU(int itask, float *data, float *numerator,
                           float *denominator, float *codebook,
                           unsigned int nSomX, unsigned int nSomY,
                           unsigned int nDimensions, unsigned int nVectors,
                           unsigned int nVectorsPerRank, float radius,
                           float scale, string mapType, int *globalBmus);
void trainOneEpochSparseCPU(int itask, svm_node **sparseData, float *numerator,
                           float *denominator, float *codebook,
                           unsigned int nSomX, unsigned int nSomY,
                           unsigned int nDimensions, unsigned int nVectors,
                           unsigned int nVectorsPerRank, float radius,
                           float scale, string mapType, int *globalBmus);
void initializeCodebook(unsigned int seed, float *codebook, unsigned int nSomX,
                        unsigned int nSomY, unsigned int nDimensions);

unsigned int *readLrnHeader(string inFilename, unsigned int &nRows, unsigned int &nColumns);
unsigned int *readWtsHeader(string inFilename, unsigned int &nRows, unsigned int &nColumns);
void getMatrixDimensions(string inFilename, unsigned int &nRows, unsigned int &nColumns);

extern "C" {
#ifdef CUDA
void setDevice(int commRank, int commSize);
void freeGpu();
void initializeGpu(float *hostData, int nVectorsPerRank, int nDimensions, int nSomX, int nSomY);
void trainOneEpochDenseGPU(int itask, float *data, float *numerator,
                           float *denominator, float *codebook,
                           unsigned int nSomX, unsigned int nSomY,
                           unsigned int nDimensions, unsigned int nVectors,
                           unsigned int nVectorsPerRank, float radius,
                           float scale, string mapType, int *globalBmus);
#endif
void my_abort(int err);
}

#define SOM_THREAD_STOP				9999

class ESOMThreadStop: public QEvent
{
public:
    ESOMThreadStop():QEvent((QEvent::Type)SOM_THREAD_STOP)
    {
    }

    ~ESOMThreadStop()
    {
    }
    /// метка класса для которого завершено обучение
    int labelClass;
};

class class_som_thread : public QThread
{
    Q_OBJECT
public:

    class_som_thread();

    ~class_som_thread();

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
    void setSizeMap(uint nSomX, uint nSomY, string mapType)
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
    //отправлять обытия этому объекту
    void setObjectEvent(QObject* object, int labelClass)
    {
        m_objectEvent = object;
        m_labelClass = labelClass;
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

    QObject* m_objectEvent;
    int m_labelClass;

};

#endif // QSOMTHREAD_H
