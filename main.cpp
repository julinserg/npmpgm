#include "mainwindow.h"
#include <QApplication>

#include <string>
#include <algorithm>

using namespace std;

#if HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_MPI
#include <mpi.h>
#endif


#define DENSE_CPU 0
#define DENSE_GPU 1
#define SPARSE_CPU 2

/// The neighbor_fuct value below which we consider
/// the impact zero for a given node in the map
#define NEIGHBOR_THRESHOLD 0.05

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

// Default parameters
#define N_EPOCH 10
#define N_SOM_X 50
#define N_SOM_Y 50
#define KERNEL_TYPE 0
#define SNAPSHOTS 0

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    int rank = 0;
       int nProcs = 1;

   #ifdef HAVE_MPI
       ///
       /// MPI init
       ///
       MPI_Init(&argc, &argv);
       MPI_Comm_rank(MPI_COMM_WORLD, &rank);
       MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
       MPI_Barrier(MPI_COMM_WORLD);
   #endif

       unsigned int nEpoch = 0;
       unsigned int nSomX = 0;
       unsigned int nSomY = 0;
       unsigned int kernelType = 0;
       string mapType;
       unsigned int radius0 = 0;
       unsigned int radiusN = 0;
       string radiusCooling;
       float scale0 = 0.0;
       float scaleN = 0.0;
       string scaleCooling;
       unsigned int snapshots = 0;
       string inFilename;
       string initialCodebookFilename;
       string outPrefix;

       if (rank==0) {
           processCommandLine(argc, argv, &inFilename, &outPrefix,
                              &nEpoch, &radius0, &radiusN, &radiusCooling,
                              &scale0, &scaleN, &scaleCooling,
                              &nSomX, &nSomY,
                              &kernelType, &mapType, &snapshots,
                              &initialCodebookFilename);
   #ifndef CUDA
           if (kernelType == DENSE_GPU) {
               cerr << "Somoclu was compile without GPU support!\n";
               my_abort(1);
           }
   #endif
       }
   #ifdef HAVE_MPI
       MPI_Bcast(&nEpoch, 1, MPI_INT, 0, MPI_COMM_WORLD);
       MPI_Bcast(&radius0, 1, MPI_INT, 0, MPI_COMM_WORLD);
       MPI_Bcast(&nSomX, 1, MPI_INT, 0, MPI_COMM_WORLD);
       MPI_Bcast(&nSomY, 1, MPI_INT, 0, MPI_COMM_WORLD);
       MPI_Bcast(&kernelType, 1, MPI_INT, 0, MPI_COMM_WORLD);

       char *inFilenameCStr = new char[255];
       if (rank == 0) {
           strcpy(inFilenameCStr,inFilename.c_str());
       }
       MPI_Bcast(inFilenameCStr, 255, MPI_CHAR, 0, MPI_COMM_WORLD);
       inFilename = inFilenameCStr;

       char *mapTypeCStr = new char[255];
       if (rank == 0) {
           strcpy(mapTypeCStr,mapType.c_str());
       }
       MPI_Bcast(mapTypeCStr, 255, MPI_CHAR, 0, MPI_COMM_WORLD);
       mapType = mapTypeCStr;

       double profile_time = MPI_Wtime();
   #endif

       float * dataRoot = NULL;
       unsigned int nDimensions = 0;
       unsigned int nVectors = 0;
       if(rank == 0 ) {
           if (kernelType == DENSE_CPU || kernelType == DENSE_GPU) {
               dataRoot = readMatrix(inFilename, nVectors, nDimensions);
           } else {
               readSparseMatrixDimensions(inFilename, nVectors, nDimensions);
           }
       }
   #ifdef HAVE_MPI
       MPI_Barrier(MPI_COMM_WORLD);
       MPI_Bcast(&nVectors, 1, MPI_INT, 0, MPI_COMM_WORLD);
       MPI_Bcast(&nDimensions, 1, MPI_INT, 0, MPI_COMM_WORLD);
   #endif
       unsigned int nVectorsPerRank = ceil(nVectors / (1.0*nProcs));

       // Allocate a buffer on each node
       float* data = NULL;
       svm_node **sparseData;
       sparseData = NULL;

       if (kernelType == DENSE_CPU || kernelType == DENSE_GPU) {
   #ifdef HAVE_MPI
           // Dispatch a portion of the input data to each node
           data = new float[nVectorsPerRank*nDimensions];
           MPI_Scatter(dataRoot, nVectorsPerRank*nDimensions, MPI_FLOAT,
                       data, nVectorsPerRank*nDimensions, MPI_FLOAT,
                       0, MPI_COMM_WORLD);
   #else
           data = dataRoot;
   #endif
       } else {
           int currentRankProcessed = 0;
           while (currentRankProcessed < nProcs) {
               if (rank == currentRankProcessed) {
                   sparseData=readSparseMatrixChunk(inFilename, nVectors, nVectorsPerRank,
                                                    rank*nVectorsPerRank);
               }
               currentRankProcessed++;
   #ifdef HAVE_MPI
               MPI_Barrier(MPI_COMM_WORLD);
   #endif
           }
       }

       if(rank == 0) {
           // No need for root data any more if compiled with MPI
   #ifdef HAVE_MPI
           if (kernelType == DENSE_CPU || kernelType == DENSE_GPU) {
               delete [] dataRoot;
           }
   #endif
           cout << "nVectors: " << nVectors << " ";
           cout << "nVectorsPerRank: " << nVectorsPerRank << " ";
           cout << "nDimensions: " << nDimensions << " ";
           cout << endl;
       }

   #ifdef CUDA
       if (kernelType == DENSE_GPU) {
           setDevice(rank, nProcs);
           initializeGpu(data, nVectorsPerRank, nDimensions, nSomX, nSomY);
       }
   #endif

   #ifdef HAVE_MPI
       MPI_Barrier(MPI_COMM_WORLD);
   #endif
       // TRAINING
       train(rank, data, sparseData, nSomX, nSomY,
             nDimensions, nVectors, nVectorsPerRank,
             nEpoch, radius0, radiusN, radiusCooling,
             scale0, scaleN, scaleCooling,
             outPrefix, snapshots, kernelType, mapType,
             initialCodebookFilename);

       if (kernelType == DENSE_CPU || kernelType == DENSE_GPU) {
           delete [] data;
       } else {
           delete [] sparseData;
       }
   #ifdef HAVE_MPI
       profile_time = MPI_Wtime() - profile_time;
       if (rank == 0) {
           cerr << "Total Execution Time: " << profile_time << endl;
       }
   #endif
   #ifdef CUDA
       if (kernelType == DENSE_GPU) {
           freeGpu();
       }
   #endif
   #ifdef HAVE_MPI
       MPI_Finalize();
   #endif

    MainWindow w;
    w.show();

    return a.exec();
}
