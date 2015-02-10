#include "qsomthread.h"
// Default parameters
#define N_EPOCH 10
#define N_SOM_X 50
#define N_SOM_Y 50
#define KERNEL_TYPE 0
#define SNAPSHOTS 0

#define DENSE_CPU 0
#define DENSE_GPU 1
#define SPARSE_CPU 2

QSOMThread::QSOMThread()
{
   m_nEpoch = N_EPOCH;
   m_nSomX = N_SOM_X;
   m_nSomY = N_SOM_Y;
   m_kernelType = KERNEL_TYPE;
   m_snapshots = SNAPSHOTS;
   m_mapType = "planar";
   m_radius0 = 0;
   m_radiusN = 0;
   m_radiusCooling = "linear";
   m_scale0 = 0.0;
   m_scaleN = 0.01;
   m_scaleCooling = "linear";
   m_rank = 0;
   m_nProcs = 4;
   m_nThread = 4;

}

QSOMThread::~QSOMThread()
{

}

void QSOMThread::run()
{
    exec();
}

 void QSOMThread::startTrain()
 {
     int argc = 1;
     QString strComLine = "-np %1";
     strComLine = strComLine.arg(m_nThread);
     char *param = strComLine.toLatin1().data();
     char **argv = &param;
    #ifdef HAVE_MPI
        ///
        /// MPI init
        ///
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &m_nProcs);
        MPI_Barrier(MPI_COMM_WORLD);
    #endif

        if (m_rank==0) {
    #ifndef CUDA
            if (m_kernelType == DENSE_GPU) {
                cerr << "Somoclu was compile without GPU support!\n";
            }
    #endif
        }
    #ifdef HAVE_MPI
        MPI_Bcast(&m_nEpoch, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&m_radius0, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&m_nSomX, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&m_nSomY, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&m_kernelType, 1, MPI_INT, 0, MPI_COMM_WORLD);

        char *inFilenameCStr = new char[255];
        if (m_rank == 0) {
            strcpy(inFilenameCStr,m_inFilename.c_str());
        }
        MPI_Bcast(inFilenameCStr, 255, MPI_CHAR, 0, MPI_COMM_WORLD);
        m_inFilename = inFilenameCStr;

        char *mapTypeCStr = new char[255];
        if (m_rank == 0) {
            strcpy(mapTypeCStr,m_mapType.c_str());
        }
        MPI_Bcast(mapTypeCStr, 255, MPI_CHAR, 0, MPI_COMM_WORLD);
        m_mapType = mapTypeCStr;

        double profile_time = MPI_Wtime();
    #endif

        float * dataRoot = NULL;
        unsigned int nDimensions = 0;
        unsigned int nVectors = 0;
        if(m_rank == 0 ) {
            if (m_kernelType == DENSE_CPU || m_kernelType == DENSE_GPU) {
                dataRoot = readMatrix(m_inFilename, nVectors, nDimensions);
            } else {
                readSparseMatrixDimensions(m_inFilename, nVectors, nDimensions);
            }
        }
    #ifdef HAVE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&nVectors, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&nDimensions, 1, MPI_INT, 0, MPI_COMM_WORLD);
    #endif
        unsigned int nVectorsPerRank = ceil(nVectors / (1.0*m_nProcs));

        // Allocate a buffer on each node
        float* data = NULL;
        svm_node **sparseData;
        sparseData = NULL;

        if (m_kernelType == DENSE_CPU || m_kernelType == DENSE_GPU) {
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
            while (currentRankProcessed < m_nProcs) {
                if (m_rank == currentRankProcessed) {
                    sparseData=readSparseMatrixChunk(m_inFilename, nVectors, nVectorsPerRank,
                                                     m_rank*nVectorsPerRank);
                }
                currentRankProcessed++;
    #ifdef HAVE_MPI
                MPI_Barrier(MPI_COMM_WORLD);
    #endif
            }
        }

        if(m_rank == 0) {
            // No need for root data any more if compiled with MPI
    #ifdef HAVE_MPI
            if (m_kernelType == DENSE_CPU || m_kernelType == DENSE_GPU) {
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
        train(m_rank, data, sparseData, m_nSomX, m_nSomY,
              nDimensions, nVectors, nVectorsPerRank,
              m_nEpoch, m_radius0, m_radiusN, m_radiusCooling,
              m_scale0, m_scaleN, m_scaleCooling,
              m_outPrefix, m_snapshots, m_kernelType, m_mapType,
              m_initialCodebookFilename);

        if (m_kernelType == DENSE_CPU || m_kernelType == DENSE_GPU) {
            delete [] data;
        } else {
            delete [] sparseData;
        }
    #ifdef HAVE_MPI
        profile_time = MPI_Wtime() - profile_time;
        if (m_rank == 0) {
            cerr << "Total Execution Time: " << profile_time << endl;
        }
    #endif
    #ifdef CUDA
        if (m_kernelType == DENSE_GPU) {
            freeGpu();
        }
    #endif
    #ifdef HAVE_MPI
        MPI_Finalize();
    #endif

 }




