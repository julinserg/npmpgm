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

#include <cstdlib>
#include <cmath>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string.h>
#include <iomanip>
using namespace std;

float euclideanDistanceOnPlanarMap(const unsigned int som_x, const unsigned int som_y, const unsigned int x, const unsigned int y) {
    unsigned int x1 = std::min(som_x, x);
    unsigned int y1 = std::min(som_y, y);
    unsigned int x2 = std::max(som_x, x);
    unsigned int y2 = std::max(som_y, y);
    unsigned int xdist = x2-x1;
    unsigned int ydist = y2-y1;
    return sqrt(float(xdist*xdist+ydist*ydist));
}

float euclideanDistanceOnToroidMap(const unsigned int som_x, const unsigned int som_y, const unsigned int x, const unsigned int y, const unsigned int nSomX, const unsigned int nSomY) {
    unsigned int x1 = std::min(som_x, x);
    unsigned int y1 = std::min(som_y, y);
    unsigned int x2 = std::max(som_x, x);
    unsigned int y2 = std::max(som_y, y);
    unsigned int xdist = std::min(x2-x1, x1+nSomX-x2);
    unsigned int ydist = std::min(y2-y1, y1+nSomY-y2);

    //unsigned int xdist = std::min(std::abs(x1-x2), std::abs(x1+nSomX-x2));
    //unsigned int ydist = std::min(std::abs(y1-y2), std::abs(y1+nSomY-y2));
    return sqrt(float(xdist*xdist+ydist*ydist));
}

float gaussianNeighborhood(float distance, float radius, float stddevs)
{
    float norm = (2 * (radius + 1)*(radius + 1)) / (stddevs*stddevs);
    return exp((-(float) distance * distance) / norm);
}

float getWeight(float distance, float radius, float scaling)
{
    if (distance <= radius)
    {
        return scaling * gaussianNeighborhood(distance, radius, 2);
    }
    else
    {
        return 0.0;
    }
}

/** Shut down MPI cleanly if something goes wrong
 * @param err - error code to print
 */
void my_abort(int err)
{
    cerr << "Aborted\n";
#ifdef HAVE_MPI
    MPI_Abort(MPI_COMM_WORLD, err);
#else
#endif
}

/** Save a SOM codebook
 * @param cbFileName - name of the file to save
 * @param codebook - the codebook to save
 * @param nSomX - dimensions of SOM map in the x direction
 * @param nSomY - dimensions of SOM map in the y direction
 * @param nDimensions - dimensions of a data instance
 */
int saveCodebook(string cbFilename, float *codebook, unsigned int nSomX, unsigned int nSomY, unsigned int nDimensions)
{
    FILE* file = fopen(cbFilename.c_str(), "wt");
    cout << "    Saving Codebook " << cbFilename << endl;
    fprintf(file, "%%%d %d\n", nSomY, nSomX);
    fprintf(file, "%%%d\n", nDimensions);
    if (file!=0) {
        for (unsigned int som_y = 0; som_y < nSomY; som_y++) {
            for (unsigned int som_x = 0; som_x < nSomX; som_x++) {
                for (unsigned int d = 0; d < nDimensions; d++) {
                    fprintf(file, "%0.10f ", codebook[som_y*nSomX*nDimensions+som_x*nDimensions+d]);
                }
                fprintf(file, "\n");
            }
        }
        fclose(file);
        return 0;
    } else {
        return 1;
    }
}

/** Save best matching units
 * @param filename - name of the file to save
 * @param bmus - the best matching units to save
 * @param nSomX - dimensions of SOM map in the x direction
 * @param nSomY - dimensions of SOM map in the y direction
 * @param nVectors - the number of vectors
 */
int saveBmus(string filename, int *bmus, unsigned int nSomX, unsigned int nSomY, unsigned int nVectors)
{
    FILE* file = fopen(filename.c_str(), "wt");
    cout << "    Saving best matching units " << filename << endl;
    fprintf(file, "%%%d %d\n", nSomY, nSomX);
    fprintf(file, "%%%d\n", nVectors);
    if (file!=0) {
        for (unsigned int i = 0; i < nVectors; ++i) {
            // ESOM Tools swaps x and y!
            fprintf(file, "%d %d %d\n", i, bmus[2*i+1], bmus[2*i]);
        }
        fclose(file);
        return 0;
    } else {
        return 1;
    }
}


/** Save u-matrix
 * @param fname
 * @param codebook - the codebook to save
 * @param nSomX - dimensions of SOM map in the x direction
 * @param nSomY - dimensions of SOM map in the y direction
 * @param nDimensions - dimensions of a data instance
 */

int saveUMatrix(string fname, float *uMatrix, unsigned int nSomX,
             unsigned int nSomY)
{

    FILE* fp = fopen(fname.c_str(), "wt");
    fprintf(fp, "%%");
    fprintf(fp, "%d %d", nSomY, nSomX);
    fprintf(fp, "\n");
    if (fp != 0) {
        for (unsigned int som_y1 = 0; som_y1 < nSomY; som_y1++) {
            for (unsigned int som_x1 = 0; som_x1 < nSomX; som_x1++) {
                fprintf(fp, " %f", uMatrix[som_y1*nSomX+som_x1]);
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
        return 0;
    }
    else
        return -2;
}


void getMatrixDimensions(string inFilename, unsigned int &nRows, unsigned int &nColumns)
{
    ifstream file;
    file.open(inFilename.c_str());
    if (file.is_open()) {
        string line;
        float tmp;
        while(getline(file,line)) {
            if (line.substr(0,1) == "#") {
                continue;
            }
            std::istringstream iss(line);
            if (nRows == 0) {
                while (iss >> tmp) {
                    nColumns++;
                }
            }
            nRows++;
        }
        file.close();
    } else {
        std::cerr << "Input file could not be opened!\n";
        my_abort(-1);
    }
}

unsigned int *readLrnHeader(string inFilename, unsigned int &nRows, unsigned int &nColumns)
{
    ifstream file;
    file.open(inFilename.c_str());
    string line;
    unsigned int currentColumn = 0;
    unsigned int nAllColumns = 0;
    unsigned int *columnMap = NULL;
    while(getline(file,line)) {
        if (line.substr(0,1) == "#") {
            continue;
        }
        if (line.substr(0,1) == "%") {
            std::istringstream iss(line.substr(1,line.length()));
            if (nRows == 0) {
                iss >> nRows;
            } else if (nAllColumns == 0) {
                iss >> nAllColumns;
            } else if (columnMap == NULL) {
                columnMap = new unsigned int[nAllColumns];
                unsigned int itmp = 0;
                currentColumn = 0;
                while (iss >> itmp) {
                    columnMap[currentColumn++] = itmp;
                    if (itmp == 1) {
                        ++nColumns;
                    }
                }
            } else {
                break;
            }
        }
    }
    file.close();
    return columnMap;
}

unsigned int *readWtsHeader(string inFilename, unsigned int &nRows, unsigned int &nColumns)
{
    ifstream file;
    file.open(inFilename.c_str());
    string line;
    while(getline(file,line)) {
        if (line.substr(0,1) == "#") {
            continue;
        }
        if (line.substr(0,1) == "%") {
            std::istringstream iss(line.substr(1,line.length()));
            if (nRows == 0) {
                iss >> nRows;
                unsigned int nSomY = 0;
                iss >> nSomY;
                nRows = nRows*nSomY;
            } else if (nColumns == 0) {
                iss >> nColumns;
            } else {
                break;
            }
        }
    }
    file.close();
    unsigned int *columnMap = new unsigned int[nColumns];
    for (unsigned int i = 0; i < nColumns; ++i) {
        columnMap[i] = 1;
    }
    return columnMap;
}

/** Reads a matrix
 * @param inFilename
 * @param nRows - returns the number of rows
 * @param nColumns - returns the number of columns
 * @return the matrix
 */
float *readMatrix(string inFilename, unsigned int &nRows, unsigned int &nColumns)
{
    float *data = NULL;
    unsigned int *columnMap = NULL;
    if (inFilename.compare(inFilename.size()-3, 3, "lrn") == 0) {
        columnMap = readLrnHeader(inFilename, nRows, nColumns);
    } else if (inFilename.compare(inFilename.size()-3, 3, "wts") == 0) {
        columnMap = readWtsHeader(inFilename, nRows, nColumns);
    } else {
        getMatrixDimensions(inFilename, nRows, nColumns);
        columnMap = new unsigned int[nColumns];
        for (unsigned int i = 0; i < nColumns; ++i) {
            columnMap[i] = 1;
        }
    }
    ifstream file;
    file.open(inFilename.c_str());
    string line;
    float tmp;
    unsigned int j = 0;
    unsigned int currentColumn = 0;

    while(getline(file,line)) {
        if ( (line.substr(0,1) == "#") | (line.substr(0,1) == "%") ) {
            continue;
        }
        if (data == NULL) {
            data = new float[nRows*nColumns];
        }
        std::istringstream iss(line);
        currentColumn = 0;
        while (iss >> tmp) {
            if (columnMap[currentColumn++] != 1) {
                continue;
            }
            data[j++] = tmp;
        }
    }
    file.close();
    delete [] columnMap;
    return data;
}

void readSparseMatrixDimensions(string filename, unsigned int &nRows,
                                unsigned int &nColumns) {
    ifstream file;
    file.open(filename.c_str());
    if (file.is_open()) {
        string line;
        int max_index=-1;
        while(getline(file,line)) {
            if (line.substr(0,1) == "#") {
                continue;
            }
            stringstream linestream(line);
            string value;
            int dummy_index;
            while(getline(linestream,value,' ')) {
                int separator=value.find(":");
                istringstream myStream(value.substr(0,separator));
                myStream >> dummy_index;
                if(dummy_index > max_index) {
                    max_index = dummy_index;
                }
            }
            ++nRows;
        }
        nColumns=max_index+1;
        file.close();
    } else {
        std::cerr << "Input file could not be opened!\n";
        my_abort(-1);
    }
}

svm_node** readSparseMatrixChunk(string filename, unsigned int nRows,
                                 unsigned int nRowsToRead,
                                 unsigned int rowOffset) {
    ifstream file;
    file.open(filename.c_str());
    string line;
    for (unsigned int i=0; i<rowOffset; i++) {
        getline(file, line);
    }
    if (rowOffset+nRowsToRead >= nRows) {
        nRowsToRead = nRows-rowOffset;
    }
    svm_node **x_matrix = new svm_node *[nRowsToRead];
    for(unsigned int i=0; i<nRowsToRead; i++) {
        getline(file, line);
        if (line.substr(0,1) == "#") {
            --i;
            continue;
        }
        stringstream tmplinestream(line);
        string value;
        int elements = 0;
        while(getline(tmplinestream,value,' ')) {
            elements++;
        }
        elements++; // To account for the closing dummy node in the row
        x_matrix[i] = new svm_node[elements];
        stringstream linestream(line);
        int j=0;
        while(getline(linestream,value,' ')) {
            int separator=value.find(":");
            istringstream myStream(value.substr(0,separator));
            myStream >> x_matrix[i][j].index;
            istringstream myStream2(value.substr(separator+1));
            myStream2 >> x_matrix[i][j].value;
            j++;
        }
        x_matrix[i][j].index = -1;
    }
    file.close();
    return x_matrix;
}

/** Distance b/w a feature vector and a weight vector
 * = Euclidean
 * @param som_y
 * @param som_x
 * @param r - row number in the input feature file
  */

float get_distance(float* codebook, svm_node **sparseData,
                   unsigned int som_y, unsigned int som_x, unsigned int nSomX,
                   unsigned int nDimensions, unsigned int r)
{
    float distance = 0.0f;
    unsigned int j=0;
    for ( unsigned int d=0; d < nDimensions; d++ ) {
        if ( (int) d == sparseData[r][j].index ) {
            distance += (codebook[som_y*nSomX*nDimensions+som_x*nDimensions+d]-
                         sparseData[r][j].value) *
                        (codebook[som_y*nSomX*nDimensions+som_x*nDimensions+d]-
                         sparseData[r][j].value);
            ++j;
        } else {
            distance += codebook[som_y*nSomX*nDimensions+som_x*nDimensions+d]*
                        codebook[som_y*nSomX*nDimensions+som_x*nDimensions+d];
        }
    }
    return distance;
}

/** Get node coords for the best matching unit (BMU)
 * @param coords - BMU coords
 * @param n - row num in the input feature file
 */
void get_bmu_coord(float* codebook, svm_node **sparseData,
                   unsigned int nSomY, unsigned int nSomX,
                   unsigned int nDimensions, int* coords, unsigned int n)
{
    float mindist = 9999.99;
    float dist = 0.0f;

    /// Check nSomX * nSomY nodes one by one and compute the distance
    /// D(W_K, Fvec) and get the mindist and get the coords for the BMU.
    ///
    for (unsigned int som_y = 0; som_y < nSomY; som_y++) {
        for (unsigned int som_x = 0; som_x < nSomX; som_x++) {
            dist = get_distance(codebook, sparseData, som_y, som_x, nSomX,
                                nDimensions, n);
            if (dist < mindist) {
                mindist = dist;
                coords[0] = som_x;
                coords[1] = som_y;
            }
        }
    }
}

void trainOneEpochSparseCPU(int itask, svm_node **sparseData, float *numerator,
                            float *denominator, float *codebook,
                            unsigned int nSomX, unsigned int nSomY,
                            unsigned int nDimensions, unsigned int nVectors,
                            unsigned int nVectorsPerRank, float radius,
                            float scale, string mapType, int *globalBmus)
{
    int p1[2] = {0, 0};
    int *bmus = new int[nVectorsPerRank*2];
#ifdef _OPENMP
    #pragma omp parallel default(shared) private(p1)
#endif
    {
#ifdef _OPENMP
        #pragma omp for
#endif
#ifdef _WIN32
      for (int n = 0; n < nVectorsPerRank; n++) {
#else
      for (unsigned int n = 0; n < nVectorsPerRank; n++) {
#endif
            if (itask*nVectorsPerRank+n<nVectors) {
                /// get the best matching unit
                get_bmu_coord(codebook, sparseData, nSomY, nSomX,
                              nDimensions, p1, n);
                bmus[2*n] = p1[0]; bmus[2*n+1] = p1[1];
            }
        }
    }

    float *localNumerator = new float[nSomY*nSomX*nDimensions];
    float *localDenominator = new float[nSomY*nSomX];
#ifdef _OPENMP
    #pragma omp parallel default(shared)
#endif
    {
#ifdef _OPENMP
        #pragma omp for
#endif
#ifdef _WIN32
      for (int som_y = 0; som_y < nSomY; som_y++) {
#else
      for (unsigned int som_y = 0; som_y < nSomY; som_y++) {
#endif
            for (unsigned int som_x = 0; som_x < nSomX; som_x++) {
                localDenominator[som_y*nSomX + som_x] = 0.0;
                for (unsigned int d = 0; d < nDimensions; d++)
                    localNumerator[som_y*nSomX*nDimensions + som_x*nDimensions + d] = 0.0;
            }
        }

    /// Accumulate denoms and numers
#ifdef _OPENMP
        #pragma omp for
#endif
#ifdef _WIN32
      for (int som_y = 0; som_y < nSomY; som_y++) {
#else
      for (unsigned int som_y = 0; som_y < nSomY; som_y++) {
#endif

            for (unsigned int som_x = 0; som_x < nSomX; som_x++) {
                for (unsigned int n = 0; n < nVectorsPerRank; n++) {
                    if (itask*nVectorsPerRank+n<nVectors) {
                        float dist = 0.0f;
                        if (mapType == "planar") {
                            dist = euclideanDistanceOnPlanarMap(som_x, som_y, bmus[2*n], bmus[2*n+1]);
                        } else if (mapType == "toroid") {
                            dist = euclideanDistanceOnToroidMap(som_x, som_y, bmus[2*n], bmus[2*n+1], nSomX, nSomY);
                        }
                        float neighbor_fuct = getWeight(dist, radius, scale);
                        unsigned int j=0;
                        while ( sparseData[n][j].index!=-1 ) {
                            localNumerator[som_y*nSomX*nDimensions +
                                           som_x*nDimensions +
                                           sparseData[n][j].index] +=
                                               1.0f * neighbor_fuct * sparseData[n][j].value;
                            j++;
                        }
                        localDenominator[som_y*nSomX + som_x] += neighbor_fuct;
                    }
                }
            }
        }
    }
#ifdef HAVE_MPI
    MPI_Reduce(localNumerator, numerator,
               nSomY*nSomX*nDimensions, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(localDenominator, denominator,
               nSomY*nSomX, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Gather(bmus, nVectorsPerRank*2, MPI_INT, globalBmus, nVectorsPerRank*2, MPI_INT, 0, MPI_COMM_WORLD);
#else
    for (unsigned int i=0; i < nSomY*nSomX*nDimensions; ++i) {
        numerator[i] = localNumerator[i];
    }
    for (unsigned int i=0; i < nSomY*nSomX; ++i) {
        denominator[i] = localDenominator[i];
    }
    for (unsigned int i=0; i < 2*nVectorsPerRank; ++i) {
      globalBmus[i]=bmus[i];
    }
#endif
    delete [] bmus;
    delete [] localNumerator;
    delete [] localDenominator;
}







/** Distance b/w a feature vector and a weight vector
 * = Euclidean
 * @param som_y
 * @param som_x
 * @param r - row number in the input feature file
  */

float get_distance(float* codebook, float* data,
                   unsigned int som_y, unsigned int som_x, unsigned int nSomX,
                   unsigned int nDimensions, unsigned int r)
{
    float distance = 0.0f;
    for (unsigned int d = 0; d < nDimensions; d++)
        distance += (codebook[som_y*nSomX*nDimensions+som_x*nDimensions+d] -
                     *(data + r*nDimensions + d))
                    *
                    (codebook[som_y*nSomX*nDimensions+som_x*nDimensions+d] -
                     *(data + r*nDimensions + d));
    return distance;
}

/** Get node coords for the best matching unit (BMU)
 * @param coords - BMU coords
 * @param n - row num in the input feature file
 */
void get_bmu_coord(float* codebook, float* data,
                   unsigned int nSomY, unsigned int nSomX,
                   unsigned int nDimensions, unsigned int* coords, unsigned int n)
{
    float mindist = 9999.99;
    float dist = 0.0f;

    /// Check nSomX * nSomY nodes one by one and compute the distance
    /// D(W_K, Fvec) and get the mindist and get the coords for the BMU.
    ///
    for (unsigned int som_y = 0; som_y < nSomY; som_y++) {
        for (unsigned int som_x = 0; som_x < nSomX; som_x++) {
            dist = get_distance(codebook, data, som_y, som_x, nSomX,
                                nDimensions, n);
            if (dist < mindist) {
                mindist = dist;
                coords[0] = som_x;
                coords[1] = som_y;
            }
        }
    }
}

void trainOneEpochDenseCPU(int itask, float *data, float *numerator,
                           float *denominator, float *codebook,
                           unsigned int nSomX, unsigned int nSomY,
                           unsigned int nDimensions, unsigned int nVectors,
                           unsigned int nVectorsPerRank, float radius,
                           float scale, string mapType, int *globalBmus)
{
    unsigned int p1[2] = {0, 0};
    unsigned int *bmus = new unsigned int[nVectorsPerRank*2];
#ifdef _OPENMP
    #pragma omp parallel default(shared) private(p1)
#endif
    {
#ifdef _OPENMP
        #pragma omp for
#endif
#ifdef _WIN32
        for (int n = 0; n < nVectorsPerRank; n++) {
#else
        for (unsigned int n = 0; n < nVectorsPerRank; n++) {
#endif
            if (itask*nVectorsPerRank+n<nVectors) {
                /// get the best matching unit
                get_bmu_coord(codebook, data, nSomY, nSomX,
                              nDimensions, p1, n);
                bmus[2*n] = p1[0]; bmus[2*n+1] = p1[1];
              }
        }
    }

    float *localNumerator = new float[nSomY*nSomX*nDimensions];
    float *localDenominator = new float[nSomY*nSomX];
#ifdef _OPENMP
    #pragma omp parallel default(shared)
#endif
    {
#ifdef _OPENMP
        #pragma omp for
#endif
#ifdef _WIN32
      for (int som_y = 0; som_y < nSomY; som_y++) {
#else
      for (unsigned int som_y = 0; som_y < nSomY; som_y++) {
#endif
            for (unsigned int som_x = 0; som_x < nSomX; som_x++) {
                localDenominator[som_y*nSomX + som_x] = 0.0;
                for (unsigned int d = 0; d < nDimensions; d++)
                    localNumerator[som_y*nSomX*nDimensions + som_x*nDimensions + d] = 0.0;
            }
        }
        /// Accumulate denoms and numers
#ifdef _OPENMP
        #pragma omp for
#endif
#ifdef _WIN32
      for (int som_y = 0; som_y < nSomY; som_y++) {
#else
      for (unsigned int som_y = 0; som_y < nSomY; som_y++) {
#endif
            for (unsigned int som_x = 0; som_x < nSomX; som_x++) {
                for (unsigned int n = 0; n < nVectorsPerRank; n++) {
                    if (itask*nVectorsPerRank+n<nVectors) {
                        float dist = 0.0f;
                        if (mapType == "planar") {
                            dist = euclideanDistanceOnPlanarMap(som_x, som_y, bmus[2*n], bmus[2*n+1]);
                        } else if (mapType == "toroid") {
                            dist = euclideanDistanceOnToroidMap(som_x, som_y, bmus[2*n], bmus[2*n+1], nSomX, nSomY);
                        }
                        float neighbor_fuct = getWeight(dist, radius, scale);

                        for (unsigned int d = 0; d < nDimensions; d++) {
                            localNumerator[som_y*nSomX*nDimensions + som_x*nDimensions + d] +=
                                1.0f * neighbor_fuct
                                * (*(data + n*nDimensions + d));
                        }
                        localDenominator[som_y*nSomX + som_x] += neighbor_fuct;
                    }
                }
            }
        }
    }
#ifdef HAVE_MPI
    MPI_Reduce(localNumerator, numerator,
               nSomY*nSomX*nDimensions, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(localDenominator, denominator,
               nSomY*nSomX, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Gather(bmus, nVectorsPerRank*2, MPI_INT, globalBmus, nVectorsPerRank*2, MPI_INT, 0, MPI_COMM_WORLD);

#else
    for (unsigned int i=0; i < nSomY*nSomX*nDimensions; ++i) {
        numerator[i] = localNumerator[i];
    }
    for (unsigned int i=0; i < nSomY*nSomX; ++i) {
        denominator[i] = localDenominator[i];
    }
    for (unsigned int i=0; i < 2*nVectorsPerRank; ++i) {
      globalBmus[i]=bmus[i];
    }
#endif
    delete [] bmus;
    delete [] localNumerator;
    delete [] localDenominator;
}

float* get_wvec(float *codebook, unsigned int som_y, unsigned int som_x,
                unsigned int nSomX, unsigned int nDimensions)
{
    float* wvec = new float[nDimensions];
    for (unsigned int d = 0; d < nDimensions; d++)
        wvec[d] = codebook[som_y*nSomX*nDimensions+som_x*nDimensions+d]; /// CAUTION: (y,x) order
    return wvec;
}

/** Euclidean distance between vec1 and vec2
 * @param vec1
 * @param vec2
 * @param nDimensions
 * @return distance
 */

float get_distance(const float* vec1, const float* vec2,
                   unsigned int nDimensions) {
    float distance = 0.0f;
    float x1 = 0.0f;
    float x2 = 0.0f;
    for (unsigned int d = 0; d < nDimensions; d++) {
        x1 = std::min(vec1[d], vec2[d]);
        x2 = std::max(vec1[d], vec2[d]);
        distance += std::abs(x1-x2)*std::abs(x1-x2);
    }
    return sqrt(distance);
}


/** Calculate U-matrix
 * @param codebook - the codebook
 * @param nSomX - dimensions of SOM map in the x direction
 * @param nSomY - dimensions of SOM map in the y direction
 * @param nDimensions - dimensions of a data instance
 */

float *calculateUMatrix(float *codebook, unsigned int nSomX,
             unsigned int nSomY, unsigned int nDimensions, string mapType)
{
    float *uMatrix = new float[nSomX*nSomY];
    float min_dist = 1.5f;
    for (unsigned int som_y1 = 0; som_y1 < nSomY; som_y1++) {
        for (unsigned int som_x1 = 0; som_x1 < nSomX; som_x1++) {
            float dist = 0.0f;
            unsigned int nodes_number = 0;

            for (unsigned int som_y2 = 0; som_y2 < nSomY; som_y2++) {
                for (unsigned int som_x2 = 0; som_x2 < nSomX; som_x2++) {

                    if (som_x1 == som_x2 && som_y1 == som_y2) continue;
                    float tmp = 0.0f;
                    if (mapType == "planar") {
                        tmp = euclideanDistanceOnPlanarMap(som_x1, som_y1, som_x2, som_y2);
                    } else if (mapType == "toroid") {
                        tmp = euclideanDistanceOnToroidMap(som_x1, som_y1, som_x2, som_y2, nSomX, nSomY);
                    }
                    if (tmp <= min_dist) {
                        nodes_number++;
                        float* vec1 = get_wvec(codebook, som_y1, som_x1, nSomX, nDimensions);
                        float* vec2 = get_wvec(codebook, som_y2, som_x2, nSomX, nDimensions);
                        dist += get_distance(vec1, vec2, nDimensions);
                        delete [] vec1;
                        delete [] vec2;
                    }
                }
            }
            dist /= (float)nodes_number;
            uMatrix[som_y1*nSomX+som_x1] = dist;
        }
    }
    return uMatrix;
}

float linearCooling(float start, float end, float nEpoch, float epoch) {
  float diff = (start - end) / (nEpoch-1);
  return start - (epoch * diff);
}

float exponentialCooling(float start, float end, float nEpoch, float epoch) {
  float diff = 0;
  if (end == 0.0)
  {
      diff = -log(0.1) / nEpoch;
  }
  else
  {
      diff = -log(end / start) / nEpoch;
  }
  return start * exp(-epoch * diff);
}



/** Initialize SOM codebook with random values
 * @param seed - random seed
 * @param codebook - the codebook to fill in
 * @param nSomX - dimensions of SOM map in the currentEpoch direction
 * @param nSomY - dimensions of SOM map in the y direction
 * @param nDimensions - dimensions of a data instance
 */

void initializeCodebook(unsigned int seed, float *codebook, unsigned int nSomX,
                        unsigned int nSomY, unsigned int nDimensions)
{
    ///
    /// Fill initial random weights
    ///
#ifdef HAVE_R
    GetRNGstate();
#else
    srand(seed);
#endif
    for (unsigned int som_y = 0; som_y < nSomY; som_y++) {
        for (unsigned int som_x = 0; som_x < nSomX; som_x++) {
            for (unsigned int d = 0; d < nDimensions; d++) {
#ifdef HAVE_R
                int w = 0xFFF & (int) (RAND_MAX*unif_rand());
#else
                int w = 0xFFF & rand();
#endif
                w -= 0x800;
                codebook[som_y*nSomX*nDimensions+som_x*nDimensions+d] = (float)w / 4096.0f;
            }
        }
    }
#ifdef HAVE_R
    PutRNGstate();
#endif
}

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
           unsigned int kernelType, string mapType){

    float N = (float)nEpoch;
    float *numerator;
    float *denominator;
    float scale = scale0;
    float radius = radius0;
    if (itask == 0) {
        numerator = new float[nSomY*nSomX*nDimensions];
        denominator = new float[nSomY*nSomX];
        for (unsigned int som_y = 0; som_y < nSomY; som_y++) {
            for (unsigned int som_x = 0; som_x < nSomX; som_x++) {
                denominator[som_y*nSomX + som_x] = 0.0;
                for (unsigned int d = 0; d < nDimensions; d++) {
                    numerator[som_y*nSomX*nDimensions + som_x*nDimensions + d] = 0.0;
                }
            }
        }

        if (radiusCooling == "linear") {
          radius = linearCooling(float(radius0), radiusN, N, currentEpoch);
        } else {
          radius = exponentialCooling(radius0, radiusN, N, currentEpoch);
        }
        if (scaleCooling == "linear") {
          scale = linearCooling(scale0, scaleN, N, currentEpoch);
        } else {
          scale = exponentialCooling(scale0, scaleN, N, currentEpoch);
        }
//        cout << "Epoch: " << currentEpoch << " Radius: " << radius << endl;
    }
#ifdef HAVE_MPI
    MPI_Bcast(&radius, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&scale, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(coreData.codebook, nSomY*nSomX*nDimensions, MPI_FLOAT,
              0, MPI_COMM_WORLD);
#endif

    /// 1. Each task fills localNumerator and localDenominator
    /// 2. MPI_reduce sums up each tasks localNumerator and localDenominator to the root's
    ///    numerator and denominator.
    switch (kernelType) {
    default:
    case DENSE_CPU:
        trainOneEpochDenseCPU(itask, data, numerator, denominator,
                              coreData.codebook, nSomX, nSomY, nDimensions,
                              nVectors, nVectorsPerRank, radius, scale,
                              mapType, coreData.globalBmus);
        break;
#ifdef CUDA
    case DENSE_GPU:
        trainOneEpochDenseGPU(itask, data, numerator, denominator,
                              coreData.codebook, nSomX, nSomY, nDimensions,
                              nVectors, nVectorsPerRank, radius, scale,
                              mapType, coreData.globalBmus);
        break;
#endif
    case SPARSE_CPU:
        trainOneEpochSparseCPU(itask, sparseData, numerator, denominator,
                               coreData.codebook, nSomX, nSomY, nDimensions,
                               nVectors, nVectorsPerRank, radius, scale,
                               mapType, coreData.globalBmus);
        break;
    }

    /// 3. Update codebook using numerator and denominator
#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if (itask == 0) {
        #pragma omp parallel for
#ifdef _WIN32
        for (int som_y = 0; som_y < nSomY; som_y++) {
#else
        for (unsigned int som_y = 0; som_y < nSomY; som_y++) {
#endif
            for (unsigned int som_x = 0; som_x < nSomX; som_x++) {
                float denom = denominator[som_y*nSomX + som_x];
                for (unsigned int d = 0; d < nDimensions; d++) {
                    float newWeight = numerator[som_y*nSomX*nDimensions
                                                + som_x*nDimensions + d] / denom;
                    if (newWeight > 0.0) {
                        coreData.codebook[som_y*nSomX*nDimensions+som_x*nDimensions+d] = newWeight;
                    }
                }
            }
        }
    }
    if (enableCalculatingUMatrix) {
        coreData.uMatrix = calculateUMatrix(coreData.codebook, nSomX, nSomY, nDimensions, mapType);
    }
    if (itask == 0) {
        delete [] numerator;
        delete [] denominator;
    }
    return coreData;
}

/** Main training loop
 * @param itask - number of work items
 * @param kv
 * @param ptr
 */

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
           string initialCodebookFilename)
{
    ///
    /// Codebook
    ///
    core_data coreData;
    coreData.codebook = new float[nSomY*nSomX*nDimensions];
    coreData.globalBmus = NULL;
    coreData.uMatrix = NULL;
    if (itask == 0) {
        coreData.globalBmus = new int[nVectorsPerRank*int(ceil(nVectors/(double)nVectorsPerRank))*2];

        if (initialCodebookFilename.empty()){
            initializeCodebook(0, coreData.codebook, nSomX, nSomY, nDimensions);
        } else {
            unsigned int nSomXY = 0;
            unsigned int tmpNDimensions = 0;
            delete [] coreData.codebook;
            coreData.codebook = readMatrix(initialCodebookFilename, nSomXY, tmpNDimensions);
            if (tmpNDimensions != nDimensions) {
                cerr << "Dimension of initial codebook does not match data!\n";
                my_abort(5);
            } else if (nSomXY / nSomY != nSomX) {
                cerr << "Dimension of initial codebook does not match specified SOM grid!\n";
                my_abort(6);
            }
            cout << "Read initial codebook: " << initialCodebookFilename << "\n";
        }
    }
    ///
    /// Parameters for SOM
    ///
    if (radius0 == 0) {
        unsigned int minDim = min(nSomX, nSomY);
        radius0 = minDim / 2.0f;              /// init radius for updating neighbors
    }
    if (radiusN == 0) {
        radiusN = 1;
    }
    if (scale0 == 0) {
      scale0 = 0.1;
    }

    unsigned int currentEpoch = 0;             /// 0...nEpoch-1

    ///
    /// Training
    ///
#ifdef HAVE_MPI
    double training_time = MPI_Wtime();
#endif

    while ( currentEpoch < nEpoch ) {

#ifdef HAVE_MPI
        double epoch_time = MPI_Wtime();
#endif

        coreData = trainOneEpoch(itask, data, sparseData,
                                 coreData, nEpoch, currentEpoch,
                                 snapshots > 0,
                                 nSomX, nSomY,
                                 nDimensions, nVectors,
                                 nVectorsPerRank,
                                 radius0, radiusN,
                                 radiusCooling,
                                 scale0, scaleN,
                                 scaleCooling,
                                 kernelType, mapType);

        if (snapshots > 0 && itask == 0) {
            cout << "Saving interim U-Matrix..." << endl;
            stringstream sstm;
            sstm << outPrefix << "." << currentEpoch + 1;
            saveUMatrix(sstm.str() + string(".umx"), coreData.uMatrix, nSomX, nSomY);
            if (snapshots == 2){
                saveBmus(sstm.str() + string(".bm"), coreData.globalBmus, nSomX, nSomY, nVectors);
                saveCodebook(sstm.str() + string(".wts"), coreData.codebook, nSomX, nSomY, nDimensions);
            }
        }
        currentEpoch++;
#ifdef HAVE_MPI
        if (itask == 0) {
            epoch_time = MPI_Wtime() - epoch_time;
            cerr << "Epoch Time: " << epoch_time << endl;
            if ( (currentEpoch != nEpoch) && (currentEpoch % (nEpoch/100+1) != 0) ){}
            else{
              float ratio  =  currentEpoch/(float)nEpoch;
              int   c      =  ratio * 50 + 1;
              cout << std::setw(7) << (int)(ratio*100) << "% [";
              for (int x=0; x<c; x++) cout << "=";
              for (int x=c; x<50; x++) cout << " ";
              cout << "]\n" << flush;
            }
        }
#endif
    }
#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    training_time = MPI_Wtime() - training_time;
    if (itask == 0) {
        cerr << "Total training Time: " << training_time << endl;
    }
#endif
    ///
    /// Save SOM map and u-mat
    ///
    if (itask == 0) {
        ///
        /// Save U-mat
        ///
        coreData.uMatrix = calculateUMatrix(coreData.codebook, nSomX, nSomY, nDimensions, mapType);
        int ret =  saveUMatrix(outPrefix + string(".umx"), coreData.uMatrix, nSomX, nSomY);
        if (ret < 0)
            cout << "    Failed to save u-matrix. !" << endl;
        else {
            cout << "    Done!" << endl;
        }
        saveBmus(outPrefix + string(".bm"), coreData.globalBmus, nSomX, nSomY, nVectors);
        ///
        /// Save codebook
        ///
        saveCodebook(outPrefix + string(".wts"), coreData.codebook, nSomX, nSomY, nDimensions);
    }
#ifdef HAVE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    delete [] coreData.codebook;
    delete [] coreData.globalBmus;
    delete [] coreData.uMatrix;
}



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
   m_nProcs = 1;
   m_nThread = 1;

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




