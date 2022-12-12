#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <vector>
#include "linearSystem.h"
#include <mpi.h>

/*
    mpicc main.c -o gauss & mpirun -np 3 gauss
*/

using namespace std;

int *pParallelPivotPos; // Number of rows selected as the pivot ones
int *pProcPivotIter; // Number of iterations, at which the processor
// rows were used as the pivot ones
int processRank, processNumber;

void loadSystem(vector<vector<double>> &matrxi, vector<double> &);

/*
    Função sequencial do metodo de gauss
    retorna o vetor resultado
*/
vector<double> sequencialGauss(vector <vector <double> > matrix, vector<double> rightVector);

/*
    Função paralela do metodo de gauss com MPI
    retorna o vetor resultado
*/
vector<double> parallelMPIGauss(vector <vector <double> > matrix, vector<double> rightVector);

/*
    Função paralela do metodo de gauss com OpenMP
    retorna o vetor resultado
*/
vector<double> parallelOMPGauss(vector <vector <double> > matrix, vector<double> rightVector);

int main(int argc, char **argv)
{
    int sistemSize; // Tamanho da matriz e dos vetores
    vector <vector <double>> matrix; // matriz do sistema linear
    vector <double> rightVec; // Right parts of the linear system
    vector <double> vecResult; // vetor resultado
    double *pProcRows; // Rows of the matrix A
    double *pProcVector; // Block of the vector b
    double *pProcResult; // Block of the vector x
    int Size; // Size of the matrix and vectors

    MPI_Init ( &argc, &argv );
    MPI_Comm_rank ( MPI_COMM_WORLD, &processRank );
    MPI_Comm_size ( MPI_COMM_WORLD, &processNumber );
    matrix = {
        {3,2,0,1},
        {9,8,-3,4},
        {-6,4,-8,0},
        {3,-8,3,-12}
    };
    rightVec = {3,6,-16,22};

    LinearSystem linSys("linearSys1.txt");

    if (processRank == 0)
    {
        cout << "Test" << endl;
        
        linSys.printSystem();
    }

    MPI_Finalize();

    return 0;
}

/*void ParallelGaussianElimination (double* pProcRows, double* pProcVector, int Size, int RowNum) 
{
    double MaxValue; // Value of the pivot element of thе process
    int PivotPos; // Position of the pivot row in the process stripe
    // Structure for the pivot row selection
    struct { double MaxValue; int processRank; } ProcPivot, Pivot;
    // pPivotRow is used for storing the pivot row and the corresponding
    // element of the vector b
    double* pPivotRow = new double [Size+1];
    // The iterations of the Gaussian elimination stage
    for (int i=0; i<Size; i++) 
    {
        // Calculating the local pivot row
        double MaxValue = 0;
        for (int j=0; j<RowNum; j++) 
        {
            if ((pProcPivotIter[j] == -1) &&
                (MaxValue < fabs(pProcRows[j*Size+i]))) 
            {
                MaxValue = fabs(pProcRows[j*Size+i]);
                PivotPos = j;
            }
        }

        ProcPivot.MaxValue = MaxValue;
        ProcPivot.processRank = processRank;
        // Finding the pivot process (process with the maximum value of MaxValue)
        MPI_Allreduce(&ProcPivot, &Pivot, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
        // Broadcasting the pivot row
        if ( processRank == Pivot.processRank )
        {
            pProcPivotIter[PivotPos]= i; //iteration number
            pParallelPivotPos[i]= pProcInd[processRank] + PivotPos;
        }

        MPI_Bcast(&pParallelPivotPos[i], 1, MPI_INT, Pivot.processRank, MPI_COMM_WORLD);
        if ( processRank == Pivot.processRank )
        {
        // Fill the pivot row
            for (int j=0; j<Size; j++) 
            {
                pPivotRow[j] = pProcRows[PivotPos*Size + j];
            }
            pPivotRow[Size] = pProcVector[PivotPos];
        }
        MPI_Bcast(pPivotRow, Size+1, MPI_DOUBLE, Pivot.processRank, MPI_COMM_WORLD);
        ParallelEliminateColumns(pProcRows, pProcVector, pPivotRow, Size,RowNum, i);
    }
}

void ParallelBackSubstitution (double* pProcRows, double* pProcVector, double* pProcResult, int Size, int RowNum) {
    int IterProcRank; // Rank of the process with the current pivot row
    int IterPivotPos; // Position of the pivot row of the process
    double IterResult; // Calculated value of the current unknown
    double val;
    
    // Iterations of the back substitution stage
    for (int i=Size-1; i>=0; i--) {
        // Calculating the rank of the process, which holds the pivot row
        FindBackPivotRow(pParallelPivotPos[i], Size, IterProcRank, IterPivotPos);
        // Calculating the unknown
        if (ProcRank == IterProcRank) 
        {
            IterResult = pProcVector[IterPivotPos]/pProcRows[IterPivotPos*Size+i];
            pProcResult[IterPivotPos] = IterResult;
        }
        // Broadcasting the value of the current unknown
        MPI_Bcast(&IterResult, 1, MPI_DOUBLE, IterProcRank, MPI_COMM_WORLD);
        // Updating the values of the vector b
        for (int j=0; j<RowNum; j++)
            if ( pProcPivotIter[j] < i ) 
            {
                val = pProcRows[j*Size + i] * IterResult;
                pProcVector[j]=pProcVector[j] - val;
            }
    }
}
*/

vector<double> sequencialGauss(vector<vector<double>> matrix, vector<double> rightVector)
{
    return vector<double>();
}

vector<double> parallelGauss(vector<vector<double>> matrix, vector<double> rightVector)
{
    return vector<double>();
}
