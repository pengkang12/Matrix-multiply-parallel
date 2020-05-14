#include <mkl.h>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
using namespace std;

// To Run: 
// mpicxx matrix-multiple-normal.c  -DMKL_ILP64 -m64 -I${MKLROOT}/include -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl
 
// mpirun -n 1 ./a.out 64

void matGene(double *A, int size, int actual_size){
    // actual size: the matrix we use may have a larger dimension than n * n
    for (int i = 0; i < actual_size; i++){
        for (int j = 0; j < actual_size; j++){
            if(i < size && j < size) A[i * actual_size + j] = 1; //A[i][j]
            else A[i * actual_size + j] = 0;
        }
    }
}

void matMulti(double *A, double *B, double *C, int m, int n, int p){
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, p, 1, A, p, B, n, 0, C, n);
    return ;
}

void print_matrix(char *name, double* A, int row, int column, int n){
	return;
	cout << name << endl;
       	for(int i = 0; i < row; i++){
        	for(int j = 0; j < column; j++)
                    cout << A[i * n + j] << " ";
                cout << endl;
        }       
 
}

int main(int argc, char *argv[]){
    // Only Deal With Square Matrixs

    // Calculate Parameters Definition
    int n = atoi(argv[1]); // matrix dimension
    // int beginRow, endRow; // the range of rows calculating in certain process
    double beginTime, endTime; // time record
    srand(time(NULL));

    double *A = NULL;
    double *B = NULL;
    double *C = NULL;

    // Prepare data
    A = new double[n * n + 2];
    B = new double[n * n + 2];
    C = new double[n * n + 2];
    int saveN = n;
    matGene(A, saveN, n);
    matGene(B, saveN, n);

    // Calculate C[i][j] & Time
    std::clock_t    start, end;

    start = std::clock();
    matMulti(A, B, C, n, n, n);
    end = std::clock();
 
    cout << "Time: " << (end- start) / (double)(CLOCKS_PER_SEC ) << endl;

    // Output
    print_matrix("A", A, saveN, saveN, n);
    print_matrix("B", B, saveN, saveN, n);
    print_matrix("C", C, saveN, saveN, n);
    
    delete[] A;
    delete[] B;
    delete[] C;
    return 0;
}
