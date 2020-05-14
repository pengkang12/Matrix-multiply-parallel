#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
using namespace std;

// To Run: 
// mpicxx matrix_multi.cpp  ../../scalapack-2.1.0/libscalapack.a -llapack -lblas


// mpirun -n 4 ./a.out 64

// Cannon Algorithm + Block Multiplication
// Main process: process 0, data distribution & collect calculate results, no calculation
// Others: calculate for a block of A * B
extern "C" {
	void pcgemm_(char, char, int*, int*, int*, float*, float*, int*, int*, int*, float*, int*, int*, int*, float*, float*, int*, int *, int*);
	int* IA, *JA,* IB, *JB, *IC, *JC, *M, *N, *K, *DESC,* DESA, *DESB; 
        float* ALPHA, *BETA, *AA, *BB, *CC;
        char TRANSA, TRANSB;
	    
} 
void matGene(float *A, int size, int actual_size){
    // actual size: the matrix we use may have a larger dimension than n * n
    for (int i = 0; i < actual_size; i++){
        for (int j = 0; j < actual_size; j++){
            if(i < size && j < size) A[i * actual_size + j] = 1; //A[i][j]
            else A[i * actual_size + j] = 0;
        }
    }
}


void print_matrix(char *name, float* A, int row, int column, int n){
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
    float beginTime, endTime; // time record
    srand(time(NULL));

    // MPI Common Head
    int my_rank = 0, comm_sz = 0;
    int a = sqrt(comm_sz);
        if(my_rank == 0 && comm_sz != a * a){
            cout << "Not Full Square" << endl;
            return 0;
        }

        int saveN = n;
        // must equal scatter: actual n is bigger than input
        if(n % a != 0){
            n -= n % a;
            n += a;
        }   


    if (my_rank == 0){  
            // Prepare data
            cout << "n = " << n << endl;
        AA = new float[n * n + 2];
        BB = new float[n * n + 2];
        CC = new float[n * n + 2];
 
//            matGene(AA, saveN, n);
 //           matGene(BB, saveN, n); 
            for(int ii = 0; ii < n; ii++){
                for(int jj = 0; jj < n; jj++){
//                    CC[ii * n + jj] = 0;
                }
            }  
    }

    return 0;         
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Status status;

    beginTime = MPI_Wtime();   
    float *A = NULL;
    float *B = NULL;
    float *C = NULL;

    if (comm_sz == 1){ // no parallel
    }

    else{ // parallel: main process collect the result and also involve in calculation

        int each_row = n / a;
        int beginRow, beginColumn;
            
        // Data generation
        for(int k = 0; k < a; k++){
            int begin_part = k * each_row;
            TRANSA='N', TRANSB='T';
	    
	    for(int i = 0; i < comm_sz; i++){
               	beginRow = (i / a) * each_row;
                beginColumn = (i % a) * each_row;
             
           	*M = each_row;
	        *N = each_row;
	        *K = each_row;
	        *IA = beginRow;
	        *JA = begin_part;
                *IB = begin_part;
	        *JB = beginColumn;
	        *IC = beginRow;
	        *JC = beginColumn;
               
    		//pcgemm_(TRANSA, TRANSB, M, N, K, ALPHA, AA, IA, JA, DESA, BB, IB,JB,DESB, BETA, CC, IC, JC, DESC);
            	print_matrix("C", CC, saveN, saveN, n);
            
	    } 
        }

        if (my_rank == 0){
            endTime = MPI_Wtime();
            // Output   
	    cout << "Time: " << endTime - beginTime << endl; 
	    print_matrix("A", AA, saveN, saveN, n);
            print_matrix("B", BB, saveN, saveN, n);
            print_matrix("C", CC, saveN, saveN, n);
            
        }
    }           
    MPI_Finalize();
    return 0;
}
