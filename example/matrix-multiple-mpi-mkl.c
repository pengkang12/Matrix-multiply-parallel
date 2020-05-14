#include <mpi.h>
#include <mkl.h>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
using namespace std;

// To Run: 
// mpicxx matrix-multiple-normal.c  -DMKL_ILP64 -m64 -I${MKLROOT}/include -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl
 
// mpirun -n 4 ./a.out 64

// Cannon Algorithm + Block Multiplication
// Main process: process 0, data distribution & collect calculate results, no calculation
// Others: calculate for a block of A * B

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

    // MPI Common Head
    int my_rank = 0, comm_sz = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Status status;

    double *A = NULL;
    double *B = NULL;
    double *C = NULL;

    if (comm_sz == 1){ // no parallel
        // Prepare data
        A = new double[n * n + 2];
        B = new double[n * n + 2];
        C = new double[n * n + 2];
        int saveN = n;
        matGene(A, saveN, n);
        matGene(B, saveN, n);

        // Calculate C[i][j] & Time
        beginTime = MPI_Wtime();
        matMulti(A, B, C, n, n, n);
        endTime = MPI_Wtime();
        cout << "Time: " << endTime - beginTime << endl;

        // Output
        print_matrix("A", A, saveN, saveN, n);
        print_matrix("B", B, saveN, saveN, n);
        print_matrix("C", C, saveN, saveN, n);
        
        delete[] A;
        delete[] B;
        delete[] C;
    }

    else{ // parallel: main process collect the result and also involve in calculation

        int a = sqrt(comm_sz);
        if(my_rank == 0 && comm_sz != a * a){
            cout << "Not Full Square" << endl;
            MPI_Finalize();
            return 0;
        }



        int saveN = n;
        // must equal scatter: actual n is bigger than input
        if(n % a != 0){
            n -= n % a;
            n += a;
        }   

        int each_row = n / a;
        double* partA = new double[each_row * each_row + 2]; // A[beginRow:endRow, :]
        double* partB = new double[each_row * each_row + 2]; // B[:, beginColumn:endColumn]
        double* partC = new double[each_row * each_row + 2]; // C[beginRow:endRow, beginColumn:endColumn]
        int beginRow, beginColumn;
            
        // Data generation
        if (my_rank == 0){  
            // Prepare data
            cout << "n = " << n << endl;
            A = new double[n * n + 2];
            B = new double[n * n + 2];
            C = new double[n * n + 2];
            matGene(A, saveN, n);
            matGene(B, saveN, n); 
            for(int ii = 0; ii < n; ii++){
                for(int jj = 0; jj < n; jj++){
                    C[ii * n + jj] = 0;
                }
            }  
            beginTime = MPI_Wtime();   
        }

        for(int k = 0; k < a; k++){
            int begin_part = k * each_row;
            // k th
            // Data Distributing
            if (my_rank == 0){
                for(int i = 0; i < comm_sz; i++){
                    // A[beginRow:beginRow+each_row, begin_part:begin_part+each_row]
                    // B[begin_part:begin_part+each_row, beginColumn:beginColumn+each_row]
                    beginRow = (i / a) * each_row;
                    beginColumn = (i % a) * each_row;
                    if(i == 0){
                        // Copy Straightly
                        for(int ii = 0; ii < each_row; ii++){
                            for(int jj = 0; jj < each_row; jj++){
                                partA[ii * each_row + jj] = A[(beginRow + ii) * n + (begin_part + jj)];
                                partB[ii * each_row + jj] = A[(begin_part + ii) * n + (beginColumn + jj)];
                            }
                        }
                    }
                    else{
                        for(int ii = 0; ii < each_row; ii++){
                            MPI_Send(&A[(beginRow + ii) * n + begin_part], each_row, MPI_DOUBLE, i, i * each_row + ii, MPI_COMM_WORLD);
                            MPI_Send(&B[(begin_part + ii) * n + beginColumn], each_row, MPI_DOUBLE, i, (i + comm_sz) * each_row + ii, MPI_COMM_WORLD);
                        }
                    }
                }
            }
            // Data Receive
            if (my_rank != 0){
                for(int ii = 0; ii < each_row; ii++){
                    MPI_Recv(&partA[ii * each_row + 0], each_row, MPI_DOUBLE, 0, my_rank * each_row + ii, MPI_COMM_WORLD, &status);
                    MPI_Recv(&partB[ii * each_row + 0], each_row, MPI_DOUBLE, 0, (my_rank + comm_sz) * each_row + ii, MPI_COMM_WORLD, &status);
                }

            }
            // Calculation
            matMulti(partA, partB, partC, each_row, each_row, each_row);
            // Return Result
            if (my_rank != 0){
                for(int ii = 0; ii < each_row; ii++){
		    MPI_Send(&partC[ii * each_row + 0], each_row, MPI_DOUBLE, 0, (my_rank + 2 * comm_sz) * each_row + ii, MPI_COMM_WORLD);
                }
            }
            // Data Collection & add
            if (my_rank == 0){
                // C[beginRow:beginRow+each_row, beginColumn:beginColumn+each_row]
                for(int i = 0; i < comm_sz; i++){
                    beginRow = (i / a) * each_row;
                    beginColumn = (i % a) * each_row;
                    if(i == 0){
                        // Copy Straightly
                        for(int ii = 0; ii < each_row; ii++){
                            for(int jj = 0; jj < each_row; jj++){
                                C[(beginRow + ii) * n + (beginColumn + jj)] += partC[ii * each_row + jj];
                            }
                        }
                    }
                    else{  
                        for(int ii = 0; ii < each_row; ii++){
                            double* tmp_partC = new double[each_row + 2];
		            MPI_Recv(&tmp_partC[0], each_row, MPI_DOUBLE, i, (i + 2 * comm_sz) * each_row + ii, MPI_COMM_WORLD, &status);
                            for(int jj = 0; jj < each_row; jj++){
                                C[(beginRow + ii) * n + (beginColumn + jj)] += tmp_partC[jj];
                            }
                            delete[] tmp_partC;
                        }
                    }
                }
            }
        }

        if (my_rank == 0){
            endTime = MPI_Wtime();
            // Output   
	    cout << "Time: " << endTime - beginTime << endl; 
	    print_matrix("A", A, saveN, saveN, n);
            print_matrix("B", B, saveN, saveN, n);
            print_matrix("C", C, saveN, saveN, n);
            
            delete[] A;
            delete[] B;
            delete[] C;
        }
        delete[] partA;
        delete[] partB;
        delete[] partC;
    }           
    MPI_Finalize();
    return 0;
}
