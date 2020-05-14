#include <mkl.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

// To Run: 
// gcc -o matrix matrix.c  -DMKL_ILP64 -m64 -Iinclude -Wl,--start-group libmkl_intel_ilp64.a libmkl_gnu_thread.a libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl
void print_matrix(double *array, int height, int width){
	int ii, kk;
    printf("%d %d\n",height, width);
 
    for(ii=0;ii<height;ii++)
    {
        for(kk=0;kk<width;kk++)
        {
            printf("%.2lf ", array[ii*width+kk]);
        }
        printf("\n");
    }
 
}

void write_file(char* filename, double *array, int m, int n){
    FILE* f = fopen(filename, "w");
    int ii=0,jj;

    for(jj=0; jj<m; jj++){

        for(ii=0; ii<n; ii++)
		fprintf(f, "%.2lf ", array[jj*m+ii]);
 	fprintf(f, "\n");	
    }
    fclose(f);
}

double *read_file(char* filename, int height, int width){
    FILE* f = fopen(filename, "r");
    int ii=0,jj;
    double *array = (double *)mkl_malloc( height*width*sizeof( double ), 64 );

    for(jj=0; jj<height; jj++)
        for(ii=0; ii<width; ii++)
            if(fscanf(f, "%lf ", &array[jj*width + ii]) != 1)
		exit(1);
    print_matrix(array, height, width);    
    fclose(f);
    return array;
}
void matrix_multiply(char *filename, char *filename2,char *filename3, int m, int k, int n)
{
	double *A = read_file(filename, m, k);
	double *B = read_file(filename2, k, n);
    	double *C = (double *)mkl_malloc( m*n*sizeof( double ), 64 );

   	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, A, k, B, n, 0.0, C, n);
	print_matrix(C, m, n); 
	write_file(filename3, C, m, n);
	mkl_free(A);
	mkl_free(B);
	mkl_free(C);	

}

int main(int argc, char* argv[]){
	char *filename = argv[1];
	char *filename2 = argv[2];
	int m = atoi(argv[3]);
	int k = atoi(argv[4]);
	int n = atoi(argv[5]);
	char *filename3 = argv[6];
	matrix_multiply(filename, filename2,filename3, m, k,n);
	return 0;
}
