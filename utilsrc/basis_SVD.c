#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "lapacke.h"

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, lapack_int m, lapack_int n, float* a, lapack_int lda );
extern void print_int_vector( char* desc, lapack_int n, lapack_int* a );

#define min(A,B) ((A<B)? A:B)


static float u[16260004];  // u[M][M] NEVER MAKE THIS LOCAL!! WILL OVERFLOW THE STACK!!!!! ALWAYS USE LARGE ARRAYS AS GLOBAL VARIABLES

int main(int argc, char *argv[]) {
        FILE *ftraj = fopen("trajectory_matrix.txt", "r");
	FILE *fleft = fopen("u_left.txt","w");
	FILE *fright = fopen("v_right.txt","w");
	FILE *fsing = fopen("singular.txt","w");
	
	int i, j, info;
	int M = atoi(argv[2]);
	int N = 3*atoi(argv[3]);

	float s[N];
	float vt[N*N];
	float superb[min(M,N) - 1];
	float a[M*N];

	for (i = 0; i < M; i++)
		for (j = 0; j < N; j++)
			fscanf(ftraj, "%f", &a[i*N+j]);
	
        printf( "LAPACKE_sgesvd SVD with single precision floating point numbers\n" );
        info = LAPACKE_sgesvd( LAPACK_ROW_MAJOR, 'A', 'A', M, N, a, N, s, u, M, vt, N, superb);
                    
        if( info > 0 ) {
                printf( "Failed to converge\n" );
                exit( 1 );
        }
	
	for (i = 0; i < M; i++)
		for (j = 0; j < N; j++)
			fprintf(fleft, "%f ", u[j*M + i]);
	
	for(i = 0; i < N; i++)
		for( j = 0; j < N; j++)
			fprintf(fright, "%f ", vt[i*N + j]);

	for(i = 0; i < N; i++)
		fprintf(fsing, "%f ", s[i]);
	
	fclose(ftraj);
	fclose(fleft);
	fclose(fright);
	fclose(fsing);
        exit( 0 );
}

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, lapack_int m, lapack_int n, float* a, lapack_int lda ) {
        lapack_int i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %6.2f", a[i*lda+j] );
                printf( "\n" );
        }
}

/* Auxiliary routine: printing a vector of integers */
void print_int_vector( char* desc, lapack_int n, lapack_int* a ) {
        lapack_int j;
        printf( "\n %s\n", desc );
        for( j = 0; j < n; j++ ) printf( " %6i", a[j] );
        printf( "\n" );
}
