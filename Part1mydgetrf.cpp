#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include "lapacke.h"
#include "blas.h"
using namespace std;
int  main(){
        void mydgetrf(double *, int *, int);
        void mydtrsm(int *pvt, double *	B, double *A, int n, char UPLO);
        int n=1000,i=0,j=0;
        printf("The matrix size is: %d\n",n);
	      double *A=new double[n*n];
        double *AA=new double[n*n];
	      double *B=new double[n];
        double *C=new double[n];
        int *pvt=new int[n];
        double *x=new double[n];
	      srand(217);
        
        char    TRANS = 'N';
        int     INFO = n;
        int     LDA = n;
        int     LDB = n;
        int     N = n;
        int     NRHS = 1;
        int     IPIV[n];
     //initialize the input matrix A(column dominant)/AA(row dominant) and vector B/C
        for (i = 0; i < n; i++){
           for(j=0;j<n;j++){
	            A[i*n+j] = (double)rand()/(double)RAND_MAX*99.0+1.0;
	            AA[j*n+i]=A[i*n+j];
           }
      	}
	      for (i=0;i<n;i++){
	         B[i] = (double)rand()/(double)RAND_MAX*99.0+1.0;
           C[i]=B[i];
           pvt[i]=i;
	      }
        clock_t start0=clock(); 
        //LU fractorization by calling function dgetrf() in LAPACK
        LAPACK_dgetrf(&N,&N,A,&LDA,IPIV,&INFO);
        clock_t end0=clock();
        printf("Runing time of dgetrf is: %f\n", (double)(end0-start0)/CLOCKS_PER_SEC);
        char     SIDE = 'L';
        char     UPLO = 'L';
        char     DIAG = 'U';
        int      M    = 1;
        double   a    = 1.0;
       
        for(i = 0; i < N; i++)
        {
           double tmp = B[IPIV[i]-1];
           B[IPIV[i]-1] = B[i];
           B[i] = tmp;
        }


    // forward  L(Ux) = B => y = Ux
        dtrsm_(&SIDE,&UPLO,&TRANS,&DIAG,&N,&M,&a,A, &N, B, &N);

        UPLO = 'U';
        DIAG = 'N';
    // backward Ux = y
       
        dtrsm_(&SIDE,&UPLO,&TRANS,&DIAG,&N,&M,&a,A, &N, B, &N);
       
        clock_t start1=clock();
    //call the self implemented LU fractorization function mydgetrf()
        mydgetrf(AA,pvt,n);
        clock_t end1=clock();
        printf("Runing time of mydgetrf is: %f\n", (double)(end1-start1)/CLOCKS_PER_SEC);
    // forward  L(Ux) = B => y = Ux
        mydtrsm(pvt,C,AA,n,'L');
    // backward Ux = y
        mydtrsm(pvt,C,AA,n,'U');

       for(i=0;i<n;i++){

          if(fabs(C[i]-B[i])>0.01) {
              printf("The result is wrong!\n");
              return 0;
            }
       }

       printf("The result is correct!\n");

}

//self implemented LU fractorization of the coefficient matrix A
void mydgetrf(double *A,int *pvt, int n){
	int i=0;
	int t=0,temps=0,k=0,j=0;
	double tempv=0.0;

	for(i=0;i<n;i++){
		int maxind=i;
		double max=fabs(A[i*n+i]);
		for(t=i+1;t<n;t++){
			if(fabs(A[t*n+i])>max){
				maxind=t;
			    max=fabs(A[t*n+i]);}
		}
		if(max==0) printf("LUfactoration failed: coefficient matrix is singlular");
		else if(maxind!=i){
				temps=pvt[i];
				pvt[i]=pvt[maxind];
				pvt[maxind]=temps;
				for(k=0;k<n;k++){
					tempv=A[i*n+k];
				    A[i*n+k]=A[maxind*n+k];
				    A[maxind*n+k]=tempv;
				}
			}
                
		for(j=i+1;j<n;j++){
			A[j*n+i]=A[j*n+i]/A[i*n+i];
			for(k=i+1;k<n;k++){
				A[j*n+k]=A[j*n+k]-A[j*n+i]*A[i*n+k];
			}
		}
	}
}

void mydtrsm(int *pvt, double *B, double *A, int n, char UPLO){
    double *x=new double[n];
    int i=0,j=0;
    double sum=0.0;
    if(UPLO=='L'){
	    x[0]=B[pvt[0]];
	     for(i=1;i<n;i++){
          sum=0.0;
		    for(j=0;j<i;j++){
		       sum+=x[j]*A[i*n+j];
	      }
           x[i]=B[pvt[i]]-sum;
	     }
    } else{

	      x[n-1]=B[n-1]/A[(n-1)*n+n-1];

	      for(i=n-2;i>=0;i--){
                sum=0.0;
		      for(j=i+1;j<n;j++){
		        sum+=x[j]*A[i*n+j];
		      }
		      x[i]=(B[i]-sum)/A[i*n+i];
	      }
    }

    for(i=0;i<n;i++){
       B[i]=x[i];
    }
}
