#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include "lapacke.h"
#include "blas.h"
using namespace std;
void optimize(int i, int l, int k, int b, int n, double *A[]);
int  main(){
        void mydgetrfBlock(double *A[], int *pvt, int n, int b);
        void mydtrsm(int *pvt, double *B, double *A1[], int n, char UPLO);
       // void optimize(int i, int l, int k, int b, int n, double *A[]);

        int n=1000,i=0,j=0,b=20;

        printf("The matrix size is: %d\n",n);
	double *A=new double[n*n];
        double *A1[n];
	double *B=new double[n];
        double *C=new double[n];
        int *pvt=new int[n];
	srand(217);
        
        char    TRANS = 'N';
        int     INFO = n;
        int     LDA = n;
        int     LDB = n;
        int     N = n;
        int     NRHS = 1;
        int     IPIV[n];
       
        for(i=0;i<n;i++){
            A1[i]=(double *) malloc(sizeof(double)*n);
        }
        for (i = 0; i < n; i++){
           for(j=0;j<n;j++){
	    A[i*n+j] = (double)rand()/(double)RAND_MAX*99.0+1.0;
	    A1[j][i]=A[i*n+j];
           }
      	}
	for (i=0;i<n;i++){
	   B[i] = (double)rand()/(double)RAND_MAX*99.0+1.0;
           C[i]=B[i];
            pvt[i]=i;
	}
        clock_t start0=clock();
        LAPACK_dgetrf(&N,&N,A,&LDA,IPIV,&INFO);
        clock_t end0=clock();
        printf("Runnin time of LAPACK_dgetrf is: %f\n", (double)(end0-start0)/CLOCKS_PER_SEC);

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
        mydgetrfBlock(A1,pvt,n,b);
        clock_t end1=clock();
        printf("Running time of mydgetrfBlock is: %f\n", (double)(end1-start1)/CLOCKS_PER_SEC);
        mydtrsm(pvt,C,A1,n,'L');

        mydtrsm(pvt,C,A1,n,'U');

       for(i=0;i<n;i++){

          if(fabs(C[i]-B[i])>0.01) {
              printf("The result is wrong!");
              return 0;
            }
       }

       printf("The result is correct!\n");

}
void mydgetrfBlock(double *A[],int *pvt, int n, int b){
	int i=0;
	int t=0,temps=0,k=0,j=0,l=0;
	double *tempv;
       
	for(i=0;i<n;i+=b){
           for(j=i;j<i+b&&j<n;j++){
		int maxind=j;
		double max=fabs(A[j][j]);
		for(t=j+1;t<n;t++){
			if(fabs(A[t][j])>max){
				maxind=t;
			    max=fabs(A[t][j]);}
		}
		if(max==0) printf("LUfactoration failed: coefficient matrix is singlular");
		else if(maxind!=j){
				temps=pvt[j];
				pvt[j]=pvt[maxind];
				pvt[maxind]=temps;

                                tempv=A[j];
                                A[j]=A[maxind];
                                A[maxind]=tempv;
				}
                
		for(l=j+1;l<n;l++){
			A[l][j]=A[l][j]/A[j][j];
                        register double R=A[l][j];
			for(k=j+1;k<i+b&&k<n;k++){
				A[l][k]-=R*A[j][k];
			}
		}
             }
          

                for(l=i;l<i+b&&l<n;l++){
                    for(k=i+b;k<n;k++){
                       double sum=0.0;
                       for(j=i;j<l;j++){
                          sum+=A[l][j]*A[j][k];
                       }
                       A[l][k]-=sum;
                    }

                }

 
               for(l=i+b;l<n;l+=2){
                   for(k=i+b;k<n;k+=2){
                      register double A00=A[l][k]; register double A01=A[l][k+1];
                      register double A10=A[l+1][k]; register double A11=A[l+1][k+1];
                      for(j=i;j<i+b&&j<n;j+=2){
                         register double B0=A[l][j]; register double B1=A[l+1][j];
                         register double C0=A[j][k]; register double C1=A[j][k+1];
                         A00-=B0*C0;
                         A01-=B0*C1;
                         A10-=B1*C0;
                         A11-=B1*C1;
                         B0=A[l][j+1]; B1=A[l+1][j+1];
                         C0=A[j+1][k]; C1=A[j+1][k+1];
                         A00-=B0*C0;
                         A01-=B0*C1;
                         A10-=B1*C0;
                         A11-=B1*C1;
                      }

                      A[l][k]=A00;
                      A[l][k+1]=A01;
                      A[l+1][k]=A10;
                      A[l+1][k+1]=A11;
                     // double sum=0.0;
                     
                     //sum+=A[l][j]*A[j][k];
                     // A[l][k]=A[l][k]-sum;                    
                  }
               }
/*
              for(l=i+b;l<n;l+=b){
                  for(k=i+b;k<n;k+=b){

                     optimize(i,l,k,b,n,A);
                   }

               }*/
/*
               for(l=i+b;l<n;l++){
                  for(k=i+b;k<n;k++){
                     double sum=0.0;
                     for(j=i;j<i+b;j++){

                        sum+=A[l][j]*A[j][k];
                     }
                     A[l][k]-=sum;
                  }
               }*/

	}
}
/*
void optimize(int i, int l, int k, int b, int n, double *A[]){
  int q,e,r,p;
  
      for(q=l;l<l+b;l+=2){
         for(e=k;e<k+b;k+=2){
            register double A00=A[q][e]; register double A01=A[q][e+1];
            register double A10=A[q+1][e];register double A11=A[q+1][e+1];
            for(r=i;r<i+b;i+=2){
               for(p=0;p<2;p++){
                  register double B0=A[q][r+p];
                  register double B1=A[q+1][r+p];
                  register double C0=A[r+p][e];
                  register double C1=A[r+p][e+1];
                  
                  A00=A00-B0*C0;
                  A01=A01-B0*C1;
                  A10=A10-B1*C0;
                  A11=A11-B1*C1;
                }
            }
                
      
         A[q][e]=A00;
         A[q][e+1]=A01;
         A[q+1][e]=A10;
         A[q+1][e+1]=A11;
         }
      }

}*/

void mydtrsm(int *pvt, double *B, double *A[], int n, char UPLO){
    double *x=new double[n];
    int i=0,j=0;
   double sum=0.0;
     if(UPLO=='L'){
	x[0]=B[pvt[0]];
	for(i=1;i<n;i++){
               sum=0.0;
		for(j=0;j<i;j++){
		    sum+=x[j]*A[i][j];
	         }
           x[i]=B[pvt[i]]-sum;
	}
     } else{

	x[n-1]=B[n-1]/A[n-1][n-1];

	for(i=n-2;i>=0;i--){
                sum=0.0;
		for(j=i+1;j<n;j++){
		   sum+=x[j]*A[i][j];
		}
		x[i]=(B[i]-sum)/A[i][i];
	}
      }

    for(i=0;i<n;i++){
       B[i]=x[i];
    }
}
