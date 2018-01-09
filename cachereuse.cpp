#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <iostream>

int main(){
int n=1000;
int i=0,j=0,k=0;
double *a=new double[n*n];
double *b=new double[n*n];
double *c0=new double[n*n];
double *c1=new double[n*n];
double *c2=new double[n*n];
double *c3=new double[n*n];
double *c4=new double[n*n];
double *c5=new double[n*n];
double *c6=new double[n*n];
double *c7=new double[n*n];
double *c8=new double[n*n];
double *c9=new double[n*n];
double *c10=new double[n*n];
double *c11=new double[n*n];
double *c12=new double[n*n];


srand(217);

for (i = 0; i < n*n; i++){
	a[i] = rand()%RAND_MAX;
	b[i] = rand()%RAND_MAX;
	c0[i] = 0.0;
	c1[i] = 0.0;
	c2[i] = 0.0;
	c3[i] = 0.0;
	c4[i] = 0.0;
	c5[i] = 0.0;
	c6[i] = 0.0;
	c7[i] = 0.0;
	c8[i] = 0.0;
	c9[i] = 0.0;
	c10[i] = 0.0;
	c11[i] = 0.0;
	c12[i] = 0.0;



}

clock_t start0=clock();
for(i=0;i<n;i++){
	for(j=0;j<n;j++){
		register double r0=c0[i*n+j];
		for(k=0;k<n;k++){
			r0+=a[i*n+k]*b[k*n+j];}
		c0[i*n+j]=r0;
	}
}
clock_t end0=clock();
printf("Runing time of simple triple loop (ijk version) algorithm is: %f", (double)(end0-start0)/CLOCKS_PER_SEC);

clock_t start1=clock();
for(i=0;i<n;i++){
	for(k=0;k<n;k++){
		register double r1=a[i*n+k];
		for(j=0;j<n;j++){
			c1[i*n+j]+=r1*b[k*n+j];}
	}
}
clock_t end1=clock();
printf("Runing time of simple triple loop (ikj version) algorithm is: %f", (double)(end1-start1)/CLOCKS_PER_SEC);

clock_t start2=clock();
for(j=0;j<n;j++){
	for(i=0;i<n;i++){
		register double r2=c2[i*n+j];
		for(k=0;k<n;k++){
			r2+=a[i*n+k]*b[k*n+j];}
		c2[i*n+j]=r2;
	}
}
clock_t end2=clock();
printf("Runing time of simple triple loop (jik version) algorithm is: %f", (double)(end2-start2)/CLOCKS_PER_SEC);

clock_t start3=clock();
for(j=0;j<n;j++){
	for(k=0;k<n;k++){
		register double r3=b[k*n+j];
		for(i=0;i<n;i++){
		c3[i*n+j]+=a[i*n+k]*r3;}
	}
}
clock_t end3=clock();
printf("Runing time of simple triple loop (jki version) algorithm is: %f", (double)(end3-start3)/CLOCKS_PER_SEC);

clock_t start4=clock();
for(k=0;k<n;k++){
	for(i=0;i<n;i++){
		register double r4=a[i*n+k];
		for(j=0;j<n;j++){
		c4[i*n+j]+=r4*b[k*n+j];}
	}
}
clock_t end4=clock();
printf("Runing time of simple triple loop (kij version) algorithm is: %f", (double)(end4-start4)/CLOCKS_PER_SEC);

clock_t start5=clock();
for(k=0;k<n;k++){
	for(j=0;j<n;j++){
		register double r5=b[k*n+j];
		for(i=0;i<n;i++){
		c5[i*n+j]+=a[i*n+k]*r5;}
	}
}
clock_t end5=clock();
printf("Runing time of simple triple loop (kji version) algorithm is: %f", (double)(end5-start5)/CLOCKS_PER_SEC);


int B = 10;
int i1, j1, k1;

clock_t start6 = clock();
for (i = 0; i < n; i+=B)
	for (j = 0; j < n; j+=B)
		for (k = 0; k < n; k+=B)
			for (i1 = i; i1 < i+B; i1++)
				for (j1 = j; j1<j+B; j1++){
					register double r6 = c1[i1*n+j1];
					for (k1 = k; k1 < k+B; k1++)
						r6 += a[i1*n+k1]*b[k1*n+j1];
					c6[i1*n+j1] = r6;
				}
clock_t end6=clock();
printf("Runing time of blocked version (ijk) algorithm is: %f", (double)(end6-start6)/CLOCKS_PER_SEC);
std::cout << "here" << std::endl;

clock_t start7 = clock();
for (i = 0; i < n; i+=B)
	for (k = 0; k < n; k+=B)
		for (j = 0; j < n; j+=B)
			for (i1 = i; i1 < i+B; i1++)
				for (k1 = k; k1<k+B; k1++){
					register double r7 = a[i1*n+k1];
					for (j1 = j; j1 < j+B; j1++)
						c7[i1*n+j1] += r7*b[k1*n+j1];
				}
clock_t end7=clock();
printf("Runing time of blocked version (ikj) algorithm is: %f", (double)(end7-start7)/CLOCKS_PER_SEC);
std::cout << "here" << std::endl;

clock_t start8 = clock();
for (j = 0; j < n; j+=B)
	for (i = 0; i < n; i+=B)
		for (k = 0; k < n; k+=B)
			for (j1 = j; j1 < j+B; j1++)
				for (i1 = i; i1<i+B; i1++){
					register double r8 = c1[i1*n+j1];
					for (k1 = k; k1 < k+B; k1++)
						r8 += a[i1*n+k1]*b[k1*n+j1];
					c8[i1*n+j1] = r8;
				}
clock_t end8=clock();
printf("Runing time of blocked version (jik) algorithm is: %f", (double)(end8-start8)/CLOCKS_PER_SEC);
std::cout << "here" << std::endl;

clock_t start9 = clock();
for (j = 0; j < n; j+=B)
	for (k = 0; k < n; k+=B)
		for (i = 0; i < n; i+=B)
			for (j1 = j; j1 < j+B; j1++)
				for (k1 = k; k1<k+B; k1++){
					register double r9 = b[k1*n+j1];
					for (i1 = i; i1 < i+B; i1++)
						c9[i1*n+j1] += a[i1*n+k1]*r9;
				}
clock_t end9=clock();
printf("Runing time of blocked version (jki) algorithm is: %f", (double)(end9-start9)/CLOCKS_PER_SEC);
std::cout << "here" << std::endl;

clock_t start10 = clock();
for (k = 0; k < n; k+=B)
	for (i = 0; i < n; i+=B)
		for (j = 0; j < n; j+=B)
			for (k1 = k; k1 < k+B; k1++)
				for (i1 = i; i1<i+B; i1++){
					register double r10 = a[i1*n+k1];
					for (j1 = j; j1 < j+B; j1++)
						c10[i1*n+j1] += r10*b[k1*n+j1];
				}
clock_t end10=clock();
printf("Runing time of blocked version (kij) algorithm is: %f", (double)(end10-start10)/CLOCKS_PER_SEC);
std::cout << "here" << std::endl;

clock_t start11 = clock();
for (k = 0; k < n; k+=B)
	for (j = 0; j < n; j+=B)
		for (i = 0; i < n; i+=B)
			for (k1 = k; k1 < k+B; k1++)
				for (j1 = j; j1<j+B; j1++){
					register double r11 = b[k1*n+j1];
					for (i1 = i; i1 < i+B; i1++)
						c11[i1*n+j1] += a[i1*n+k1]*r11;
				}
clock_t end11=clock();
printf("Runing time of blocked version (kji) algorithm is: %f", (double)(end11-start11)/CLOCKS_PER_SEC);
std::cout << "here" << std::endl;
clock_t start12 = clock();
for (i = 0; i < n; i+=B)
	for (j = 0; j < n; j+=B)
		for (k = 0; k < n; k+=B)
			for (i1 = i; i1 < i+B; i1+=2)
				for (j1 = j; j1<j+B; j1+=2){
					register int t=i1*B+j1; register int tt=t+B;
		            register double c00=c12[t]; register double c01=c12[t+1]; register double c10=c12[tt]; register double c11=c12[tt+1];
					for (k1 = k; k1 < k+B; k1+=2){
						register int ta=i1*B+k1; register int tta=ta+B; register int tb=k1*B+j1; register int ttb=tb+B;
			            register double a00=a[ta]; register double a01=a[ta+1]; register double a10=a[tta]; register double a11=a[tta+1];
			            register double b00=b[tb]; register double b01=b[tb+1]; register double b10=b[ttb]; register double b11=b[ttb+1];

			            c00+=a00*b00+a01*b10;
			            c01+=a00*b01+a01*b11;
			            c10+=a10*b00+a11*b10;
			            c11+=a10*b01+a11*b11;}
			        c12[t]=c00;
		            c12[t+1]=c01;
		            c12[tt]=c10;
		            c12[tt+1]=c11;
				}
clock_t end12=clock();
printf("Runing time of blocked version (kij) algorithm is: %f", (double)(end12-start12)/CLOCKS_PER_SEC);
std::cout << "here" << std::endl;
}