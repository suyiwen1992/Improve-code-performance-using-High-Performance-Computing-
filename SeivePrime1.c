#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//#include "MyMPI.h"
#define MIN(a,b) ((a)<(b)?(a):(b))
#define BLOCK_LOW(id,p,n) ((id)*((n-1)/2)/p)*2
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-2)
#define BLOCK_SIZE(id,p,n) (BLOCK_LOW((id)+1,p,n)-BLOCK_LOW((id),p,n))/2
//#define BLOCK_OWNER(index,p,n) ((((p)*(index)+1)-1)/(n))

int main(int argc, char *argv[]){

    long long int n,low_value,high_value,size,index,prime,first;
    double  elapsed_time;
    int id,p,count,global_count;
    MPI_Init (&argc,&argv);
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time=-MPI_Wtime();
    MPI_Comm_rank(MPI_COMM_WORLD,&id);
    MPI_Comm_size (MPI_COMM_WORLD,&p);
    if(argc!=2){
         if(!id) printf("Command line: %s <m>\n",argv[0]);
         MPI_Finalize();
         exit(1);
    }

    n=atoi(argv[1]);
    n=10000000000;
    low_value=3+BLOCK_LOW(id,p,n-1);

    high_value=3+BLOCK_HIGH(id,p,n-1);
  //  printf("The low_value is %d, the high_value is %d\n",low_value,high_value);
    size=BLOCK_SIZE(id,p,n-1);
    int proc0_size=((n-1)/2-1)/p;
    if((3+proc0_size)<(int)sqrt((double) n)){
         if(!id) printf("Too many processes\n");
         MPI_Finalize();
         exit(1);
     }
    char * marked = (char *) malloc (size);
       if (marked == NULL) {
          printf ("Cannot allocate enough memory\n");
          MPI_Finalize();
          exit (1);
        }

    for (long long int i = 0; i < size; i++) marked[i] = 0;
    if (!id) index = 0;
    prime = 3;
    do {
      if (prime * prime > low_value)
         first = (prime * prime - low_value)/2;
      else {
         if (!(low_value % prime)) first = 0;
         else {
             if((prime - (low_value % prime)+low_value)%2==0) first = (2*prime - (low_value % prime))/2;
             else first=(prime-(low_value%prime))/2;} 
      }
      for (long long int i = first; i < size; i += prime) marked[i] = 1;

      if (!id) {
         while (marked[++index]);
         prime = 2*index + 3;
      }
      MPI_Bcast (&prime,  1, MPI_INT, 0, MPI_COMM_WORLD);
   } while (prime * prime <= n);
   count = 0;
   for (long long int i = 0; i < size; i++)
      if (!marked[i]) count++;
   MPI_Reduce (&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
   elapsed_time += MPI_Wtime();
   if (!id) {
      if(p==32) {
         printf("Part1\n");
         printf("The total number of prime within %lld: %d, total time: %10.6f, total node 1\n",n,global_count+1,elapsed_time);}
      if(p==64) printf("The total number of prime within %lld: %d, total time: %10.6f, total node 2\n",n,global_count+1,elapsed_time);
      if(p==128) printf("The total number of prime within %lld: %d, total time: %10.6f, total node 4\n",n,global_count+1,elapsed_time);
   
      if(p==256) printf("The total number of prime within %lld: %d, total time: %10.6f, total node 8\n",n,global_count+1,elapsed_time);
   }
   MPI_Finalize ();
   return 0;
}



