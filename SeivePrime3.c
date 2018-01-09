#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define MIN(a,b) ((a)<(b)?(a):(b))
#define BLOCK_LOW(id,p,n) ((id)*((n-1)/2)/p)*2
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-2)
#define BLOCK_SIZE(id,p,n) (BLOCK_LOW((id)+1,p,n)-BLOCK_LOW((id),p,n))/2

int main(int argc, char *argv[]){

    long long int i,n,low_value,high_value,size,index,prime;
    double  elapsed_time;
    int count,global_count,id,p,step=640000;
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
    size=BLOCK_SIZE(id,p,n-1);
    char * marked = (char *) malloc (size);
       if (marked == NULL) {
          printf ("Cannot allocate enough memory\n");
          MPI_Finalize();
          exit (1);
        }
    char * SevingPrime = (char *) malloc ((int)sqrt((double) n));
    int * tem=(int *) malloc ((int)sqrt((double)n));
     for (i=0;i<(int)sqrt((double) n);i++) SevingPrime[i]=0;
    

    for (i = 0; i < size; i++) marked[i] = 0;
    index = 0;
    prime = 3;
    long long int prm=1;
    tem[0]=3;
    long long int first0;  
    do {
     first0=(prime*prime-3)/2;
     for(i=first0;i<=(int)sqrt((double) n)/2;i+=prime) SevingPrime[i]=1;
     while(SevingPrime[++index]);
     prime=2*index+3;
     tem[prm++]=prime;

    }while(prime*prime<=n);
    
    for(i=0;i<(size-1)/step+1;i++){
        index=0;
        prime=3;
       do {
         if (prime * prime > low_value) 
              first0 = (prime * prime - low_value)/2;
         else {
             if (!(low_value % prime)) first0 = 0;
             else {
                if((prime - (low_value % prime)+low_value)%2==0) first0 = (2*prime - (low_value % prime))/2;
                else first0=(prime-(low_value%prime))/2;} 
         }
         
         for (prm = first0+step*i; prm < size && prm<step*i+step; prm += prime) marked[prm] = 1;
         index++;
         prime=tem[index];
      } while (prime * prime <= low_value+2*step && prime*prime<=n);
      low_value+=2*step;
   }


   count = 0;
   
   for (i = 0; i < size; i++)
      if (!marked[i]) count++;
   MPI_Reduce (&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
   elapsed_time += MPI_Wtime();
   if (!id) {
      if(p==32) {
          printf("Part3\n");
          printf("The total number of prime within %lld: %d, total time: %10.6f, total node 1\n",n,global_count+1,elapsed_time);}
      if(p==64) printf("The total number of prime within %lld: %d, total time: %10.6f, total node 2\n",n,global_count+1,elapsed_time);      
      if(p==128) printf("The total number of prime within %lld: %d, total time: %10.6f, total node 4\n",n,global_count+1,elapsed_time);
   
      if(p==256) printf("The total number of prime within %lld: %d, total time: %10.6f, total node 8\n",n,global_count+1,elapsed_time);
   }
   MPI_Finalize ();
   return 0;
}



