#include <iostream>
#include <math.h>
#include <algorithm>
#include <stdio.h>
#include<float.h>


#define THREADS_PER_BLOCK 1024           //max of the threads in one block is 1024

// Kernel function to add the elements of two arrays
__global__ 
void iteration(double *d_A,double *d_B,int n)
{
    int i=blockIdx.x*blockDim.x+threadIdx.x; 
    
    if(i<n*n){
      if((0<i&&i<n)||i%n==0||(i+1)%n==0||(n*n-n<i&&i<n*n-1)){
        d_B[i]=d_A[i];
      }
      else{
        double local[4];
        double temp;
        double first_small;
        double secnd_small;


        local[0]=d_A[i+n-1];
        local[1]=d_A[i+n+1];
        local[2]=d_A[i-n-1];
        local[3]=d_A[i-n+1];

        if (local[0]>local[1]){
          first_small = local[1];
          secnd_small = local[0];
        }
        else{
          first_small = local[0];
          secnd_small = local[1];
        }

        if(local[2]<first_small){
          secnd_small = first_small;
          first_small = local[2];
        }
        else if(local[2]<secnd_small){
          secnd_small = local[2];
        }

        if(local[3]<first_small){
          secnd_small = first_small;
          first_small = local[3];
        }
        else if(local[3]<secnd_small){
          secnd_small = local[3];
        }
        d_B[i]=secnd_small+d_A[i];
      }
    }    
}

__global__ 
void sumblock(double *d_A, int size, double*sum_temp)
{

  extern __shared__ double sum_block[];

  int i=blockIdx.x*blockDim.x+threadIdx.x;
  int k;
  int bound;

  double temp=0;
  if(i<size){
    temp=d_A[i]+temp;
  }

  sum_block[threadIdx.x]=temp;
  __syncthreads();

  k = blockDim.x;
  while(k>1){
    if(k%2==0){
      bound=k/2;
      if (threadIdx.x<bound){
        sum_block[threadIdx.x]=sum_block[threadIdx.x]+sum_block[threadIdx.x+bound];
      }
      __syncthreads();
      k=k/2;
    }
    else{
      bound=k/2;
      if (threadIdx.x<=bound||threadIdx.x!=0){
        sum_block[threadIdx.x]=sum_block[threadIdx.x]+sum_block[threadIdx.x+bound];
      }
      __syncthreads();
      k=k/2+1;
    }             
  } 

  if(threadIdx.x==0){
    sum_temp[blockIdx.x] =sum_block[0];
  }
}

int main(int argc, char **argv)
{
  int n = atoi(argv[1]);          //The size of the matrix
  
  double *A;                      //The definition of matrix A
  double *d_A;                    //The definition of matrix A in gpu
  double *d_B;                    //The definition of iterated matrix A in gpu

  cudaEvent_t start;              // the start time of gpu calculation
  cudaEvent_t end;                // the end time of gpu calculation
  float elapsedTime=0;            // the elapsed time of gpu calculation

  double sum=0;                   //The sum of the matrix A after 10 iterations
  double *sum_temp;               //The sum of the matrix A after 10 iterations for each block on gpu
  double center=0;                //The center of the matrix A after 10 iterations
  double verification=0;          //The A(17,31) of the matrix A after 10 iterations

  int count=0;                    //To count the number of the iterations

  int grid;                       //The dimension of the grid on the gpu
  int block;                      //The dimension of the block on the gpu
  
  int size;                       //The size of elements in each iteration for sum
  double* sum_ptr;                //The pointer of the first element in each iteration for sum

  
  cudaEventCreate(&start);  
  cudaEventCreate(&end);

  // Allocate Unified Memory on the CPU
  A = (double*)malloc(n*n*sizeof(double));
  
  //initialize the matrix A on the host
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      A[i*n+j]=(1+cos(2*i)+sin(j))*(1+cos(2*i)+sin(j));
    }
  }

  // Allocate Memory on the GPU
  cudaMalloc(&d_A, n*n*sizeof(double)); 
  cudaMalloc(&d_B, n*n*sizeof(double));  

  //copy the data to the gpu from cpu
  cudaMemcpy(d_A, A, n*n*sizeof(double), cudaMemcpyHostToDevice);

  block=THREADS_PER_BLOCK;
  grid=(n*n%block==0)?n*n/block:(n*n/block+1);

  cudaMalloc(&sum_temp, grid*sizeof(double));

  cudaEventRecord(start);  
  
  while(count<10){
    // Run kernel on the GPU
    iteration<<<grid,block>>>(d_A,d_B,n);
    cudaDeviceSynchronize();
    double *temp;
    temp=d_B;
    d_B=d_A;
    d_A=temp;   
    count++;
    std::cout<<"The numbers of iterations is "<<count<<std::endl;
  }

  size=n*n;
  sum_ptr=d_A;

  while(grid!=0){
    sumblock<<<grid,block,block*sizeof(double)>>>(sum_ptr,size,sum_temp);
    cudaDeviceSynchronize();
    size=grid;
    if(grid==1){
      grid=0;
    }
    else{
      grid=(grid%block==0)?grid/block:grid/block+1;
    }    
    sum_ptr=sum_temp;
  }
  
  cudaEventRecord(end);
  
  
  cudaEventSynchronize(end);
  cudaEventElapsedTime(&elapsedTime, start, end);

  cudaMemcpy(&sum, sum_temp, sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(&center, (d_A+n/2*n+n/2), sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(&verification, (d_A+37*n+47), sizeof(double), cudaMemcpyDeviceToHost);


  std::cout <<"The elapsed_time of the cuda program is  "<< elapsedTime << std::endl;
  std::cout <<"The sum of the matrix A after 10 iterations is  "<< sum << std::endl;
  std::cout <<"The center of the matrix A after 10 iterations is  "<< center << std::endl;
  std::cout <<"The A(17,31) of the matrix A after 10 iterations is  "<< verification << std::endl;

  cudaFree(d_A);
  cudaFree(d_B);
  cudaFree(sum_temp);

  free(A);

  return 0;
}
