#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mpi.h"

void merge_sort(int n, double * a);
double * merge_array(int n, double * a, int m, double * b);
void swap(double * a, double * b);

int MPI_Sort_shell(int n, double * a, double max, int root, MPI_Comm comm, double * compTime, double * commTime);

int main(int argc, char * argv[]) {
    double commT, compT = 0;

    int rank, size;
    int n = 10000000, i;
    double m = 10.0;
    double * scattered_array, * array;

    // initialise MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    array = (double *) calloc(n, sizeof(double));

    if(rank == 0) {
        srand(((unsigned)time(NULL)+rank));
        for(i=0;i<n;i++) array[i] = ((double)rand()/RAND_MAX)*m;
    }

    MPI_Sort_shell(n, array, m, 0, MPI_COMM_WORLD, &compT, &commT);

    printf("Rank %d: Comm.T: %lfs, Comp.T: %lfs\n", rank, commT, compT);

    MPI_Finalize();
}

int MPI_Sort_shell(int n, double * array, double max, int root, MPI_Comm comm, double * compTime, double * commTime) {
    int rank, size;
    double commT, compT, time = 0;

    commT = (*commTime);
    compT = (*compTime);

    time = MPI_Wtime();
    
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    
    time = MPI_Wtime() - time;
    commT += time;
    time = MPI_Wtime();

    // scatter sub-array locally based on rank and size
    double * localArray = (double *)calloc(n/size, sizeof(double));
    
    time = MPI_Wtime() - time;
    compT += time;
    time = MPI_Wtime();
    
    MPI_Scatter(array, n/size, MPI_DOUBLE, localArray, n/size, MPI_DOUBLE, root, comm);

    time = MPI_Wtime() - time;
    commT += time;
    time = MPI_Wtime();

    // implement shell sort
    for (int gap=(n/size)/2; gap>0; gap/=2 ) {
        for(int i=gap;i<n;i++) {
            int temp = localArray[i];
            int j;

            for(j=i; j>=gap && localArray[j-gap] > temp; j -= gap) {
                localArray[j] = localArray[j-gap]; 
            }
            localArray[j] = temp;
        }
    }
    time = MPI_Wtime() - time;
    compT += time;
    time = MPI_Wtime();

    // gather sorted localArray into array
    MPI_Gather(localArray, n/size, MPI_DOUBLE, array, n/size, MPI_DOUBLE, root, comm);

    time = MPI_Wtime() - time;
    commT += time;
    time = MPI_Wtime();

    // in the root processor, merge the sorted chunks into one sorted array
    if(rank == 0) {
        for(int i=1; i<size; i++) {
            
            // merge chunks into tmpArray, and move it to array
            double * tmpArray = merge_array(i*n/size, array, n/size, array+i*n/size);
            for(int j=0;j<(i+1)*n/size;j++) array[j] = tmpArray[j];
        }
    }

    time = MPI_Wtime() - time;
    compT += time;
    time = MPI_Wtime();

    (*compTime) = compT;
    (*commTime) = commT;    

    return MPI_SUCCESS;
}

// merge sort an array a with n elements
void merge_sort(int n, double * a) {
    
    // index and pointer for sorted array
    int i;
    double * c;

    // deal with edge cases
    if(n<=1) return;

    // do comparison swap on pair of elements
    if(n==2) {
        if(n==2) swap(&a[0], &a[1]);
        return;
    }

    // split array in half and call sort recursively
    merge_sort(n/2, a);
    merge_sort(n-n/2, a+n/2);

    // merge array with sorted subarrays into temp array
    c=merge_array(n/2, a, n-n/2, a+n/2);

    // move sorted temp array into correct real location
    for(i=0;i<n;i++) a[i]=c[i];

    return;
}

// merge array a with n elements, with array b with m elements
// return the merged array
double * merge_array(int n, double * a, int m, double * b) {
    
    // indexes for two arrays and result array
    int i, j, k;

    // create array with size n+m elements
    double * c = (double *) calloc(n+m, sizeof(double));

    // iterate through a and b, adding smallest elem to c
    for(i=j=k=0;(i<n)&&(j<m);) {
        if(a[i]<=b[j]) c[k++]=a[i++];
        else c[k++] = b[j++];
    }
    // if reach the end of one array, add the rest of the other
    if(i==n) for(;j<m;) c[k++]=b[j++];
    else for(;i<n;) c[k++]=a[i++];

    return c;
}

// in-place swap two doubles
void swap(double * a, double * b) {
    // temp double for swapping
    double temp;

    // move a to temp during b to a, then temp to b
    temp = *a;
    *a = *b;
    *b = temp;
}