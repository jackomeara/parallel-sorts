#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mpi.h"

void merge_sort(int n, double * a);
double * merge_array(int n, double * a, int m, double * b);
void swap(double * a, double * b);

int MPI_Sort_bucket(int n, double * a, double max, int root, MPI_Comm comm, double * compTime, double * commTime);

int main(int argc, char * argv[]) {
    // time counters
    double commT, compT = 0;

    // initialise variables
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

    MPI_Sort_bucket(n, array, m, 0, MPI_COMM_WORLD, &compT, &commT);

    printf("Rank %d: Comm.T: %lfs, Comp.T: %lfs\n", rank, commT, compT);

    MPI_Finalize();
}

int MPI_Sort_bucket(int n, double * array, double max, int root, MPI_Comm comm, double * compTime, double * commTime) {
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

    double * localBucket = (double *)calloc(n, sizeof(double));

    time = MPI_Wtime() - time;
    compT += time;
    time = MPI_Wtime();


    // bcast array
    int rc = MPI_Bcast(array, n, MPI_DOUBLE, root, comm);
    if(rc!=MPI_SUCCESS) return rc;

    time = MPI_Wtime() - time;
    commT += time;
    time = MPI_Wtime();

    int count = 0;
    for(int i=0;i<n;i++) {
        if(array[i] >= rank*max/size && array[i] < (rank+1)*max/size) {
            localBucket[count++] = array[i];
        }
    }

    merge_sort(count, localBucket);
    int * counts = (int *)calloc(size, sizeof(int));
    int * displacements = (int *)calloc(size, sizeof(int));

    time = MPI_Wtime() - time;
    compT += time;
    time = MPI_Wtime();

    MPI_Gather(&count, 1, MPI_INT, counts, 1, MPI_INT, root, comm);

    time = MPI_Wtime() - time;
    commT += time;
    time = MPI_Wtime();


    if(rank == root){    
        displacements[0] = 0;
        for(int i=1;i<size;i++){
            displacements[i] = displacements[i-1]+counts[i];
            }
    }

    time = MPI_Wtime() - time;
    compT += time;
    time = MPI_Wtime();
    
    rc = MPI_Gatherv(localBucket, count, MPI_DOUBLE, array, counts, displacements,
    MPI_DOUBLE, root, comm);
    if(rc != MPI_SUCCESS)return rc;

    time = MPI_Wtime() - time;
    commT += time;
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