#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "mpi.h"

void merge_sort(int n, double * a);
void swap(double * a, double * b);
double * merge_array(int n, double * a, int m, double * b);

int MPI_Exchange(int n, double * array, int rank1, int rank2, MPI_Comm comm, double * compTime, double * commTime);
int MPI_Sort_oddeven(int n, double * array, int root, MPI_Comm comm, double * compTime, double * commTime);
int isArraySorted(int n, double * array, int root, MPI_Comm comm, double * compTime, double * commTime);

int main(int argc, char * argv[]) {
    int size, rank;
    int n = 10000000, i;
    double m = 10.0;
    double commT, compT = 0;
    double * array;

    MPI_Status status;
    
    // initialise MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // allocate space for array
    array = (double *)calloc(n, sizeof(double));

    if(rank == 0) {
        srand(((unsigned)time(NULL)+rank));
        for(i=0;i<n;i++) array[i] = ((double)rand()/RAND_MAX)*m;
    }

    MPI_Sort_oddeven(n, array, 0, MPI_COMM_WORLD, &compT, &commT);

    printf("Rank %d: Comm.T: %lfs, Comp.T: %lfs\n", rank, commT, compT);

    MPI_Finalize();
}


int MPI_Sort_oddeven(int n, double * array, int root, MPI_Comm comm, double * compTime, double * commTime) {
    int size, rank;
    double compT, commT, time = 0;

    commT = (*commTime);
    compT = (*compTime);

    time = MPI_Wtime();

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    time = MPI_Wtime() - time;
    commT += time;
    time = MPI_Wtime();

    double * localA = (double *)calloc(n/size, sizeof(double));

    time = MPI_Wtime() - time;
    compT += time;
    time = MPI_Wtime();

    MPI_Scatter(array, n/size, MPI_DOUBLE, localA, n/size, MPI_DOUBLE, root, comm);

    time = MPI_Wtime() - time;
    commT += time;
    time = MPI_Wtime();

    merge_sort(n/size, localA);

    time = MPI_Wtime() - time;
    compT += time;
    time = MPI_Wtime();

    (*compTime) = compT;
    (*commTime) = commT;

    for(int step=0;step<size;step++){
        if((step+rank)%2==0){
            if(rank<size-1)MPI_Exchange(n/size, localA, rank, rank+1, comm, compTime, commTime);
            } 
            else {
                if(rank>0)MPI_Exchange(n/size, localA, rank-1, rank, comm, compTime, commTime);
            }

        // test is localA sorted
        if(isArraySorted(n/size, localA, root, comm, compTime, commTime)){
            break;
        }
    }

    MPI_Gather(localA, n/size, MPI_DOUBLE, array, n/size, MPI_DOUBLE, root, comm);

    time = MPI_Wtime() - time;
    compT += time;
    time = MPI_Wtime();

    (*compTime) = compT;
    (*commTime) = commT;

    return MPI_SUCCESS;
}


int isArraySorted(int n, double * array, int root, MPI_Comm comm, double * compTime, double * commTime){
    // get rank and size
    int size, rank;
    double compT, commT, time = 0;

    commT = (*commTime);
    compT = (*compTime);

    time = MPI_Wtime();

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    time = MPI_Wtime() - time;
    commT += time;
    time = MPI_Wtime();

    // gather the first and last elements of array
    double * first = (double *)calloc(size, sizeof(double));
    double * last = (double *)calloc(size, sizeof(double));

    time = MPI_Wtime() - time;
    compT += time;
    time = MPI_Wtime();

    MPI_Gather(&array[0], 1, MPI_DOUBLE, first, 1, MPI_DOUBLE, root, comm);
    MPI_Gather(&array[n-1], 1, MPI_DOUBLE, last, 1, MPI_DOUBLE, root, comm);

    time = MPI_Wtime() - time;
    commT += time;
    time = MPI_Wtime();

    // if root then test array is sorted
    int answer = 1;
    if(rank==root){
        for(int i=0;i<size-1;i++){
            if(last[i]>first[i+1]){
                answer = 0;
            }
        }
    }

    time = MPI_Wtime() - time;
    compT += time;
    time = MPI_Wtime();

    // Bcast answer
    MPI_Bcast(&answer, 1, MPI_INT, root, comm);

    time = MPI_Wtime() - time;
    compT += time;
    time = MPI_Wtime();

    (*compTime) = compT;
    (*commTime) = commT;

    return answer;
}


int MPI_Exchange(int n, double * array, int rank1, int rank2, MPI_Comm comm, double * compTime, double * commTime) {
    int rank, size, result, i, tag1 = 0, tag2 = 1;
    double commT, compT, time = 0;

    commT = (*commTime);
    compT = (*compTime);

    time = MPI_Wtime();

    double * b = ( double * ) calloc( n, sizeof( double ) );
    double * c;

    time = MPI_Wtime() - time;
    compT += time;
    time = MPI_Wtime();
    
    MPI_Status status;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    
    time = MPI_Wtime() - time;
    commT += time;
    time = MPI_Wtime();
    
    if(rank == rank1) {
        result = MPI_Send(&array[0], n, MPI_DOUBLE, rank2, tag1, comm);
        result = MPI_Recv(&b[0], n, MPI_DOUBLE, rank2, tag2, comm, &status);

        time = MPI_Wtime() - time;
        commT += time;
        time = MPI_Wtime();

        c = merge_array(n, array, n, b);
        for(i = 0; i < n; i++) {
            array[ i ] = c[ i ];
        }

        time = MPI_Wtime() - time;
        compT += time;
        time = MPI_Wtime();
    }
    else if( rank == rank2 ) {
        result = MPI_Recv(&b[ 0 ], n, MPI_DOUBLE, rank1, tag1, comm, &status);
        result = MPI_Send(&array[ 0 ], n, MPI_DOUBLE, rank1, tag2, comm);

        time = MPI_Wtime() - time;
        commT += time;
        time = MPI_Wtime();

        c = merge_array(n, array, n, b);
        for(i =0; i < n; i++) {
            array[i] = c[i + n];
        }

        time = MPI_Wtime() - time;
        compT += time;
        time = MPI_Wtime();
    }

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