# Parallel Sorting Using C and MPICH

This repo contains implementations of various parallel sorting algorithms using the MPI standard (MPICH library). Each algorithm is timed for sorting 10,000,000 random numbers, and prints the time taken for each processor - distinguishing between communication time and computation time.

### Running the Algorithms

To run the algorithms, the following steps should be followed:

1. Clone the repository

```bash
git clone https://github.com/jackomeara/parallel-sorts
```

2. Ensure required programs are installed

MPICH and C must be installed on your computer.
On windows make must be installed additionally.

3. Navigate to a directory and make the file

```bash
cd bucket_sort
make
```

4. Run the program (for example with 4 processors)

```bash
mpirun -np 4 sort
```

5. Optional - add other machines to use their processors

If you have other machines on the network that you can use, create a machines file with their addresses and add this as an argument. Step 4 will now look like

```bash
mpirun -np 4 -machinefile machines sort
```