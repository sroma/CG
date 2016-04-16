# CG
Conjugate gradient method: 
solve linear system Ax=b, A=A^T, (Au,u) > 0 for any vector u from R^N, N - dimension of the matrix A.
A - symmetric and positive-definite matrix.
x - unknown vector.
b - known vector.

Files:
1. matrix.operations.omp.f90 - module for different matrix and vector operations.
2. main.omp.f90 - CG method (OpenMP).
3. main.mpi.f90 - CG method (MPI).
4. CG.py - python code.
