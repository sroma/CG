# CG
Conjugate gradient method: 
solve linear system Ax=b, A=A^T, (Au,u) > 0 for any vector u from R^N, N - dimension of the matrix A.
A - symmetric and positive-definite matrix.
x - unknown vector.
b - known vector.

Theory: please read paper by Jonathan Richard Shewchuk "An Introduction to the Conjugate Gradient Method Without the Agonizing Pain" on http://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf

Files:
1. matrix.operations.omp.f90 - module for different matrix and vector operations.
2. main.omp.f90 - CG method (OpenMP).
3. main.mpi.f90 - CG method (MPI).
4. CG.py - python code.
