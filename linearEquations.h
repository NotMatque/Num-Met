#pragma once
#include <math.h>
#include <iostream>
// TODO: Gauss-CROUT

// WIP
unsigned int findIndexOfMax(double* const arr, unsigned int const size);


// WIP
double* gaussElim(double** A, double* b, unsigned int const size);


// Doolitle's algorithm for LU decomposition
// (matrix A, size of matrix A, empty matrix L, empty matrix U)
void DoolittleLU(double** A, unsigned int size, double** L, double** U);
