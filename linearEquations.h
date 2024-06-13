#pragma once
#include <math.h>
#include <iostream>
#include "integrals.h"
// TODO: Gauss-CROUT

// WIP
unsigned int findIndexOfMax(double* const arr, unsigned int const size);


// WIP
double* gaussElim(double** A, double* b, unsigned int const size);


// Doolitle's algorithm for LU decomposition
// (matrix A, size of matrix A, empty matrix L, empty matrix U)
void DoolittleLU(double** A, unsigned int size, double** L, double** U);

// Gram-Schmidt process for orthogonalisation
// Orthogonalises matrix F
// (matrix F, integral bound begining, integral bound end, size of matrix F)
double** gramSchmidt(double** F, double a, double b, unsigned int n);