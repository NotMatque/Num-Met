#pragma once

// Horner's Polynomial Evaluation
// (arr a_i, polynomial degree, value X)
double hornerEval(double* const arr, unsigned int n, double const X);

// Horner's Polynomial Evaluation for Newton polynomial form
// (arr b_i, arr x_i, polynomial degree, value X)
double hornerNewtonEval(double* const b, double* const x, unsigned int n, double const X);

// Allows change from Newton polynomial form to natural form
// (arr b_i, arr x_i, polynomial degree)
double* NewtonToNatural(double* const b, double* const x, unsigned int n);

// Evaluation by Lagrange Interpolation
// (arr x_i, arr y_i, no. of points, value X)
double LagrInterpolEval(double* const x, double* const y, unsigned int const n, double const X);