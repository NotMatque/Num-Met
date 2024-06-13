#pragma once

// Newton-Cotes formula: Trapezoidal rule
// (pointer to a func, begin, end, step)
double quadTrap(double (*func)(double), double a, double b, double h);

// Newton-Cotes formula: Simpson's rule
// (pointer to a func, begin, end, step)
double quadSimp(double (*func)(double), double a, double b, double h);

// Gauss-Legendre quadrature
// (pointer to a func, begin end, no. of steps)
double quadGaus4(double (*func)(double), double const a, double const b, unsigned int const steps);

// Gauss-Legendre quadrature for 2 polynomials
double quadGaus4_2xHorner(double* poly1, unsigned int const deg1, double* poly2, unsigned int const deg2, double a, double b, unsigned int steps);

// Gauss-Legendre quadrature for a func and a polynomial
double quadGaus4_FuncxHorner(double (*func)(double), double* poly, unsigned int const deg, double a, double b, unsigned int steps);
