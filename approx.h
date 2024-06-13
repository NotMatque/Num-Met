#pragma once
#include "linearEquations.h"
#include "helpful.h"

double* approx(double (*func)(double), double a, double b, unsigned int size)
{
	double* result = new double[size] {0};
	double** matrixbeg = new double* [size]; // Identity matrix
	for (int i = 0; i < size; i++)
	{
		matrixbeg[i] = new double[size] {0};
		matrixbeg[i][i] = 1;
	}
	
	double** matrix = gramSchmidt(matrixbeg, 0, 1, size); // Orthogonal matrix
	double norm;
	for (int i = 0; i < size; i++)
	{
		norm = sqrt(quadGaus4_2xHorner(matrix[i], size, matrix[i], size, a, b, 1.e5));
		for (int j = 0; j < size; j++)
			matrix[i][j] /= norm; // Orthonormal matrix
	}

	double* alpha = new double[size];
	for (int i = 0; i < size; i++)
		alpha[i] = quadGaus4_FuncxHorner(func, matrix[i], size, a, b, 1.e5);

	for(int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
			result[i] += alpha[j] * matrix[j][i];
	}

	return result;
}
