#include <iostream>
#include <iomanip>

#include "interpolation.h"
#include "linearEquations.h"

void printPoly(double* const arr, unsigned int const n)
{
	std::cout << std::setprecision(5);
	for (int i = 0; i <= n; i++)
		std::cout << arr[i] << "x^" << i << " ";
	std::cout << "\n";
}

// TODO: change into pointer arythmetic
void printMatrix(double** const arr, unsigned int rows, unsigned int columns)
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < columns; j++)
			printf("%.3lf\t", arr[i][j]);
		std::cout << "\n";
	}
}

int main()
{   
	double** A = new double*[3];
	double** L = new double* [3];
	double** U = new double* [3];
	for (int i = 0; i < 3; i++) 
	{
		A[i] = new double[3];
		L[i] = new double[3];
		U[i] = new double[3];
	}

	A[0][0] = 3;	A[0][1] = 0;	A[0][2] = 1;
	A[1][0] = 0;	A[1][1] = -1;	A[1][2] = 3;
	A[2][0] = 1;	A[2][1] = 3;	A[2][2] = 0;

	DoolittleLU(A, 3, L, U);

	std::cout << "A:\n";
	printMatrix(A, 3, 3);
	std::cout << "L:\n";
	printMatrix(L, 3, 3);
	std::cout << "U:\n";
	printMatrix(U, 3, 3);
	return 0;
}

