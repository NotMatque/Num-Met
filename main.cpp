//TODO: linearEquations.h: Gauss-Crout

#include <iostream>
#include <iomanip>

#include "interpolation.h"
#include "linearEquations.h"
#include "integrals.h"

void printArr(double* const arr, unsigned int const size)
{
	for (int i = 0; i < size; i++)
		printf("%.6lf\t", arr[i]);
	std::cout << "\n";
}

void printArr2D(double** const arr, unsigned int rows, unsigned int columns)
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < columns; j++)
			printf("%.6lf\t", arr[i][j]);
		std::cout << "\n";
	}
}

double** gramSchmidt(double** F, double a, double b, unsigned int n)
{
	double** G = (double**)malloc(5 * sizeof(double*));
	for (int i = 0; i < n; i++)
		G[i] = (double*)malloc(5 * sizeof(double));

	for (int i = 0; i < n; i++)
	{
		// Copy row
		for (int k = 0; k < n; k++)
			G[i][k] = F[i][k];

		for (int j = 0; j < i; j++)
		{
			double ratio = quad_Horner4(F[i], n, G[j], n, a, b, 1.e5) / quad_Horner4(G[j], n, G[j], n, a, b, 1.e5);
			for (int k = 0; k < n; k++)
			{
				G[i][k] -= ratio * G[j][k];
			}
		}
	}
	
	return G;
}

int main()
{   
	double** F = new double* [5];
	for (int i = 0; i < 5; i++)
		F[i] = new double[5];

	F[0][0] = 1; F[0][1] = 0; F[0][2] = 0; F[0][3] = 0; F[0][4] = 0;
	F[1][0] = 0; F[1][1] = 1; F[1][2] = 0; F[1][3] = 0; F[1][4] = 0;
	F[2][0] = 0; F[2][1] = 0; F[2][2] = 1; F[2][3] = 0; F[2][4] = 0;
	F[3][0] = 0; F[3][1] = 0; F[3][2] = 0; F[3][3] = 1; F[3][4] = 0;
	F[4][0] = 0; F[4][1] = 0; F[4][2] = 0; F[4][3] = 0; F[4][4] = 1;

	printArr2D(F, 5, 5);
	std::cout << "\n";
	double** G = gramSchmidt(F, 0, 1, 5);
	printArr2D(G, 5, 5);

	return 0;
}

