//TODO: linearEquations.h: Gauss-Crout

#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>

#include "interpolation.h"
#include "linearEquations.h"
#include "integrals.h"

void printArr(double* const arr, unsigned int const size)
{
	for (int i = 0; i < size; i++)
		printf("%.3lf\t", arr[i]);
	std::cout << "\n";
}

void printArr2D(double** const arr, unsigned int rows, unsigned int columns)
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < columns; j++)
			printf("%.3lf\t", arr[i][j]);
		std::cout << "\n";
	}
}


double f1(double X) // 0 -> 4.5
{
	return pow(X, 2) * pow(sin(X), 3);
}

double f2(double X) // -2 -> 2
{
	return (X - 1) * exp(pow(X,2));
}

int main()
{   
	//std::cout << quadSimp(&f1, 0, 4.5, 1.e-5) << "\n";
	//std::cout << quadGaus4(&f1, 0, 4.5, 1.e5) << "\n";

	//std::cout << quadSimp(&f2, -2, 2, 1.e-5) << "\n";
	//std::cout << quadGaus4(&f2, -2, 2, 1.e5) << "\n";

	double poly1[] = { -2, 1, 2, 1 , 1}; double poly2[] = {10, -1, 1, 2};
	printArr(poly1, sizeof(poly1)/sizeof(double));

	//std::cout << quad_Horner4(poly1, 3, poly2, 3, -2.f, 2.f, 1.e3);
	return 0;
}

