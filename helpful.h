#pragma once
#include <iostream>
void inline printArr2D(double** const arr, unsigned int rows, unsigned int columns)
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < columns; j++)
			printf("%.6lf\t", arr[i][j]);
		std::cout << "\n";
	}
}

void inline printArr(double* const arr, unsigned int const size)
{
	for (int i = 0; i < size; i++)
		printf("%.6lf\t", arr[i]);
	std::cout << "\n";
}