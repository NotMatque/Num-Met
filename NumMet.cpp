#include <iostream>
#include "interpolation.h"
#include "linearEquations.h"

void printPoly(double* const arr, unsigned int const n)
{
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
			std::cout << arr[i][j] << "\t";
		std::cout << "\n";
	}
}

unsigned int findIndexOfMax(double* arr, unsigned int size)
{
	int indexMax = 0;
	for(int i = 1; i < size; i++)
		if (fabs(arr[i]) > fabs(arr[indexMax]))
			indexMax = i;
	return indexMax;
}


int indexOfMax(double** A, unsigned int n, int row)
{
	int maxPos = 0;
	double max = fabs(A[row][0]);
	for (int j = 0; j < n; j++)
		if (fabs(A[row][j]) > fabs(max))
		{
			max = fabs(A[row][j]);
			maxPos = j;
		}
	return maxPos;
}

double* gaussElim(double** A, double* b, unsigned int const n)
{
	double* x = new double[n] {0.0}; // result
	int* change = new int[n]; //additional array to record changes in matrix A
	for (int i = 0; i < n; i++)
	{
		change[i] = i;
		std::cout << change[i] << " ";
	}
	std::cout << "\n\nOriginal matrix A:\n";
	printMatrix(A, n, n);

	// Elimination
	for (int row = 0; row < n - 1; row++) // no. of last constant row
	{
		int indexMax = indexOfMax(A, n, row); // finding max element in this row 
		
		//int temp = change[row]; 
		//change[row] = change[indexMax];
		//change[indexMax] = temp;
		
		for (int j = row + 1; j < n; j++) // no. of current column
		{
			double a_ji = A[j][change[row]];
			double a_ii = A[row][change[row]];
			b[j] -= b[row] * a_ji / a_ii;
			for (int k = row; k < n; k++) // current element
			{
				A[j][change[k]] -= A[row][change[k]] * a_ji / a_ii;
			}
		}
	}
	
	std::cout << "\n\nMatrix A after elimination:\n";
	printMatrix(A, n, n);

	// Calculating values of every x
	for (int i = n - 1; i >= 0; i--)
	{
		x[change[i]] = b[i];
		for (int j = i + 1; j < n; j++)
			x[change[i]] -= A[i][change[j]] * x[change[j]];
		x[change[i]] /= A[i][change[i]];
	}

	return x;
}

/*
void LUDecomp(double** A, double** L, double** U, unsigned int const n)
{
	
	for (int row = 0; row < n; row++)
	{
		
		// L creation
		for (int col = 0; col < n; col++)
		{
			if (row > col)
				L[row][col] = 0;
			else if (row == col)
				L[row][col] = 1;
			else
			{
				L[row][col] = A[row][col];
				for (int i = 0; i < col; i++)
					L[row][col] -= L[row][i] * U[i][col];
				L[row][col] /= U[col][col];
			}
		}

		// U creation
		for (int col = 0; col < n; col++)
		{
			if (row < col)
				U = 0;
			else
			{
				U[row][col] = A[row][col];
				for (int i = 0; i < col; i++)
					U[row][col] -= L[row][i] * U[i][col];
			}
		}
		
		
		
	}
}*/


void LUDecomp(double** A, double** L, double** U, unsigned int const n) 
{
	for (int i = 0; i < n; i++) // row
	{
		for (int j = 0; j < n; j++) // L creation
		{
			if (j < i)
				L[j][i] = 0;
			else 
			{
				L[j][i] = A[j][i];
				for (int k = 0; k < i; k++) 
				{
					L[j][i] = L[j][i] - L[j][k] * U[k][i];
				}
			}
		}

		for (int j = 0; j < n; j++) // U creation
		{
			if (j < i)
				U[i][j] = 0;
			else if (j == i)
				U[i][j] = 1;
			else 
			{
				U[i][j] = A[i][j] / L[i][i];
				for (int k = 0; k < i; k++) 
				{
					U[i][j] = U[i][j] - ((L[i][k] * U[k][j]) / L[i][i]);
				}
			}
		}
	}
}

int main()
{   
	double** A = new double* [3];
	for (int i = 0; i < 3; i++)
		A[i] = new double[3];

	A[0][0] = 60;	A[0][1] = 30;	A[0][2] = 20;
	A[1][0] = 30;	A[1][1] = 20;	A[1][2] = 15;
	A[2][0] = 20;	A[2][1] = 15;	A[2][2] = 12;

	double** L = new double* [3];
	for (int i = 0; i < 3; i++)
		L[i] = new double[3];

	double** U = new double* [3];
	for (int i = 0; i < 3; i++)
		U[i] = new double[3];

	LUDecomp(A, L, U, 3);
	
	std::cout << "A:\n";
	printMatrix(A, 3, 3);

	std::cout << "L:\n";
	printMatrix(L, 3, 3);

	std::cout << "U:\n";
	printMatrix(U, 3, 3);

	return 0;
}

/*
double** A = new double* [5];
	for (int i = 0; i < 5; i++)
		A[i] = new double[5];

	A[0][0] = 1;    A[0][1] = -3;	A[0][2] = 4;    A[0][3] = 6.8;	A[0][4] = 9;
	A[1][0] = -3;	A[1][1] = 2;	A[1][2] = 4.6;	A[1][3] = 6.3;	A[1][4] = -10;
	A[2][0] = 2;    A[2][1] = -1;	A[2][2] = 2.8;	A[2][3] = -8.4;	A[2][4] = -5;
	A[3][0] = -6;	A[3][1] = 2;	A[3][2] = 7;	A[3][3] = -0.5;	A[3][4] = -0.9;
	A[4][0] = 5;    A[4][1] = -2;	A[4][2] = -0.5;	A[4][3] = 12;	A[4][4] = -4;

	double* b = new double[5] {74.64, -40.26, -2.32, 12.6, -8.9};

	double* result = gaussElim(A, b, 5);

	for (int i = 0; i < 5; i++)
	{
		std::cout << result[i] << "\n";
	}
*/