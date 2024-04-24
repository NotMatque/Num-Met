#pragma once
// TODO: Gauss-CROUT

// WIP
unsigned int findIndexOfMax(double* const arr, unsigned int const size)
{
	unsigned int iMax = 0;

	for (int i = 1; i < size; i++)
		if (fabs(arr[i]) > fabs(arr[iMax]))
			iMax = i;

	return iMax;
}

// WIP
double* gaussElim(double** A, double* b, unsigned int const size)
{
	double* x = new double[size];

	unsigned int* change = new unsigned int[size];
	for (int i = 0; i < size; i++)
		change[i] = i;

	// Elimination
	for (int i = 0; i < size - 1; i++) // Last const row
	{
		int iMax = findIndexOfMax(A[i], size);

		//unsigned int temp = change[i];
		//change[i] = iMax;
		//change[iMax] = temp;

		for (int j = i + 1; j < size; j++) // Current row
		{
			double temp = A[j][change[i]] / A[i][change[i]];
			for (int k = i; k < size; k++) // Current column
				A[j][change[k]] -= A[i][change[k]] * temp;
			b[j] -= b[i] * temp;
		}
	}

	// Calculating values of x
	for (int i = size - 1; i >= 0; i--)
	{
		x[change[i]] = b[i];
		for (int j = i + 1; j < size; j++)
			x[change[i]] -= A[i][change[j]] * x[change[j]];
		x[change[i]] /= A[i][change[i]];
	}

	for (int i = 0; i < size; i++)
	{
		std::cout << change[i] << " ";
	}

	delete[] change;
	return x;
}

// Doolitle's algorithm for LU decomposition
// (matrix A, size of matrix A, empty matrix L, empty matrix U)
void DoolittleLU(double** A, unsigned int size, double** L, double** U)
{
	// Preparing L and U
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			if (i == j)
				L[i][j] = 1;
			else
				L[i][j] = 0;

			U[i][j] = 0;
		}
	}

	// Checking for 0s on the diagonal of A
	/*
	unsigned int* change = new unsigned int[size];
	for (int i = 0; i < size; i++)
		change[i] = i; 
	
	for (int i = 0; i < size; i++)
	{
		if (A[change[i]][i] == 0)
		{
			int temp = change[i];
			change[i] = size % (i + 1);
			change[size % (i + 1)] = temp;
			i = 0;
		}
	}
	*/
	for (int i = 0; i < size; i++)
	{
		if (A[i][i] == 0)
		{
			// Copy values to next row
			for (int j = 0; j < size; j++)
			{
				double temp = A[i][j];
				A[i][j] = A[size % (i + 1)][j];
				A[size % (i + 1)][j] = temp;
			}
			i = 0;
		}
	}

	// Doolittle's algorithm
	for (int j = 0; j < size; j++) // j = col; i = row
	{
		if (A[j][j] == 0)
			exit(0);

		for (int i = 0; i <= j; i++) // U
		{
			U[i][j] = A[i][j];
			for (int k = 0; k < i; k++)
				U[i][j] -= L[i][k] * U[k][j];
		}
		for (int i = j + 1; i < size; i++) // L
		{
			L[i][j] = A[i][j] / U[j][j];
			for (int k = 0; k < j; k++)
				L[i][j] -= (L[i][k] * U[k][j]) / U[j][j];
		}
	}
}