#pragma once

// Horner's Polynomial Evaluation
// (arr a_i, polynomial degree, value X)
double hornerEval(double* const arr, unsigned int n, double const X)
{
	double result = arr[n]; // b_n = a_n
	while (n-- > 0)
		result = arr[n] + X * result; // b_i = a_i + b_{i+1}
	return result;
}

// Horner's Polynomial Evaluation for Newton polynomial form
// (arr b_i, arr x_i, polynomial degree, value X)
double hornerNewtonEval(double* const b, double* const x, unsigned int n, double const X)
{
	double result = b[n]; // w_n = b_n
	while (n-- > 0)
		result = b[n] + (X - x[n]) * result; // w_i = b_i + (X - x_i) * w_{i+1}
	return result;
}

// Allows change from Newnton polynomial form to natural form
// (arr b_i, arr x_i, polynomial degree)
double* NewtonToNatural(double* const b, double* const x, unsigned int n)
{
	double* a = new double[n + 1];
	a[n] = b[n]; // a_n = b_n
	for (int i = n - 1; i >= 0; i--)
	{
		a[i] = b[i]; // max value of a_i
		for (int j = i; j < n; j++)
			a[j] -= x[i] * a[j + 1];
	}
	return a;
}

// Lagrange Interpolation Evaluation
// (arr x_i, arr y_i, no. of points, value X)
double LagrInterpolEval(double* const x, double* const y, unsigned int const n, double const X)
{
	double result = 0;
	for (int i = 0; i < n; i++) // y_i * l_i
	{
		double l = 1;
		for (int j = 0; j < n; j++) // l_i = (x - x_1)...(x - x_n) / (x_i - x_1)...(x_i-x_n)
			if(i!=j)
				l *= (X - x[j]) / (x[i] - x[j]);
		result += y[i] * l;
	}
	return result;
}