#pragma once

// Newton-Cotes formula: Trapezoidal rule
// (pointer to a func, begin, end, step)
double quadTrap(double (*fun)(double), double a, double b, double h)
{
	double sum = 0;

	for (double i = a; i <= b; i += h)
		sum += h / 2 * (fun(i) + fun(i + h));

	return sum;
}

// Newton-Cotes formula: Simpson's rule
// (pointer to a func, begin, end, step)
double quadSimp(double (*fun)(double), double a, double b, double h)
{
	double sum = 0;

	for (double i = a; i <= b; i += h)
		sum += h / 6 * (fun(i) + 4 * fun((i + i + h) / 2) + fun(i + h));

	return sum;
}

// Gauss-Legendre quadrature
// (pointer to a func, begin end, no. of steps)
double quadGaus4(double (*func)(double), double const a, double const b, unsigned int const steps)
{
	double result = 0;
	double X[5] = { -0.90618f, -0.538469f, 0, 0.538469f, 0.90618f};
	double A[5] = {  0.236927, 0.478629, 0.568889, 0.478629, 0.236927};

	double dx = (b - a) / steps; // Size of 1 step

	for (int i = 0; i < steps; i++)
	{
		double x1 = i * dx + a; // begining
		double x2 = (i + 1) * dx + a; // end

		for (int j = 0; j < 5; j++)
			result += A[j] * func((x2 - x1) / 2 * X[j] + ((x1 + x2) / 2));
	}
	return result * dx / 2;
}

// Temporarily
double quad_Horner4(double* poly1, unsigned int const deg1, double* poly2, unsigned int const deg2, double a, double b, unsigned int steps)
{
	double result = 0;
	double X[5] = { -0.90618f, -0.538469f, 0, 0.538469f, 0.90618f };
	double A[5] = { 0.236927, 0.478629, 0.568889, 0.478629, 0.236927 };

	double dx = (b - a) / steps; // Size of 1 step

	for (int i = 0; i < steps; i++)
	{
		double x1 = i * dx + a; // begining
		double x2 = (i + 1) * dx + a; // end
		
		for (int j = 0; j < 5; j++)
		{
			double arg = (x2 - x1) / 2 * X[j] + ((x1 + x2) / 2);
			result += A[j] * hornerEval(poly1, deg1, arg) * hornerEval(poly2, deg2, arg);
		}
			
	}
	return result * dx / 2;
}