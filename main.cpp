//TODO: linearEquations.h: Gauss-Crout

#include <iostream>

#include "helpful.h"
#include "interpolation.h"
#include "linearEquations.h"
#include "integrals.h"
#include "approx.h"

double func1(double X)
{
	return pow(X, 3) + pow(X, 2) - 3 * X - 3;
}
double func2(double X)
{
	return pow(X, 2) - 2;
}
double func3(double X)
{
	return sin(pow(X,2)) - pow(X,2);
}
double func4(double X)
{
	return sin(pow(X, 2)) - pow(X, 2) + 0.5;
}

double bissec(double(*func)(double), double beg, double end, double tolerance)
{
	if (func(beg) * func(end) >= 0)
		exit(-1);

	double temp; unsigned int loop = 0;
	while (loop < 1000)
	{
		temp = (beg + end) / 2;
		if (fabs(func(temp)) < tolerance)
			return temp;

		if (func(beg) * func(temp) < 0)
			end = temp;
		else
			beg = temp;

		loop++;
	}

	return beg;
}
double newtonRaphson(double(*func)(double), double X, double tolerance)
{
	// How to numerical dif?
	// f'(x) = (f(x+h) - f(x-h)) / 2h
	unsigned int loop = 0;
	double h = 1.e-12;
	while (loop++ < 1000)
	{
		if (fabs(func(X)) < tolerance)
			return X;
		else
			X -= func(X) / ((func(X + h) - func(X - h)) / (2 * h));
	}

	return X + 1;
}

int main()
{
	return 0;
}

