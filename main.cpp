//TODO: linearEquations.h: Gauss-Crout

#include <iostream>

#include "helpful.h"
#include "interpolation.h"
#include "linearEquations.h"
#include "integrals.h"
#include "approx.h"

double funcTest(double X)
{
	return sin(-X) + exp(-X) - pow(X, 3);
}

double func(double T)
{
	double alpha = -1.e-12;
	double beta = 0;
	return alpha * (pow(T,4) - beta);
}

// (pointer to a func, step, beg value, end case)
double me(double(*func)(double), double step, double T_0, double tau)
{
	double result = T_0;
	for (int i = 0; i < tau / step; i++)
	{
		result = result + step * func(result);
	}
	return result;
}

// (pointer to a func, step, beg value, end case)
double me2(double(*func)(double), double step, double T_0, double tau)
{
	double result = T_0;
	for (int i = 0; i < tau / step; i++)
	{
		result = result + step * func(result + 0.5 * step * func(result));
	}
	return result;
}

// (pointer to a func, step, beg value, end case)
double me3(double(*func)(double), double step, double T_0, double tau)
{
	double result = T_0;
	for (int i = 0; i < tau / step; i++)
	{
		result = result + 0.5 * step * (func(result) + func(result + step * func(result))) ;
	}
	return result;
}

double k1(double Y, double step)
{
	return func(Y);
}
double k2(double Y, double step)
{
	return func(Y + 0.5 * step * k1(Y, step));
}
double k3(double Y, double step)
{
	return func(Y + 0.5 * step * k2(Y, step));
}
double k4(double Y, double step)
{
	return func(Y + step * k3(Y, step));
}
double phi(double Y, double step)
{
	return 1/6 * (k1(Y, step) + 2 * k2(Y, step) + 2 * k3(Y, step) + k4(Y, step));
}

double rangeKutt(double (*func)(double), double step, double T_0, double tau)
{
	double result = T_0;
	for (int i = 0; i < tau / step; i++)
	{
		result = result + step * phi(result, step);
	}
	return result;
}

// 877.7543 <- wynik poprawny
int main()
{
	double value = 877.7542611;
	printf("%.7lf\n", value);
	double temp = me(&func, 0.1, 1200, 300);
	printf("%.7lf, delta = %.7lf\n", temp, temp - value);
	temp = me2(&func, 0.1, 1200, 300);
	printf("%.7lf, delta = %.7lf\n", temp, temp - value);
	temp = me3(&func, 0.1, 1200, 300);
	printf("%.7lf, delta = %.7lf\n", temp, temp - value);
	temp = rangeKutt(&func, 0.1, 1200, 300);
	printf("%.7lf, delta = %.7lf\n", temp, temp - value);
	return 0;
}

