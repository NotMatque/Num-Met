#include <iostream>
#include "interpolation.h"

void printPoly(double* const arr, unsigned int const n)
{
    for (int i = 0; i <= n; i++)
        std::cout << arr[i] << "x^" << i << " ";
    std::cout << "\n";
}

void printMatrix(double** const arr, unsigned int columns, unsigned int rows)
{
    for (int i = 0; i < columns; i++)
    {
        for (int j = 0; j < rows; j++)
        {
            std::cout << arr[i][j] << " ";
        }
        std::cout << "\n";
    }
}

int main()
{
    
    return 0;
}