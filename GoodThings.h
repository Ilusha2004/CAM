#pragma once
#ifndef GOODTHINGS_H
#define GOODTHINGS_H

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <random>
#include <algorithm>
#include <string>
#include <sstream>

typedef std::vector<long double> vector_d;

// Функция для вычисления нормы вектора
long double vectorNorm(const vector_d& x)
{
    long double norm = 0.0;
    for (auto value : x)
    {
        norm += std::abs(value);
    }
    return norm;
}

long double MatrixNorm(const std::vector<vector_d>& matrix)
{
    long double sum = 0.0;

    for (auto& row : matrix)
    {
        for (auto element : row)
        {
            sum += std::abs(element);
        }
    }
    return sum;
}

// Построение расширенной матрицы
std::vector<vector_d> ExtendedMatrix(const std::vector<vector_d>& matrix, const vector_d& resultVector)
{
    std::vector<vector_d> matrix_k(matrix.size(), vector_d(matrix.size() + 1, 0.0));

    for (size_t i = 0; i < matrix_k.size(); i++)
    {
        for (size_t j = 0; j < matrix_k.size() + 1; j++)
        {
            if (j == matrix_k.size())
            {
                matrix_k[i][j] = resultVector[i];
            }
            else
            {
                matrix_k[i][j] = matrix[i][j];
            }
        }
    }

    return matrix_k;
}

// Функция для вычисления вектора невязки
vector_d calculateResidualVector(const std::vector<vector_d>& matrix, const vector_d& solutionVector, const vector_d& resultVector)
{
    int size = matrix.size();
    vector_d residualVector(size, 0.0);

    for (size_t i = 0; i < size; i++)
    {
        double sum = 0.0;

        for (size_t j = 0; j < size; j++)
        {
            sum += matrix[i][j] * solutionVector[j];
        }

        residualVector[i] = resultVector[i] - sum;
    }

    return residualVector;
}

void printMatrix(const std::vector<vector_d>& matrix)
{
    for (auto& row : matrix)
    {
        for (auto element : row)
        {
            std::cout << element << " " << std::setprecision(20);
        }
        std::cout << std::endl;
    }
}

// Функция для вывода вектора на экран
void printVector(const vector_d& vector)
{
    int size = vector.size();

    for (int i = 0; i < size; i++)
    {
        std::cout << std::fixed << std::setprecision(10) << vector[i] << " ";
    }

    std::cout << std::endl;
}

void printSolution(const vector_d& solution)
{
    for (int i = 0; i < solution.size(); ++i)
    {
        std::cout << "x[" << i << "] = " << solution[i] << std::endl;
    }
}

long double norm(const vector_d& vec)
{
    return *std::max_element(vec.begin(), vec.end(), [](double a, double b)
    {
        return std::abs(a) < std::abs(b);
    });
}

// Функция для умножения матрицы на вектор
vector_d VectorProduct(const std::vector<vector_d>& matrix, const vector_d& vector)
{
    int size = matrix.size();
    vector_d result(size, 0.0);

    for (size_t i = 0; i < size; i++)
    {
        for (size_t j = 0; j < size; j++)
        {
            result[i] += matrix[i][j] * vector[j];
        }
    }

    return result;
}

long double dotProduct(const vector_d leftVector, const vector_d rightVector, int index = 0)
{
    long double result = 0.0;

    for (size_t i = index; i < leftVector.size(); i++)
    {
        result += leftVector[i] * rightVector[i];
    }

    return result;
}

std::vector<vector_d> multiplyMatrices(const std::vector<vector_d>& matrix1, const std::vector<vector_d>& matrix2)
{
    int rows1 = matrix1.size();
    int cols1 = matrix1[0].size();
    int cols2 = matrix2[0].size();

    std::vector<vector_d> result(rows1, vector_d(cols2, 0));

    if (cols1 != matrix2.size())
    {
        std::cout << "Невозможно выполнить перемножение матриц!" << std::endl;
        return result;
    }

    for (size_t i = 0; i < rows1; i++)
    {
        for (size_t j = 0; j < cols2; j++)
        {
            for (size_t k = 0; k < cols1; k++)
            {
                result[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }

    return result;
}

template <class Value>
Value arg(Value number)
{
    if (number > 0.)
    {
        return 0;
    }
    else
    {
        return M_PI;
    }
}

template <class Value>
Value sign(Value number)
{
    if (number == 0.)
    {
        return 0;
    }
    if (number > 0.)
    {
        return 1;
    }
    else
    {
        return -1;
    }
}


#endif //!GOODTHINGS_H