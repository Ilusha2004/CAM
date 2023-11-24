#pragma once
#ifndef MEDHODS_H
#define MEDHODS_H

#include<vector>
#include<algorithm>
#include<cmath>

#include"GoodThings.h"

typedef std::vector<double> vector_d;

// метод простой итерации(метод Якоби)
vector_d methodJakob(const std::vector<vector_d>& matrix, const vector_d& solveVector, double epsilon = 10e-4)
{
    std::vector<vector_d> B(matrix.size(), vector_d(matrix.size(), 0.0));
    vector_d d(matrix.size(), 0.0);

    for (size_t i = 0; i < matrix.size(); i++)
    {
        d[i] = solveVector[i] / matrix[i][i];
        for (size_t j = 0; j < matrix.size(); j++)
        {
            B[i][j] = -matrix[i][j] / matrix[i][i];
        }
    }

    vector_d x_k = solveVector;
    vector_d temp(matrix.size(), 0.0);

    while (true)
    {
        vector_d x_k1(matrix.size(), 0.0);

        for (size_t i = 0; i < matrix.size(); i++)
        {
            for (size_t j = 0; j < matrix.size(); j++)
            {
                if (i != j)
                    x_k1[i] += B[i][j] * x_k[j];
            }
            x_k1[i] += d[i];

            temp[i] = std::abs(x_k1[i] - x_k[i]);
        }

        if (norm(temp) < epsilon)
        {
            break;
        }

        x_k = x_k1;
        x_k1.clear();

    }

    return x_k;

}

// Градиентный спуск
vector_d gradientSolve(const std::vector<vector_d>& matrix, const vector_d& resultVector, size_t am_iter = 35U)
{
    vector_d x_k1(matrix.size(), 0.0);
    vector_d x_k = resultVector;

    for (size_t k = 0; k < am_iter; k++)
    {

        vector_d r_k = VectorProduct(matrix, x_k);

        for (size_t i = 0; i < r_k.size(); i++)
        {
            r_k[i] -= resultVector[i];
        }

        double temp = dotProduct(r_k, r_k) / dotProduct(VectorProduct(matrix, r_k), r_k);

        for (size_t i = 0; i < r_k.size(); i++)
        {
            r_k[i] *= temp;
        }

        vector_d approxiamateResult = x_k;

        for (size_t i = 0; i < r_k.size(); i++)
        {
            approxiamateResult[i] -= r_k[i];
        }

        x_k1 = approxiamateResult;

        x_k = approxiamateResult;

    }

    return x_k1;
}

void swap_rows(std::vector<vector_d>& matrix, int row1, int row2, int n) {
    for (int i = 0; i <= n; i++) {
        double temp = matrix[row1][i];
        matrix[row1][i] = matrix[row2][i];
        matrix[row2][i] = temp;
    }
}

std::vector<vector_d> forwardElimination(const std::vector<vector_d>& Matrix)
{
    std::vector<vector_d> matrix(Matrix);
    int n = matrix.size();

    for (size_t i = 0; i < n; i++) {
        // Поиск главного элемента в текущей подматрице
        int max_row = i;
        double max_val = std::abs(matrix[i][i]);
        for (size_t j = i + 1; j < n; j++)
        {
            for (size_t k = i + 1; k <= n; k++)
            {
                if (std::abs(matrix[j][k]) > max_val)
                {
                    max_val = std::abs(matrix[j][k]);
                    max_row = j;
                }
            }
        }

        // Перестановка строк, если найденный главный элемент не находится на диагонали
        if (max_row != i)
        {
            swap_rows(matrix, i, max_row, n);
        }

        // Приведение матрицы к треугольному виду
        for (int j = i + 1; j < n; j++)
        {
            double factor = matrix[j][i] / matrix[i][i];
            for (int k = i; k <= n; k++)
            {
                matrix[j][k] -= factor * matrix[i][k];
            }
        }
    }

    return matrix;
}

vector_d backSubstitution(const std::vector<vector_d>& matrix)
{
    vector_d solution(matrix.size(), 0.0);
    int n = matrix.size();
    for (int i = n - 1; i >= 0; i--)
    {
        double sum = 0.0;

        for (int j = i + 1; j < n; j++)
        {
            sum += matrix[i][j] * solution[j];
        }

        solution[i] = (matrix[i][n] - sum) / matrix[i][i];
    }

    return solution;
}

vector_d gaussianElimination(const std::vector<vector_d>& matrix, const vector_d& resultVector)
{
    return backSubstitution(forwardElimination(ExtendedMatrix(matrix, resultVector)));
}

// Метод отражений
vector_d methodReflections(const std::vector<vector_d>& matrix, const vector_d& solveVector)
{
    int size = matrix.size();
    std::vector<vector_d> A(ExtendedMatrix(matrix, solveVector));
    std::vector<vector_d> E(size, vector_d(size, 0.0));
    vector_d s(size, 0.0);
    vector_d e(size, 0.0);

    for (size_t i = 0; i < size; i++)
    {
        E[i][i] = 1;
    }

    for (size_t k = 0; k < size; k++)
    {
        e[k] = 1;

        for (size_t i = 0; i < size; i++)
        {
            if (i >= k)
            {
                s[i] = A[i][k];
            }
        }

        double alpha = sign(arg(dotProduct(s, s)) - M_PI) * std::sqrt(dotProduct(s, s));

        vector_d s_t(size, 0.0);

        for (size_t i = 0; i < size; i++)
        {
            s_t[i] = s[i] - alpha * e[i];
        }

        double ka = 1 / std::sqrt(2 * alpha * alpha + 2 * std::abs(alpha) * dotProduct(s, e));

        vector_d w = s_t;

        for (size_t i = 0; i < size; i++)
        {
            w[i] *= ka;
        }

        vector_d new_b_k(size, 0.0);
        std::vector<vector_d> U(size, vector_d(size, 0.0));

        for (size_t i = 0; i < size; i++)
        {
            for (size_t k = 0; k < size; k++)
            {
                U[i][k] = E[i][k] - 2 * w[i] * w[k];
            }
        }

        std::vector<vector_d> A_1 = multiplyMatrices(U, A);

        A = A_1;

        e.clear();
        s.clear();
        s_t.clear();
        w.clear();
        U.clear();
        A_1.clear();
        e.resize(size);
        s.resize(size);
    }

    return backSubstitution(A);
}

//Определитель матрицы
double det(const std::vector<vector_d>& matrix)
{
    vector_d temp(matrix.size(), 0.0);
    std::vector<vector_d> tr = forwardElimination(ExtendedMatrix(matrix, temp));

    double result = 1.0;

    for (size_t i = 0; i < tr.size(); i++)
    {
        result *= tr[i][i];
    }

    return result;
}

std::vector<vector_d> InvertedMatrix(const std::vector<vector_d>& matrix)
{
    std::vector<vector_d> inverted_matix(matrix.size(), vector_d(matrix.size(), 0.0));

    for (size_t i = 0; i < matrix.size(); i++)
    {
        vector_d temp(matrix.size(), 0.0);
        temp[i] = 1;
        inverted_matix[i] = methodJakob(matrix, temp);
    }

    return inverted_matix;

}

// функция для вычисления числа обусловленности
double calculateConditionNumber(const std::vector<vector_d>& matrix)
{
    return MatrixNorm(InvertedMatrix(matrix)) * MatrixNorm(matrix);
}

#endif // !MEDHODS_H