#include <iostream>

#include "Methods.h"
#include "GoodThings.h"

typedef std::vector<long double> vector_d;

//---------------------------------------------------------------------------------------------------------------------------
//
// Функция для генерации случайного вещественного числа с двумя знаками после запятой
double generateRandomDouble()
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(-100.00, 100.00);
    return std::round(dis(gen) * 100.00) / 100;
}

// Функция для генерации случайной квадратной матрицы с диагональным доминированием
std::vector<vector_d> generateRandomMatrix(int size)
{
    std::vector<vector_d> matrix(size, vector_d(size, 0.0));

    for (size_t i = 0; i < size; i++)
    {
        long double sum = 0.0;

        for (size_t j = 0; j < size; j++)
        {
            if (i != j)
            {
                matrix[i][j] = generateRandomDouble();
                sum += std::abs(matrix[i][j]);
            }
        }

        matrix[i][i] = sum + std::abs(generateRandomDouble());
    }

    return matrix;
}

// Функция для генерации случайного вектора решений
vector_d generateRandomSolutionVector(int size)
{
    vector_d vector(size, 0.0);

    for (size_t i = 0; i < size; i++)
    {
        vector[i] = generateRandomDouble();
    }

    return vector;
}
//
//---------------------------------------------------------------------------------------------------------------------------

int main(int arc, char** argv)
{
    int size = 3; // Размерность матрицы и вектора

    // Генерация матрицы с диагональным доминированием
    std::vector<vector_d> matrix = generateRandomMatrix(size);

    // Генерация вектора решений
    vector_d solutionVector = generateRandomSolutionVector(size);

    // Вычисление вектора свободных членов
    vector_d resultVector = VectorProduct(matrix, solutionVector);

    // Вывод сгенерированных данных
    std::cout << "Сгенерированная матрица:" << std::endl;
    printMatrix(matrix);

    std::cout << "Сгенерированный вектор решений:" << std::endl;
    printVector(solutionVector);

    std::cout << "Вычисленный вектор свободных членов:" << std::endl;
    printVector(resultVector);

    std::cout << "Нахождение обратной матрицы: " << std:: endl;

    // Нахождение обратной матрицы
    printMatrix(InvertedMatrix(matrix));

    std::cout << "Решение методом Якоби:" << std::endl;

    vector_d JacobVector = methodJakob(matrix, resultVector);
    printSolution(JacobVector);

    std::cout << "Решение методом отражений:" << std::endl;

    vector_d reflectVector = methodReflections(matrix, resultVector);
    printSolution(reflectVector);

    std::cout << "Решение методом градиентного спуска: " << std::endl;

    vector_d gradVector = gradientSolve(matrix, resultVector);
    printSolution(gradVector);

    std::cout << "Решение методом Гаусса с выбором главного элемента по матрице: " << std::endl;

    vector_d gauseVector = gaussianElimination(matrix, resultVector);
    printSolution(gauseVector);

    // Вычисление определителя матрицы
    std::cout << "Определитель: " << det(matrix) << std::endl;

    printMatrix(multiplyMatrices(matrix, InvertedMatrix(matrix)));
    // Нахождение числа обусловленности
    std::cout << "Число обусловленности: " << calculateConditionNumber(matrix) << std::endl;

    std::cout << "Вектор невязки: " << std::endl;
    printVector(calculateResidualVector(matrix, JacobVector, resultVector));
    printVector(calculateResidualVector(matrix, reflectVector, resultVector));
    printVector(calculateResidualVector(matrix, gradVector, resultVector));
    printVector(calculateResidualVector(matrix, gauseVector, resultVector));

}