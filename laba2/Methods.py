import math
import numpy as np
import scipy as scipy

n = 10 # размерность матрицы А
c = 10
b = 100

def powMethod(A, eps=10e-9):
    At = A.transpose() # находим Аt
    A_s = At @ A # перемножаем А на Аt, теперь А – симметрическая
    yk = np.ones(n) # начальное приближение y_0
    y = np.dot(A_s, yk) # y_1
    l = y[0] / yk[0] # начальное lambda
    k = 1

    while (True): # итерационный процесс
         yk = np.dot(A_s, y)
         lk = yk[0] / y[0]
         yk /= max(yk) # нормируем

         if abs(lk - l) <= eps:
             break

         y = yk
         l = lk
         k += 1

    return k, yk, lk, y, l, A_s

def spinningMethod(A, eps=10e-9):
     At = A.transpose() # находим Аt
     A = np.dot(At, A) # перемножаем А на Аt, теперь А – симметрическая
     E, U = np.identity(n), np.identity(n)
     k = 0
     ak = A # копия А


     print("eigen A \n", np.linalg.eig(A))

     print("Матрица A")
     for e in A:
          print(*e)


     while (True): # итерационный процесс
          L = np.tril(ak)

          temp = np.absolute(ak - L) # получаем верхний треугольник, без диагонали
          sigma = sum(sum([abs(el) ** 2 for el in temp])) # считаем сумму квадратов недиагональных

          if (sigma <= eps): # условие итерирования
               break

          i, j = np.unravel_index(temp.argmax(), temp.shape) # индексы максимального
          alpha = math.atan(2 * ak[i][j] / (ak[i][i] - ak[j][j])) / 2 # считаем угол
          uk = np.identity(n)
          uk[i][i], uk[i][j], uk[j][j], uk[j][i] = np.cos(alpha), -np.sin(alpha), np.cos(alpha), np.sin(alpha)
          ak = np.dot(np.dot(uk.transpose(), ak), uk) # Ukt*A*Uk
          U = np.dot(U, uk) # U*Uk
          k += 1

     return ak, U, k # возвращаем матрицу А, матрицу U, количество итераций k


if __name__ == '__main__':
     A = np.array([c * (b - c) * np.random.random(n) for _ in range(n)])
     k, yk, lk, y, l, As = powMethod(A)
     r = np.dot(As, yk) - lk * yk # находим вектор невязки
     rnorm = np.linalg.norm(r, 1) # находим норму невязки
     print(r)

     ak, U, k = spinningMethod(A)
     l = ak.diagonal() # lambda для каждого коэффициента матрицы А
     lmax = max(ak.diagonal()) # max lambda
     print("собственные значения ", l)
     print(ak.diagonal().argmax())
     print(*U.transpose())
     x = U.transpose()[ak.diagonal().argmax()] # получаем x
     print(x)
     x /= max(x) # нормируем
     r = np.dot(A, x) - lmax * x # находим вектор невязки
     print(r)
     rnorm = np.linalg.norm(r, 1) # находим норму невязки
     print(rnorm)