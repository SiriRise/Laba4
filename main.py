# Формируется матрица F следующим образом: скопировать в нее А
# и если в С сумма чисел, по периметру, меньше чем произведение чисел по диагонали, то поменять местами В и С симметрично,
# иначе В и Е поменять местами несимметрично. При этом матрица А не меняется.
# После чего если определитель матрицы А больше суммы диагональных элементов матрицы F, то вычисляется выражение: A-1*AT – K * F-1,
# иначе вычисляется выражение (AТ +G-FТ)*K, где G-нижняя треугольная матрица, полученная из А.
# Выводятся по мере формирования А, F и все матричные операции последовательно.

import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

try:
    k = int(input("Введите число являющееся коэффициентом при умножении: "))
    n = int(input("Введите число больше 3 которое являеться рзмером матрицы: "))
    while n <= 3:
        n = int(input("Введите число больше 3 "))
    np.set_printoptions(linewidth=1000)
    A = np.random.randint(-10.0, 10.0, (n, n))
    print("\nМатрица A\n", A)
    length = n // 2
    B = np.array(A[:length, :length])
    C = np.array(A[:length, length + n % 2:n])
    E = np.array(A[length + n % 2:n, length + n % 2:n])
    F = A.copy()
    print("\nМатрица F\n", F)
    print("\nПодматрица С\n", C)
    sumC = sum(C[0]) + sum(C[length - 1]) + sum([f[0] + f[length - 1] for f in C[1:length - 1]])
    diag = np.prod(np.diag(C)) * np.prod(np.diag(np.rot90(C)))
    print("Сумма чисел по периметру:", sumC, "\nПроизведение чисел по диагонали:", diag)
    if sumC < diag:
        print("Cумма чисел по периметру меньше чем произведение чисел по диагонали" , sumC , "<" , diag , "меняем симметрично B и C")
        F[:length, length + n % 2:n] = B[:length, ::-1]
        F[:length, :length] = C[:length, ::-1]
    else:
        print("Cумма чисел по периметру больше чем произведение чисел по диагонали" , sumC , ">" , diag , "меняем не симметрично B и Е")
        F[:length, :length] = E
        F[length + n % 2:n, length + n % 2:n] = B
    print("\nОтформатированная матрица F\n", F)
    det = np.linalg.det(A)
    trace = np.trace(F) + np.trace(np.rot90(F))
    print("Определитель матрицы А:" , int(det),"\nСумма диагональных элементов матрицы F:" , trace)
    if det > trace:
        print("Определитель матрицы А больше суммы диагональных элементов матрицы F", det ,">",trace, "вычисляем выражение A^(-1)*A^T-K*F^(-1)")
        result = np.linalg.inv(A) * np.transpose(A) - k * np.linalg.inv(F)
        print(result)
    else:
        print("Определитель матрицы А меньше суммы диагональных элементов матрицы F", det ,"<",trace, "вычисляем выражение (A^T+G-F^T)*K")
        G = np.tri(n) * A
        result = (A.transpose() + G - F.transpose()) * k
        print(result)
    print("\nМатрица, которая используется при построение графиков:\n", A)
    explode = [0] * (n - 1)
    explode.append(0.1)
    plt.title("Круговая диаграмма")
    try:
        sizes = [round(np.mean(abs(F[i, ::])) * 100, 1) for i in range(n)]
    except IndexError:
        sizes = [round(np.mean(abs(F[i, ::])) * 100, 1) for i in range(n)]
    plt.pie(sizes, labels=list(range(1, n + 1)), explode=explode, autopct='%1.1f%%', shadow=True)
    plt.show()

    plt.plot(A)
    plt.title("График")
    plt.ylabel("y axis")
    plt.xlabel("x axis")
    plt.show()

    sns.heatmap(A, cmap="Spectral", annot = True)
    plt.title("тепловая карта")
    plt.ylabel("номер строки")
    plt.xlabel("номер столбца")
    plt.show()

except ValueError:
    print("Программа завершена, введите число")
