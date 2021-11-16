import gurobipy as gp
from gurobipy import GRB
import numpy as np
import math

# константы
ROOT_2 = math.sqrt(2)


def test():
    # c^t u -> min, Au <= b, u >= 0
    A = [[1, 2], [3, 4]]
    b = [3, 7]
    c = [2, 3]

    A, b, c = np.array(A), np.array(b), np.array(c)

    m = gp.Model("solve")
    u = m.addMVar(shape=c.shape[0], vtype=GRB.CONTINUOUS, name="u")

    m.setObjective(c.T @ u, GRB.MINIMIZE)

    m.addConstr(A @ u >= b, "constr1")

    m.optimize()

    return u.X


# функция логнормального распределения
def lognorm(x, m, sigma):
    return 0.5 + 0.5 * math.erf((math.log(x) - m) / (ROOT_2 * sigma))


def main():
    # число заданий
    I = 30

    # сложности заданий
    w = [1.311, 3.254, 3.254, 3.254, 4.874, 5.368, 7.011, 7.217, 8.244, 9.636,
         4.132, 6.902, 2.121, 3.436, 2.456, 5.359, 6.902, 7.283, 7.815, 9.399,
         2, 2.418, 2.666, 3.653, 5.242, 5.547, 6.453, 7.194, 8.795, 3.657]

    # группы заданий
    A = [[1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0],
         [0, 1, 0], [0, 1, 0], [0, 1, 0], [0, 1, 0], [0, 1, 0], [0, 1, 0], [0, 1, 0], [0, 1, 0], [0, 1, 0], [0, 1, 0],
         [0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]]

    # суммарная сложность и число заданий
    c = 29.46
    k = 5

    # параметры логнормального распределения
    m = [3.51, 3.92, 3.89, 3.91, 4.22, 4.63, 5.67, 5.71, 6.13, 6.39,
         4.72, 5.87, 3.83, 3.91, 3.87, 5.13, 5.25, 5.71, 5.94, 6.27,
         3.65, 3.73, 3.87, 3.96, 4.84, 4.95, 5.53, 5.89, 6.18, 3.88]
    sigma = math.sqrt(0.31)

    # для вычисления квантилей
    quantile = 1.6449782

    # разбиение
    Lni = 5
    theta_l = []
    p_l = []

    for i in range(I):
        t = [math.exp(m[i] + ROOT_2 * sigma * (2 * j - Lni + 2) / (Lni - 2)) for j in range(Lni - 1)]

        theta_ll = [math.exp(m[i] - ROOT_2 * sigma * quantile)] + \
                   [(t[j] + t[j + 1]) * 0.5 for j in range(Lni - 2)] + \
                   [math.exp(m[i] + ROOT_2 * sigma * quantile)]
        p_ll = [lognorm(t[0], m[i], sigma)] + \
               [lognorm(t[j + 1], m[i], sigma) - lognorm(t[j], m[i], sigma) for j in range(Lni - 2)] + \
               [1 - lognorm(t[Lni - 2], m[i], sigma)]
        theta_l.append(theta_ll)
        p_l.append(p_ll)

    # реализации матрицы и их вероятности
    D = int(math.pow(Lni, I))
    print(D)
    thetas = []
    p = []


if __name__ == '__main__':
    main()
