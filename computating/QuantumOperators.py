import math

import numpy as np
import scipy as sp
from itertools import product
from Qubit import *


I = np.array([
    [1, 0],
    [0, 1]
], dtype=complex)

II = np.kron(I, I)

H = np.array([
    [1, 1],
    [1, -1]
], dtype=complex) / np.sqrt(2)  # Оператор Адамара, как указан

X = np.array([
    [0, 1],
    [1, 0]
], dtype=complex)  # Оператор Паули Х, как указан

Y = np.array([
    [0, -1j],
    [1j, 0]
], dtype=complex)  # Оператор Паули Y, как указан

Z = np.array([
    [1, 0],
    [0, -1]
], dtype=complex)  # Оператор Паули Z, как указан

CNOT = np.array([
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 0, 1],
    [0, 0, 1, 0]
], dtype=complex)

NCNOT = np.array([
    [0, 1, 0, 0],
    [1, 0, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1]
], dtype=complex)

ZZ = np.array([[1, 0, 0, 0, 0, 0, 0, 0],
               [0, 1, 0, 0, 0, 0, 0, 0],
               [0, 0, 1, 0, 0, 0, 0, 0],
               [0, 0, 0, 1, 0, 0, 0, 0],
               [0, 0, 0, 0, 1, 0, 0, 0],
               [0, 0, 0, 0, 0, -1, 0, 0],
               [0, 0, 0, 0, 0, 0, 1, 0],
               [0, 0, 0, 0, 0, 0, 0, -1]], dtype=complex)


def D(n):
    d = np.ones((2**n, 2**n), dtype=complex)
    d = d / 2 ** (n - 1)
    d = d - np.eye(2**n, 2**n)
    return d

def projector(basis):

    return basis @ basis.transpose()


def make_op(dims, upr, chd, oper):
    matrix = np.eye(2 ** dims, dtype=complex)
    combs = list(map("".join, product("01", repeat=dims - 2)))
    l = np.arange(dims).tolist()
    l.remove(upr - 1)
    l.remove(chd - 1)
    for c in combs:
        str = "0" * dims
        str = str[:upr - 1] + "1" + str[upr:]
        for i in range(len(c)):
            str = str[:l[i]] + c[i] + str[l[i] + 1:]
        str0 = str
        str0 = str0[:chd - 1] + "0" + str0[chd:]
        num0 = int(str0, 2)
        str1 = str
        str1 = str1[:chd - 1] + "1" + str1[chd:]
        num1 = int(str1, 2)
        matrix[num0][num0] = oper[0][0]
        matrix[num0][num1] = oper[0][1]
        matrix[num1][num0] = oper[1][0]
        matrix[num1][num1] = oper[1][1]
    return matrix


def turn(angle, axe):
    return sp.linalg.expm(-1j * angle * axe / 2)


def toffoli(n):
    op = np.eye(2 ** n)
    op[2 ** n - 1][2 ** n - 1] = X[1][1]
    op[2 ** n - 2][2 ** n - 1] = X[0][1]
    op[2 ** n - 1][2 ** n - 2] = X[1][0]
    op[2 ** n - 2][2 ** n - 2] = X[0][0]
    return op


def custom_CNOT(qs_size, ctrl, X):
    combinations = []

    for i in range(2 ** len(ctrl)):
        bytes = format(i, f'0{len(ctrl)}b')
        combination = []
        for j in range(len(ctrl)):
            if bytes[j] == "1":
                combination.append(projector(KET_1))
            else:
                combination.append(projector(KET_0))
        combinations.append(combination)

    bCNOT = np.zeros((2**qs_size, 2**qs_size), dtype=complex)
    for i in range(len(combinations)):
        op = 1
        comb_index = 0
        for j in range(qs_size):
            if j == qs_size - 4 and i == len(combinations) - 1:
                op = np.kron(op, X)
                break
            else:
                if j in ctrl:
                    op = np.kron(op, combinations[i][comb_index])
                    comb_index += 1
                else:
                    op = np.kron(op, I)
        bCNOT = bCNOT + op
    return bCNOT

def RQFT(n):
    matrix = np.eye(n, dtype=complex)
    for i in range(n):
        for j in range(n):
            matrix[i][j] = (1 / math.sqrt(n) * np.power((math.cos(2 * math.pi / n) + 1j * math.sin(2 * math.pi / n)), -i * j))
    return matrix
