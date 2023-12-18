import math

from Qubit import *
from QubitSystem import *
import numpy as np
from scipy.linalg import block_diag
from computating.Algorithms import *
from computating.QuantumOperators import *

import itertools

def binary_fraction_to_decimal(binary_str):
    #binary_str = binary_str[::-1]
    decimal_result = sum(int(bit) * 2**(-i-1) for i, bit in enumerate(binary_str))
    return decimal_result

def estimate_probabilities(probabilities):
    n = len(probabilities)
    all_combinations = list(itertools.product([0, 1], repeat=n))

    results = {}

    max_probability = 0.0

    for combination in all_combinations:
        probability = 1.0
        for i, bit in enumerate(combination):
            probability *= probabilities[i] if bit == 1 else (1 - probabilities[i])
        results[combination] = probability

        if probability > max_probability:
            max_probability = probability

    # Находим порог вероятности (90% от максимальной)
    threshold = max_probability

    # Выбираем все строки с вероятностью в пределах [threshold, max_probability]
    selected_rows = ["".join(map(str, combination)) for combination, probability in results.items() if threshold <= probability <= max_probability]

    return selected_rows


#
# q1 = Qubit()
# q1.operator(X)
# q2 = Qubit()
#
# qs = QubitSystem(q1, q2)
#
# qs.operator(I, I)
#
# #print(qs.measure(0, KET_1))
# #print(qs.measure(1, KET_1))
#
# qq = Qubit()
# qq.operator(X)
# #print(TeleportFIM(qq).measure(2, qq.state))
# #print(TeleportPSIP(qq).measure(2, qq.state))
# #print(TeleportPSIM(qq).measure(2, qq.state))
# #print(TeleportFIP(qq).measure(2, qq.state))
#
# #print(qs.measure(0, KET_1))
# #print(qs.measure(1, KET_1))
# #print(qs.measure(2, KET_1))
#
# #q3 = Qubit()
# #q2.operator(X)
# #q2.operator(H)
# #print(q2.measure())
#
# #print(QubitSystem(q1, q2, q3).measure(1, KET_1))
#
# #res = BB84(10)
# #q = Qubit()
# #q.turn(np.pi / 2, X)
# #print(q.state)
# #qantum_CHSH(100)
# #classic_CHSH(100)
#
# # print(np.kron(I,I))
#
# print(Deutsch(II))
# print(Deutsch(np.kron(X, I)))
# print(Deutsch(CNOT))
# print(Deutsch(np.kron(X, I) @ CNOT @ np.kron(X, I)))
#
q1 = Qubit()
q2 = Qubit()
q3 = Qubit()
q4 = Qubit()
q5 = Qubit()
q1.operator(X)
q2.operator(X)
q3.operator(X)


#qs = QubitSystem(q1, q2, q3, q4, q5)
#qs.operator(make_oracle_gr2(5, 3))
#print(qs.measure_each(KET_1, 5))
# num = 87
# print(str(format(num, 'b')))
# oracle = make_oracle_bv(num)
#
# print(bernstein_vazirani(num, oracle))
q1 = Qubit()
q1.operator(H)
q2 = Qubit()
qs = QubitSystem(q1, q2)
#print(qs.measure_each(KET_P, 2))

oracle1 = make_oracle_grov(4, 9)

# dims = 4
# HH = H
# for j in range(dims-1):
#     HH = np.kron(H, HH)
# NNN = 2 * QubitSystem(dims).state @ QubitSystem(dims).state.transpose() - np.eye(2 ** dims)
#
# G = np.kron(HH @ NNN @ HH, np.eye(2))
#
# qs = QubitSystem(9)
#
# h = H
# for j in range(9 - 1):
#     h = np.kron(h, H)
# qs.operator(np.eye(2 ** 8), X)
# qs.operator(h)
#
# print(qs.measure_each(KET_1, 9))
#
# qs.operator(custom_CNOT(9, [3], G))
# qs.operator(custom_CNOT(9, [2], G @ G))
# qs.operator(custom_CNOT(9, [1], G @ G @ G @ G))
# qs.operator(custom_CNOT(9, [0], G @ G @ G @ G @ G @ G @ G @ G))
# qs.operator(RQFT(16), np.eye(32))
#qs.operator(h)
#print(qs.measure_each(KET_1, 9))

#res = phase_evaluation(4, np.kron(HH @ NNN @ HH, np.eye(2)) @ oracle1, 5)
#print(res)
#print("Theta Kirill: ", binary_fraction_to_decimal(res))

oracle2 = make_oracle_gr3(4, 9)
qs = QubitSystem(8)

h = H
for j in range(8 - 1):
    h = np.kron(h, H)
qs.operator(h)

#G = D(4) @ oracle2
#qs.operator(custom_CNOT(8, [3], G))
#qs.operator(custom_CNOT(8, [2], G @ G))
#qs.operator(custom_CNOT(8, [1], G @ G @ G @ G))
#qs.operator(custom_CNOT(8, [0], G @ G @ G @ G @ G @ G @ G @ G))
#qs.operator(QFT(16), np.eye(16))
#print(QFT(2))

#qs.operator(h)
#print(qs.measure_each(KET_1, 8))
# #print(RQFT(4))

t = 8
m = t - 2
oracle2 = make_oracle_gr3(4, 0, 1, 2, 3, 5, 7)
res = phase_evaluation(t, D(4) @ oracle2, 4)
#print(v)
#
#res = estimate_probabilities(v)
real_M = 6
dM = math.sqrt(2 * real_M * 16 * 2 ** (- m - 1)) * 2 ** (-m)
print("real M:", real_M)
print("Delta M:", dM)
for i in res:
    theta = binary_fraction_to_decimal(i)
    print("My theta:", theta)
    M = 2 * 16 * (math.sin(theta / 2) ** 2)
    if M != 0:
        print("M:", M)

        R = round(math.pi / 4 * math.sqrt(16 / M))
        print("R:", R)


# sum_m = 0
# sum_r = 0
# it = 100
# for j in range(it):
#     i = [np.random.choice([0, 1], p=[1 - prob, prob]) for prob in v]
#     theta = int("".join(map(str, i)), 2) / 2**t
#     print(theta)
#     #print(theta)
#     #print("My theta:", theta)
#     M = 2 * 16 * (math.sin(theta / 2) ** 2)
#     #print(M)
#     sum_m += M
#     #print(M, R)
# print(sum_m / it)
