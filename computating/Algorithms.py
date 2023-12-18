import numpy as np

from Qubit import *
from QubitSystem import *
from computating.QuantumOperators import *


def QRNG():
    q = Qubit()
    q.operator(H)
    return q.measure()


def BB84(length):
    send_key = []
    recv_key = []
    n = 0
    it = 0
    while n < length:
        it += 1
        key1 = QRNG()
        send_basis = QRNG()
        q = Qubit()
        if key1:
            q.operator(X)
        if send_basis:
            q.operator(H)
        recv_basis = QRNG()
        basis = KET_P if recv_basis else KET_0
        key2 = q.basis_measure(basis)
        if send_basis == recv_basis:
            n += 1
            send_key.append(int(key1))
            recv_key.append(int(key2))
    return send_key, recv_key


def classic_CHSH(games):
    wins = 0
    for j in range(games):
        a = QRNG()
        b = QRNG()
        x = y = 0
        if (a and b) == (x != y):
            wins += 1
    print("Классическая стратегия - доля побед: " + str(100 * wins / games) + "%")


def qantum_CHSH(games):
    wins = 0
    for j in range(games):
        a = QRNG()
        b = QRNG()
        q_s = QubitSystem(Qubit(), Qubit())
        q_s.operator(H, I)
        q_s.operator(CNOT)
        op1 = I
        if a == 0:
            op1 = turn(np.pi / 2, Y)
        if b == 0:
            op2 = turn(np.pi / 4, Y)
        else:
            op2 = turn(3 * np.pi / 4, Y)
        q_s.operator(op1, op2)
        x, y = q_s.measure(0, KET_1), q_s.measure(1, KET_1)
        if ((a and b) == (x != y)):
            wins += 1
    print("Квантовая стратегия - доля побед: " + str(100 * wins / games) + "%")


def Deutsch(f):
    q1 = Qubit()
    q2 = Qubit()
    q2.operator(X)

    qs = QubitSystem(q1, q2)
    qs.operator(H, H)
    qs.operator(f)
    qs.operator(H, I)
    return qs.measure(0, KET_1)


def teleport_FIP(qubit):
    # Подготовка состояния Белла
    qs = QubitSystem(qubit, Qubit(), Qubit())
    qs.operator(I, H, I)
    qs.operator(I, CNOT)

    # Телепортация
    qs.operator(CNOT, I)
    qs.operator(H, I, I)
    qs.operator(I, CNOT)
    qs.operator(ZZ)

    return qs


def teleport_FIM(qubit):
    # Подготовка состояния Белла
    qs = QubitSystem(qubit, Qubit(), Qubit())
    qs.operator(I, X, I)
    qs.operator(I, H, I)
    qs.operator(I, CNOT)

    # Телепортация
    qs.operator(CNOT, I)
    qs.operator(H, I, Z)
    qs.operator(I, CNOT)
    qs.operator(ZZ)
    return qs


def teleport_PSIP(qubit):
    # Подготовка состояния Белла
    qs = QubitSystem(qubit, Qubit(), Qubit())
    qs.operator(I, H, I)
    qs.operator(I, NCNOT)

    # Телепортация
    qs.operator(CNOT, I)
    qs.operator(H, I, I)
    qs.operator(I, NCNOT)
    qs.operator(ZZ)
    return qs


def teleport_PSIM(qubit):
    # Подготовка состояния Белла
    qs = QubitSystem(qubit, Qubit(), Qubit())
    qs.operator(I, X, I)
    qs.operator(I, H, I)

    # Телепортация
    qs.operator(I, NCNOT)
    qs.operator(CNOT, I)
    qs.operator(Z, Z, I)
    qs.operator(H, I, I)
    qs.operator(I, NCNOT)
    qs.operator(ZZ)

    return qs


def bits_count(num):
    m = 0
    while num > 0:
        num //= 2
        m += 1
    return m


def make_oracle_bv(num):

    bits = bits_count(num)

    oracle = np.eye(2 ** (bits + 1))
    for i in range(0,bits):
        if (num & (1 << i)) != 0:
            oracle = oracle @ make_op(bits + 1, i + 1, bits + 1, X)

    return oracle


def bernstein_vazirani(number, oracle):
    bits = bits_count(number)

    qs = QubitSystem(bits + 1)

    op = I
    for i in range(1, bits + 1):
        if i == bits:
            op = np.kron(op, X)
        else:
            op = np.kron(op, I)
    qs.operator(op)

    h = H
    for j in range(bits):
        h = np.kron(H, h)
    qs.operator(h)
    qs.operator(oracle)
    qs.operator(h)

    return np.flip(qs.measure_each(KET_1, bits))


def simon(f, qs_size):

    qs = QubitSystem(qs_size)

    h = I
    for j in range(qs_size // 2 - 1):
        h = np.kron(I, h)
    for j in range(qs_size // 2):
        h = np.kron(H, h)

    qs.operator(h)
    qs.operator(f)
    qs.operator(h)

    return qs.measure_each(KET_0, qs_size // 2)


def make_oracle_gr(qs_size, *nums):
    oracle = np.eye(2 ** (qs_size + 1))
    for n in nums:
        oracle[2*n][2*n] = 0
        oracle[2*n+1][2*n+1] = 0
        oracle[2 * n + 1][2 * n] = 1
        oracle[2 * n][2 * n + 1] = 1
    return oracle


def make_oracle_gr2(qs_size, number):
    num = format(number, f'0{qs_size}b')
    bits = [index for index, bit in enumerate(reversed(num)) if bit == '1']
    num = ''.join(reversed(num))
    op = 1
    for i in range(qs_size):
        if num[i] == "1":
            op = np.kron(op, X)
        else:
            op = np.kron(op, I)

    return custom_CNOT(qs_size, bits) @ op

def make_oracle_grov(dims, *nums):
    oper = np.eye(2 ** (dims + 1))
    for n in nums:
        oper[2*n][2*n] = oper[2*n+1][2*n+1] = 0
        oper[2*n+1][2*n] = oper[2*n][2*n+1] = 1
    return oper


def make_oracle_gr3(qs_size, *numbers):
    oracle = np.eye(2 ** qs_size)
    for n in numbers:
        oracle[n][n] = -1

    return oracle


def grover(oracle, qs_size):

    qs = QubitSystem(qs_size)

    h = H
    for j in range(1, qs_size):
        h = np.kron(H, h)
    qs.operator(h)

    print("Начальное распределение")
    qs.probs("1" * qs_size, qs_size)

    it_count = round(np.pi/4*np.sqrt(qs_size ** 2))
    for n in range(1):
        qs.operator(oracle)
        qs.operator(D(qs_size))
        print("После ", n + 1, " итерации.")
        qs.probs("1"*qs_size, qs_size)


def operator_G(dims, upr, G):
    dimG = int(math.log2(G.shape[0]))
    matrix = np.eye(2 ** (dims + dimG), dtype=complex)
    combs = list(map("".join, product("01", repeat=dims - 1)))
    wombs = list(map("".join, product("01", repeat=dimG)))
    smalls = []
    for w in wombs:
        smalls.append(int(w, 2))
    l = np.arange(dims).tolist()
    l.remove(upr - 1)
    for c in combs:
        str = "0" * dims
    str = str[:upr - 1] + "1" + str[upr:]
    for i in range(len(c)):
        str = str[:l[i]] + c[i] + str[l[i] + 1:]
    bigs = []
    for w in wombs:
        bigs.append(int(str + w, 2))
    for i in range(2 ** dimG):
        for j in range(2 ** dimG):
            matrix[bigs[i]][bigs[j]] = G[smalls[i]][smalls[j]]
    return matrix

def phase_evaluation(t, oracle, o_size):

    q_size = t + o_size

    qs = QubitSystem(q_size)

    h = H
    for j in range(q_size - 1):
        h = np.kron(h, H)
    qs.operator(h)
    g = oracle
    for j in range(t):
        qs.operator(custom_CNOT(q_size, [j], g))
        g = g @ g
    qs.operator(RQFT(2 ** t), np.eye(2 ** o_size))
    # s = ""
    # for j in range(t):
    #     s = str(qs.prob_one_bit("0" * j + "1" + "0" * (t - j - 1) + "0" * o_size, t)) + s

    return qs.prob_return("1" * t + "0" * o_size, q_size)
    #return qs.measure_each(KET_1, t)
    #return s



