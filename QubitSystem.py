from computating.QuantumOperators import *
from Qubit import *


class QubitSystem:

    def __init__(self, *qubits):
        if isinstance(qubits[0], int):
            self.state = 1
            for i in range(0, qubits[0]):
                self.state = sp.kron(self.state, KET_0)
        else:
            self.state = qubits[0].state
            for i in range(1, len(qubits)):
                self.state = sp.kron(self.state, qubits[i].state)

    def operator(self, *op):

        operator = op[0]

        for i in range(1, len(op)):
            operator = sp.kron(operator, op[i])

        self.state = operator @ self.state

    def measure(self, pos, basis):
        proj = basis @ basis.transpose()
        operator = 1

        for i in range(int(np.log2(self.state.size))):
            if i == pos:
                operator = sp.kron(operator, proj)
            else:
                operator = sp.kron(operator, I)

        result = operator @ np.abs(self.state) @ np.abs(self.state).conjugate().transpose()
        p = np.trace(result)
        return round(np.abs(p), 3)


    def measure_each(self, basis, bits):
        measures = []
        for i in range(bits):
            measures.append(self.measure(i, basis))
        return measures

    def probs(self, selection, size):
        sbits = 0
        chbits = []
        for i in range(0, len(selection)):
            if (selection[i] == '1'):
                sbits += 1
                chbits.append(i)

        all_combs = list(map("".join, product("01", repeat=size)))
        s_combs = list(map("".join, product("01", repeat=sbits)))
        prob = np.zeros(len(s_combs), dtype=complex)
        for i in range(0, len(s_combs)):
            for j in range(0, len(all_combs)):
                Z = True
                for k in range(0, sbits):
                    if all_combs[j][chbits[k]] != s_combs[i][k]:
                        Z = False
                if Z:
                    prob[i] += abs(self.state[j][0]) * abs(self.state[j][0])
        for i in range(0, len(s_combs)):
            print(s_combs[i], " : ", round(prob[i].real * 100, 3), " %")

    def prob_one_bit(self, selection, q_size):
        sbits = 0
        chbits = []
        for i in range(0, len(selection)):
            if (selection[i] == '1'):
                sbits += 1
                chbits.append(i)
        all_combs = list(map("".join, product("01", repeat=q_size)))
        s_combs = list(map("".join, product("01", repeat=sbits)))
        prob = np.zeros(len(s_combs), dtype=complex)
        for i in range(0, len(s_combs)):
            for j in range(0, len(all_combs)):
                Z = True
                for k in range(0, sbits):
                    if all_combs[j][chbits[k]] != s_combs[i][k]:
                        Z = False
                if Z:
                    prob[i] += abs(self.state[j][0]) * abs(self.state[j][0])
        if prob[0] > prob[1]:
            return 0
        return 1

    def prob_return(self, selection, s):
        sbits = 0
        chbits = []
        for i in range(0, len(selection)):
            if (selection[i] == '1'):
                sbits += 1
                chbits.append(i)
        all_combs = list(map("".join, product("01", repeat=s)))
        s_combs = list(map("".join, product("01", repeat=sbits)))
        prob = np.zeros(len(s_combs), dtype=complex)
        for i in range(0, len(s_combs)):
            for j in range(0, len(all_combs)):
                Z = True
                for k in range(0, sbits):
                    if (all_combs[j][chbits[k]] != s_combs[i][k]):
                        Z = False
                if (Z):
                    prob[i] += abs(self.state[j][0]) * abs(self.state[j][0])
        mx = max(prob)
        mx = mx * 0.9
        out = []
        for i in range(len(s_combs)):
            if (prob[i] > mx):
                out.append(s_combs[i])
        return out


