import numpy as np
import math
import random
import scipy as sp

KET_0 = np.array([
    [1],
    [0]], dtype=complex)

KET_1 = np.array([
    [0],
    [1]], dtype=complex)

KET_P = 1/math.sqrt(2) * (KET_0 + KET_1)

KET_M = 1/math.sqrt(2) * (KET_0 - KET_1)

class Qubit:

    def __init__(self):
        self.state = KET_0

    def reset(self):
        self.state = KET_0

    def measure(self):
        prob = abs(self.state[0][0]) ** 2
        random_number = random.uniform(0, 1)
        return random_number < prob

    def basis_measure(self, basis):
        prob = abs((basis.transpose() @ self.state)[0][0]) ** 2
        random_number = random.uniform(0, 1)
        return random_number < prob

    def turn(self, angle, axe):
        self.state = sp.linalg.expm(-1j * angle * axe / 2) @ self.state

    def operator(self, operator):
        self.state = operator @ self.state


