from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar
import numpy as np


def qho(x, psi, E):
    hbar = m = k = 1.0
    return np.asarray([psi[1], 2.0 * m * (k * x ** 2 / 2 - E) * psi[0] / hbar ** 2])

def single_shooting_method(tise, x, psi, dpsi, E):
    objective_func = lambda _ : solve_ivp(tise, x, np.asarray([psi[0], dpsi]), args=(_, )).y[0, -1] - psi[1]
    return root_scalar(objective_func, bracket=E)


x = np.asarray([-10.0, 10.0])
psi = np.asarray([0.0, 0.0])
dpsi = 1.0e-12
E = [0.01, 1.0]
print(E0 := single_shooting_method(qho, x, psi, dpsi, E))