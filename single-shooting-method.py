from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar
import numpy as np

def func(x, y):
    return np.array([y[1], 3.0 * y[0] ** 2 / 2.0])

def single_shooting_method(func, x, y, p):
    objective_func = lambda _ : solve_ivp(func, x, np.array([y[0], _])).y[0, -1] - y[1]
    p_opt = root_scalar(objective_func, x0=p[0], x1=p[1])
    return p_opt


x = np.array([0.0, 1.0])
y = np.array([4.0, 1.0])
p = [-10, -2]

res = single_shooting_method(func, x, y, p)
print(res)