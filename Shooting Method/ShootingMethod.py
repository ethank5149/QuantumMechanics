import matplotlib.pyplot as plt
import numpy as np


hbar = 1
m = 1
k = 1
w = np.sqrt(k/m)
L = 25

x_range = [-L, L]
y_0 = [0.0, 1]
y_f = [0.0, -1]
num_points = 1000
max_iterations = 1000


def V(x):
    return 0.5*m*w**2*x**2


def f(x, y, E):
    return np.array([y[1], 2*m*(V(x)-E)*y[0]/hbar**2])


def rk4_(E, f=f, x_range=x_range, y_0=y_0, num_points=num_points):
    x = np.linspace(x_range[0], x_range[1], num_points)
    y = np.zeros((2, num_points))
    y[:, 0] = y_0
    dx = (x_range[1]-x_range[0])/float(num_points)
    for i in range(1, num_points):
        k1 = dx*f(x[i-1], y[:, i-1], E)
        k2 = dx*f(x[i-1], y[:, i-1] + 0.5 * k1, E)
        k3 = dx*f(x[i-1], y[:, i-1] + 0.5 * k2, E)
        k4 = dx*f(x[i-1], y[:, i-1] + k3, E)
        y[:, i] = y[:, i-1] + (k1 + 2.0*(k2 + k3) + k4) / 6.0
    return y[0, -1]  # x, y


def rk4(E, f=f, x_range=x_range, y_0=y_0, num_points=num_points):
    x = np.linspace(x_range[0], x_range[1], num_points)
    y = np.zeros((2, num_points))
    y[:, 0] = y_0
    dx = (x_range[1]-x_range[0])/float(num_points)
    for i in range(1, num_points):
        k1 = dx*f(x[i-1], y[:, i-1], E)
        k2 = dx*f(x[i-1], y[:, i-1] + 0.5 * k1, E)
        k3 = dx*f(x[i-1], y[:, i-1] + 0.5 * k2, E)
        k4 = dx*f(x[i-1], y[:, i-1] + k3, E)
        y[:, i] = y[:, i-1] + (k1 + 2.0*(k2 + k3) + k4) / 6.0
    return y[0, :]


def bisection(bracket, f=rk4_, max_iterations=max_iterations):
    a, b = bracket[0], bracket[1]
    if f(a)*f(b) > 0:
        print("Unique Solution Isn't Bracketed")
        return None
    if a > b:
        a_ = a
        a = b
        b = a_
    for i in range(max_iterations):
        c = (a+b)/2
        f_c = f(c)
        if f_c > 0:
            b = c
        elif f_c < 0:
            a = c
        else:
            print("Exact Solution Found")
            return c
    return (a+b)/2


def secant(bracket, f=rk4_, max_iterations=max_iterations):
    em1, em2 = bracket[0], bracket[1]
    e = (em2*f(em1)-em1*f(em2))/(f(em1)-f(em2))
    for i in range(max_iterations):
        if abs(em1-em2) < 0.00000001:
            return e
        else:
            em2 = em1
            em1 = e
            e = (em2*f(em1)-em1*f(em2))/(f(em1)-f(em2))
    return e


def brent(bracket, f=rk4_, max_iterations=max_iterations):
    a, b = bracket[0], bracket[1]
    f_a, f_b = f(a), f(b)
    if f_a*f_b >= 0:
        print("the root is not bracketed")
        return None
    if abs(f_a) < abs(f_b):
        a_ = a
        a = b
        b = a_
    c = a
    mflag = True
    s = b  # To prevent error
    while f(b) != 0 and f(s) != 0 and abs(b-a) > 0.000001:
        f_a, f_b, f_c = f(a), f(b), f(c)
        if f_a != f_c and f_b != f_c:
            s = (a*f_b*f_c)/((f_a-f_b)*(f_a-f_c)) + (b*f_a*f_c) / \
                ((f_b-f_a)*(f_b-f_c)) + \
                (c*f_a*f_b)/((f_c-f_a)*(f_c-f_b))
        else:
            s = b-f_b*(b-a)/(f_b-f_a)
        cond1 = not ((3*a+b)/4 < s < b)
        cond2 = mflag and abs(s-b) >= abs(b-c)/2
        cond3 = not mflag and abs(s-b) >= abs(c-d)/2
        cond4 = mflag and abs(b-c) < 0.0001
        cond5 = not mflag and abs(c-d) < 0.0001
        if cond1 or cond2 or cond3 or cond4 or cond5:
            s = (a+b)/2
            mflag = True
        else:
            mflag = False
        f_s = f(s)
        d = c
        c = b
        if f_a*f_s < 0:
            b = s
        else:
            a = s
        if abs(f(a)) < abs(f(b)):
            a_ = a
            a = b
            b = a_
    return b


def shooting(param_range):
    return brent(param_range)


def trueEinfsqr(n):
    return 0.5*(np.pi*n/(2*L*hbar))**2


def trueEqho(n):
    return (1+0.5*n)*hbar*w


print(trueEinfsqr(1), trueEinfsqr(2), trueEinfsqr(3), trueEinfsqr(4))
#print(trueEqho(0), trueEqho(1), trueEqho(2), trueEqho(3))
e = shooting((0.0, 0.1))
print(e)

# x = np.linspace(x_range[0], x_range[1], num_points)
# y = rk4(0.5*(np.pi/(2*L))**2-0.01+0.06)
# plt.plot(x, y, 'k-')
# plt.show()
