import numpy as np

def calc_transform_matrix_XY(fi):
    """
    param fi: angle of rotation around Z axes
    return: transform matrix
    """
    a = np.array([  [np.cos(fi), -np.sin(fi), 0.0],
                    [np.sin(fi), np.cos(fi), 0.0],
                    [0.0, 0.0, 1.0]])
    return a

def new_coord(x, matr):
    """
    point coordinate in new axis
    """
    #print(f"x={x}")
    return np.dot(matr, x)

def line_point(x0, vec, t):
    res = []
    for i in range(len(x0)):
        res.append(x0[i] + vec[i] * t)
    return res

def vec(p1, p2):
    res = []
    for i in range(len(p1)):
        res.append(p2[i] - p1[i])
    return res

def distance(p1, p2):
    res = 0.0
    for i in range(len(p1)):
        res += (p2[i] - p1[i])**2
    return res**0.5

def kvadratura(f, a, b):
    '''
    8-point gauss kvadratura
    '''
    ksi = [0.96028986, 0.79666648, 0.52553242, 0.183434464,
           -0.183434464, -0.52553242, -0.79666648, -0.96028986]
    ita = [0.10122854, 0.22238104, 0.31370664, 0.36268378,
           0.36268378, 0.31370664, 0.22238104, 0.10122854]
    res = 0
    for i in range(8):
        xi = (b + a) / 2 + (b - a) / 2 * ksi[i]
        res += ita[i] * f(xi)
    res *= (b - a) / 2
    return res

def kvadratura_N(f, a, b, n):
    '''
    8-point gauss kvadratura
    '''
    ksi = [0.96028986, 0.79666648, 0.52553242, 0.183434464,
           -0.183434464, -0.52553242, -0.79666648, -0.96028986]
    ita = [0.10122854, 0.22238104, 0.31370664, 0.36268378,
           0.36268378, 0.31370664, 0.22238104, 0.10122854]
    res = 0
    for j in range(n):
        aj = (b - a) / n * j + a
        bj = (b - a) / n * (j + 1) + a
        resj = 0
        for i in range(8):
            xi = (bj + aj) / 2 + (bj - aj) / 2 * ksi[i]
            resj += ita[i] * f(xi)
        resj *= (bj - aj) / 2
        #print(f"resj={resj}")
        res += resj
    return res


def get_value(beta, ksi):
    '''
    linear interpolation at array
    '''
    bt = beta[0][0]
    for i in range(len(beta) - 1):
        if (beta[i][0] < ksi) and (beta[i + 1][0] >= ksi):
            bt = (beta[i + 1][1] - beta[i][1]) / (beta[i + 1][0] - beta[i][0]) * (ksi - beta[i][0]) + beta[i][1]
    return bt
