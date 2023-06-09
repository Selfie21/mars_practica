#!/usr/bin/python
from math import sqrt, acos

from cagd.vec import vec2, vec3


# solves the system of linear equations Ax = res
# where A is a tridiagonal matrix with diag2 representing the main diagonal
# diag1 and diag3 represent the lower and upper diagonal respectively
# all four parameters are vectors of size n
# the first element of diag1 and the last element of diag3 are ignored
# therefore diag1[i], diag2[i] and diag3[i] are located on the same row of A
# a = diag1 b = diag2 c=diag3
def solve_tridiagonal_equation(diag1, diag2, diag3, res):
    assert (len(diag1) == len(diag2) == len(diag3) == len(res))
    v = [0] * len(diag2)
    x = [0] * len(diag2)
    y = [0] * len(diag2)
    z = [0] * len(diag2)

    dim = len(diag2) - 1
    assert (v[0] == y[-1] == diag1[0] == diag3[dim])
    for k in range(0, dim):
        z[k] = 1 / (diag2[k] - (diag1[k] * v[k]))
        v[k + 1] = z[k] * diag3[k]
        y[k] = z[k] * (res[k] - (diag1[k] * y[k - 1]))
    z[dim] = 1 / ((diag2[dim]) - (diag1[dim] * v[dim]))
    y[dim] = z[dim] * (res[dim] - (diag1[dim] * y[dim - 1]))
    x[dim] = y[dim]
    for k in reversed(range(dim)):
        x[k] = y[k] - (v[k + 1] * x[k + 1])
    solution = x
    return solution


# solves the system of linear equations Ax = res
# where A is an almost tridiagonal matrix with diag2 representing the main diagonal
# diag1 and diag3 represent the lower and upper diagonal respectively
# all four parameters are vectors of size n
# the first element of diag1 and the last element of diag3 represent the top right and bottom left elements of A
# diag1[i], diag2[i] and diag3[i] are located on the same row of A
def solve_almost_tridiagonal_equation(diag1, diag2, diag3, res):
    assert (len(diag1) == len(diag2) == len(diag3) == len(res))
    dim = len(diag1)
    v = [0] * dim
    z = [0] * dim
    y = [0] * dim
    s = [0] * dim
    t = [0] * dim
    w = [0] * dim
    x = [0] * dim
    s[-1] = 1
    dim = dim - 1  # Adjusting dim to fit the max index
    for i in range(0, dim):  # iterating from second index to second last index
        z[i] = 1 / (diag2[i] + diag1[i] * v[i - 1])
        v[i] = -z[i] * diag3[i]
        y[i] = z[i] * (res[i] - diag1[i] * y[i - 1])
        s[i] = -diag1[i] * s[i - 1] * z[i]
    t[dim] = 1
    for k in reversed(range(0, dim)):
        t[k] = v[k] * t[k + 1] + s[k]
        w[k] = v[k] * w[k + 1] + y[k]
    x[dim] = (res[dim] - diag3[dim] * w[0] - diag1[dim] * w[dim - 1]) / (
            diag3[dim] * t[0] + diag1[dim] * t[dim - 1] + diag2[dim])
    for k in reversed(range(0, dim)):
        x[k] = t[k] * x[dim] + w[k]
    return x


# Algebraic
def distance(a, b):
    if isinstance(a, vec2):
        diff = vec2(b.x - a.x, b.y - a.y)
        return euclidean_norm(diff)
    elif isinstance(a, vec3):
        diff = vec3(b.x - a.x, b.y - a.y, b.z - a.z)
        return euclidean_norm(diff)


def euclidean_norm(a):
    if isinstance(a, vec2):
        return sqrt((a.x ** 2) + (a.y ** 2))
    elif isinstance(a, vec3):
        return sqrt((a.x ** 2) + (a.y ** 2) + (a.z ** 2))


def angle(a, b):
    vector_product = a.x * b.x + a.y * b.y
    cosx = vector_product / (euclidean_norm(a) * euclidean_norm(b))
    return acos(cosx)


# matrix is 2x2 in list form for meaning each row gets iterated one after another
def determinant_m2(matrix):
    return (matrix[0] * matrix[3]) - (matrix[1] * matrix[2])
