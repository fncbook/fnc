import numpy as np
import scipy
import scipy.sparse as sp
from .chapter10 import diffmat2, diffcheb
from .chapter04 import levenberg


def tensorgrid(x, y):
    """
    tensorgrid(x, y)

    Create a tensor grid for a rectangle from its 1d projections x and y.
    Returns a function to reshape a 2d array to a vector, a function to reshape 
    a vector into a 2d array, a function to evaluate a function on the grid, 
    two arrays to give the grid coordinates, and a boolean array to identify 
    the boundary points.
    """
    m, n = len(x) - 1, len(y) - 1
    vec = lambda U: U.T.flatten()
    unvec = lambda u: np.reshape(u, (n+1, m+1)).T
    mtx = lambda h: np.array([[h(xi ,yj) for yj in y] for xi in x])
    X = mtx(lambda x, y: x)
    Y = mtx(lambda x, y: y)

    # Identify boundary points.
    is_boundary = np.tile(True, (m+1, n+1))
    is_boundary[1:-1, 1:-1] = False
    
    return mtx, X, Y, vec, unvec, is_boundary


def poissonfd(f, g, m, xspan, n, yspan):
    """
    poissonfd(f, g, m, xspan, n, yspan)

    Solve Poisson's equation on a rectangle by finite differences. Function f is the
    forcing function and function g gives the  Dirichlet boundary condition. The rectangle
    is the tensor product of intervals xspan and yspan,  and the discretization uses
    m+1 and n+1 points in the two coordinates.

    Return matrices of the solution values, and the coordinate functions, on the grid.
    """
    # Discretize the domain.
    x, Dx, Dxx = diffmat2(m, xspan)
    y, Dy, Dyy = diffmat2(n, yspan)
    mtx, X, Y, vec, unvec, is_boundary = tensorgrid(x, y)
    N = (m+1) * (n+1)    # total number of unknowns

    # Form the collocated PDE as a linear system.
    Dxx = sp.lil_matrix(Dxx)
    Dyy = sp.lil_matrix(Dyy)
    A = sp.kron(sp.eye(n+1, format="lil"), Dxx) + sp.kron(Dyy, sp.eye(m+1, format="lil"))
    b = vec(mtx(f))

    # Apply Dirichlet condition.
    idx = vec(is_boundary)
    scale = np.max(abs(A[n+1, :]))
    I = sp.eye(N, format="lil")
    A[idx, :] = scale * I[idx, :]         # Dirichet assignment
    X_bd, Y_bd = vec(X)[idx], vec(Y)[idx]
    b[idx] = scale * g(X_bd, Y_bd)    # assigned values

    # Solve the linear sytem and reshape the output.
    u = scipy.sparse.linalg.spsolve(A.tocsr(), b)
    U = unvec(u)
    return U, X, Y


def elliptic(f, g, m, xspan, n, yspan):
    """
    newtonpde(f, g, m, xspan, n, yspan)

    Newton's method with finite differences to solve the PDE f(u,x,y,disc)=0 on the
    rectangle xspan \times yspan, subject to g(x,y)=0 on the boundary. Use m+1
    points in x by n+1 points in y.

    Return matrices of the solution values, and the coordinate functions, on the grid.
    """
    from scipy.sparse.linalg import spsolve
    x, Dx, Dxx = diffcheb(m, xspan)
    y, Dy, Dyy = diffcheb(n, yspan)
    mtx, X, Y, vec, unvec, is_boundary = tensorgrid(x, y)

    # Evaluate the boundary condition at the boundary nodes.
    idx = vec(is_boundary)
    X_bd, Y_bd = vec(X)[idx], vec(Y)[idx]
    g_bd = g(X_bd, Y_bd)

    # Evaluate the PDE+BC residual.
    def residual(u):
        U = unvec(u)
        R = f(X, Y, U, Dx @ U, Dxx @ U, U @ Dy.T, U @ Dyy.T)    # PDE
        r = vec(R)
        r[idx] = u[idx] - g_bd                                  # BC
        return r
    
    # Solve the equation.
    u = levenberg(residual, vec(np.zeros(X.shape)))[-1]
    U = unvec(u)

    def evaluate(xi, eta):
        v = [chebinterp(y, u, eta) for u in U]
        return chebinterp(x, v, xi)
    
    return np.vectorize(evaluate)

def chebinterp(x, v, xi):
    "Evaluate Chebyshev interpolant with nodes x, values v, at point xi"
    n = len(x) - 1
    w = np.ones(n+1)
    w[1::2] = -1    # alternating ±1
    w[[0, n]] *= 0.5
    if xi in x:    # exactly at a node
        return v[np.where(x == xi)[0][0]]
    else:
        terms = w / (xi - x)
        return np.sum(v * terms) / np.sum(terms)
