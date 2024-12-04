from scipy.linalg import cholesky
from numpy.linalg import qr, norm
from numpy import *
eps = finfo(float).eps
from .FNC02 import forwardsub, backsub

def lsnormal(A,b):
    """
    lsnormal(A,b)
    
    Solve a linear least squares problem by the normal equations. Returns the
    minimizer of ||b-Ax||.
    """
    N = A.T @ A
    z = A.T @ b
    R = cholesky(N)
    w = forwardsub(R.T, z)                   # solve R'z=c
    x = backsub(R, w)                        # solve Rx=z
    return x

def lsqrfact(A,b):
    """
    lsqrfact(A,b)
    
    Solve a linear least squares problem by QR factorization. Returns the
    minimizer of ||b-Ax||.
    """
    Q, R = qr(A)
    c = Q.T @ b
    x = backsub(R, c)
    return x

def qrfact(A):
    """
        qrfact(A)

    QR factorization by Householder reflections. Returns Q and R.
    """
    m, n = A.shape
    Qt = eye(m)
    R = copy(A)
    for k in range(n):
        z = R[k:, k]
        w = hstack((-sign(z[0]) * norm(z) - z[0], -z[1:]))
        nrmw = norm(w)
        if nrmw < eps: continue    # skip this iteration
        v = w / nrmw
        # Apply the reflection to each relevant column of R and Q
        for j in range(k, n):
            R[k:, j] -= 2 * dot(v, R[k:, j]) * v
        for j in range(m):
            Qt[k:, j] -= 2 * dot(v, Qt[k:, j]) * v 
    return Qt.T, triu(R)
