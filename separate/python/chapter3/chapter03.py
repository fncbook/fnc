import scipy
import numpy as np
from .chapter02 import forwardsub, backsub

def lsnormal(A, b):
    """
    lsnormal(A, b)
    
    Solve a linear least squares problem by the normal equations. Returns the
    minimizer of ||b-Ax||.
    """
    N = A.T @ A
    z = A.T @ b
    R = scipy.linalg.cholesky(N)
    w = forwardsub(R.T, z)                   # solve R'z=c
    x = backsub(R, w)                        # solve Rx=z
    return x

def lsqrfact(A, b):
    """
    lsqrfact(A, b)
    
    Solve a linear least squares problem by QR factorization. Returns the
    minimizer of ||b-Ax||.
    """
    Q, R = np.linalg.qr(A)
    c = Q.T @ b
    x = backsub(R, c)
    return x

def qrfact(A):
    """
        qrfact(A)

    QR factorization by Householder reflections. Returns Q and R.
    """
    m, n = A.shape
    Qt = np.eye(m)
    R = np.copy(A)
    for k in range(n):
        z = R[k:, k]
        w = np.hstack((-np.sign(z[0]) * np.linalg.norm(z) - z[0], -z[1:]))
        nrmw = np.linalg.norm(w)
        if nrmw < np.finfo(float).eps: continue    # skip this iteration
        v = w / nrmw
        # Apply the reflection to each relevant column of R and Q
        for j in range(k, n):
            R[k:, j] -= 2 * np.dot(v, R[k:, j]) * v
        for j in range(m):
            Qt[k:, j] -= 2 * np.dot(v, Qt[k:, j]) * v 
    return Qt.T, np.triu(R)
