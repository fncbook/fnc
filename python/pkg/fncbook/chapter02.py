import numpy as np

def forwardsub(L,b):
    """
     forwardsub(L,b)

    Solve the lower-triangular linear system with matrix L and right-hand side
    vector b.
    """
    n = len(b)
    x = np.zeros(n)
    for i in range(n):
        s = L[i,:i] @ x[:i]
        x[i] = ( b[i] - s ) / L[i, i]
    return x


def backsub(U,b):
    """
    backsub(U,b)

    Solve the upper-triangular linear system with matrix U and right-hand side
    vector b.
    """
    n = len(b)
    x = np.zeros(n)
    for i in range(n-1, -1, -1):
        s = U[i, i+1:] @ x[i+1:]
        x[i] = ( b[i] - s ) / U[i, i]
    return x

def lufact(A):
    """
    lufact(A)

    Compute the LU factorization of square matrix A, returning the
    factors.
    """
    n = A.shape[0]     # detect the dimensions from the input
    L = np.eye(n)      # ones on main diagonal, np.zeros elsewhere
    U = np.zeros((n, n))
    A_k = np.copy(A)   # make a working np.copy 

    # Reduction by np.outer products
    for k in range(n-1):
        U[k, :] = A_k[k, :]
        L[:, k] = A_k[:, k] / U[k,k]
        A_k -= np.outer(L[:,k], U[k,:])
    U[n-1, n-1] = A_k[n-1, n-1]
    return L, U

def plufact(A):
    """
        plufact(A)

    Compute the PLU factorization of square matrix A, returning the
    triangular factors and a row permutation vector.
    """
    n = A.shape[0]
    L = np.zeros((n, n))
    U = np.zeros((n, n))
    p = np.zeros(n, dtype=int)
    A_k = np.copy(A)

    # Reduction by np.outer products
    for k in range(n):
        p[k] = np.argmax(abs(A_k[:, k]))
        U[k, :] = A_k[p[k], :]
        L[:, k] = A_k[:, k] / U[k, k]
        if k < n-1:
            A_k -= np.outer(L[:, k], U[k, :])
    return L[p, :], U, p