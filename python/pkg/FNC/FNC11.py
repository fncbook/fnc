import numpy as np

def diffper(n, xspan):
    """
    diffper(n,xspan)

    Construct 2nd-order differentiation matrices for functions with periodic end
    conditions, using `n` unique nodes in the interval `xspan`. Return a vector of
    nodes and the  matrices for the first and second derivatives.
    """
    a, b = xspan
    h = (b - a) / n
    x = a + h * np.arange(n)  # nodes, omitting the repeated data

    # Construct Dx by diagonals, then correct the corners.
    dp = 0.5 / h * np.ones(n - 1)  # superdiagonal
    dm = -0.5 / h * np.ones(n - 1)  # subdiagonal
    Dx = np.diag(dm, -1) + np.diag(dp, 1)
    Dx[0, -1] = -1 / (2 * h)
    Dx[-1, 0] = 1 / (2 * h)

    # Construct Dxx by diagonals, then correct the corners.
    d0 = -2 / h**2 * np.ones(n)  # main diagonal
    dp = np.ones(n - 1) / h**2  # superdiagonal and subdiagonal
    Dxx = np.diag(d0) + np.diag(dp, -1) + np.diag(dp, 1)
    Dxx[0, -1] = 1 / (h**2)
    Dxx[-1, 0] = 1 / (h**2)

    return x, Dx, Dxx