import numpy as np
from numpy.linalg import solve

def hatfun(t, k):
    """
    hatfun(t, k)

    Returns a piecewise linear "hat" function,  where t is a vector of
    n+1 interpolation nodes and k is an integer in 0:n giving the index of the node
    where the hat function equals one.
    """
    n = len(t) - 1

    def evaluate(x):
        H = np.zeros(np.array(x).shape)
        for (j, xj) in enumerate(x):
            if (k > 0) and (t[k-1] <= xj) and (xj <= t[k]):
                H[j] = (xj - t[k-1]) / (t[k] - t[k-1])
            elif (k < n) and (t[k] <= xj) and (xj <= t[k+1]):
                H[j] = (t[k+1] - xj) / (t[k+1] - t[k])
        return H
    return evaluate

def plinterp(t, y):
    """
    plinterp(t, y)

    Create a piecewise linear interpolating function for data values in y given at nodes
    in t.
    """
    n = len(t) - 1
    H = [hatfun(t, k) for k in range(n+1)]
    def evaluate(x):
        f = 0
        for k in range(n+1):
            f += y[k] * H[k](x)
        return f
    return evaluate

def spinterp(t, y):
    """
    spinterp(t, y)

    Create a cubic not-a-knot spline interpolating function for data values in y given at nodes in t.
    """
    n = len(t) - 1
    h = [t[i + 1] - t[i] for i in range(n)]

    # Preliminary definitions.
    Z = np.zeros([n, n])
    I = np.eye(n)
    E = I[: n - 1, :]
    J = np.eye(n) + np.diag(-np.ones(n - 1), 1)
    H = np.diag(h)

    # Left endpoint interpolation:
    AL = np.hstack([I, Z, Z, Z])
    vL = y[:-1]

    # Right endpoint interpolation:
    AR = np.hstack([I, H, H**2, H**3])
    vR = y[1:]

    # Continuity of first derivative:
    A1 = E @ np.hstack([Z, J, 2 * H, 3 * H**2])
    v1 = np.zeros(n - 1)

    # Continuity of second derivative:
    A2 = E @ np.hstack([Z, Z, J, 3 * H])
    v2 = np.zeros(n - 1)

    # Not-a-knot conditions:
    nakL = np.hstack([np.zeros(3 * n), np.hstack([1, -1, np.zeros(n - 2)])])
    nakR = np.hstack([np.zeros(3 * n), np.hstack([np.zeros(n - 2), 1, -1])])

    # Assemble and solve the full system.
    A = np.vstack([AL, AR, A1, A2, nakL, nakR])
    v = np.hstack([vL, vR, v1, v2, 0, 0])
    z = solve(A, v)

    # Break the coefficients into separate vectors.
    rows = np.arange(n)
    a = z[rows]
    b = z[n + rows]
    c = z[2 * n + rows]
    d = z[3 * n + rows]
    S = [np.poly1d([d[k], c[k], b[k], a[k]]) for k in range(n)]

    # This function evaluates the spline when called with a value for x.
    def evaluate(x):
        f = np.zeros(x.shape)
        for k in range(n):
            # Evaluate this piece's cubic at the points inside it.
            index = (x >= t[k]) & (x <= t[k + 1])
            f[index] = S[k](x[index] - t[k])
        return f

    return evaluate

def fdweights(t, m):
    """
    fdweights(t, m)

    Return weights for the mth derivative of a function at zero using values at the
    nodes in vector t.
    """
    # This is a compact implementation, not an efficient one.

    def weight(t, m, r, k):
        # Recursion for one weight.
        # Input:
        #   t   nodes (vector)
        #   m   order of derivative sought
        #   r   number of nodes to use from t (<= length(t))
        #   k   index of node whose weight is found

        if (m < 0) or (m > r):  # undefined coeffs must be zero
            c = 0
        elif (m == 0) and (r == 0):  # base case of one-point interpolation
            c = 1
        else:  # generic recursion
            if k < r:
                denom = t[r] - t[k]
                c = (t[r] * weight(t, m, r-1, k) - m * weight(t, m-1, r-1, k)) / denom
            else:
                beta = np.prod(t[r-1] - t[:r-1]) / np.prod(t[r] - t[:r])
                c = beta * (m * weight(t, m-1, r-1, r-1) - t[r-1] * weight(t, m, r-1, r-1))
        return c

    r = len(t) - 1
    w = np.zeros(t.shape)
    return np.array([ weight(t, m, r, k) for k in range(r+1) ])

def trapezoid(f, a, b, n):
    """
    trapezoid(f, a, b, n)

    Apply the trapezoid integration formula for integrand f over interval [a,b], broken up into n equal pieces. Returns estimate, vector of nodes, and vector of integrand values at the nodes.
    """
    h = (b - a) / n
    t = np.linspace(a, b, n + 1)
    y = f(t)
    T = h * (np.sum(y[1:-1]) + 0.5 * (y[0] + y[-1]))
    return T, t, y

def intadapt(f, a, b, tol):
    """
    intadapt(f, a, b, tol)

    Do adaptive integration to estimate the integral of f over [a,b] to desired
    error tolerance tol. Returns estimate and a vector of evaluation nodes used.
    """

    # Use error estimation and recursive bisection.
    def do_integral(a, fa, b, fb, m, fm, tol):
        # These are the two new nodes and their f-values.
        xl = (a + m) / 2
        fl = f(xl)
        xr = (m + b) / 2
        fr = f(xr)
        t = np.array([a, xl, m, xr, b])  # all 5 nodes at this level

        # Compute the trapezoid values iteratively.
        h = b - a
        T = np.zeros(3)
        T[0] = h * (fa + fb) / 2
        T[1] = T[0] / 2 + (h / 2) * fm
        T[2] = T[1] / 2 + (h / 4) * (fl + fr)

        S = (4 * T[1:] - T[:-1]) / 3  # Simpson values
        E = (S[1] - S[0]) / 15  # error estimate

        if abs(E) < tol * (1 + abs(S[1])):  # acceptable error?
            Q = S[1]  # yes--done
        else:
            # Error is too large--bisect and recurse.
            QL, tL = do_integral(a, fa, m, fm, xl, fl, tol)
            QR, tR = do_integral(m, fm, b, fb, xr, fr, tol)
            Q = QL + QR
            t = np.hstack([tL, tR[1:]])  # merge the nodes w/o duplicate
        return Q, t

    m = (b + a) / 2
    Q, t = do_integral(a, f(a), b, f(b), m, f(m), tol)
    return Q, t
