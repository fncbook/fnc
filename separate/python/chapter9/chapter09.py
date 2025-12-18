import warnings
import numpy as np
from numpy.linalg import eig
from .chapter05 import intadapt


def polyinterp(t, y):
    """
    polyinterp(t, y)

    Return a callable polynomial interpolant through the points in vectors t, y. Uses
    the barycentric interpolation formula.
    """
    n = len(t) - 1
    C = (t[-1] - t[0]) / 4  # scaling factor to ensure stability
    tc = t / C

    # Adding one node at a time, compute inverses of the weights.
    omega = np.ones(n + 1)
    for m in range(n):
        d = tc[: m + 1] - tc[m + 1]  # vector of node differences
        omega[: m + 1] = omega[: m + 1] * d  # update previous
        omega[m + 1] = np.prod(-d)  # compute the new one
    w = 1 / omega  # go from inverses to weights

    def p(x):
        # Compute interpolant.
        z = np.where(x == t)[0]
        if len(z) > 0:  # avoid dividing by zero
            # Apply L'Hopital's Rule exactly.
            f = y[z[0]]
        else:
            terms = w / (x - t)
            f = np.sum(y * terms) / np.sum(terms)
        return f

    return np.vectorize(p)


def triginterp(t, y):
    """
        triginterp(t, y)

    Return trigonometric interpolant for points defined by vectors t and y.
    """
    N = len(t)

    def trigcardinal(x):
        if x == 0:
            tau = 1.0
        elif np.mod(N, 2) == 1:  # odd
            tau = np.sin(N * np.pi * x / 2) / (N * np.sin(np.pi * x / 2))
        else:  # even
            tau = np.sin(N * np.pi * x / 2) / (N * np.tan(np.pi * x / 2))
        return tau

    def p(x):
        return np.sum([y[k] * trigcardinal(x - t[k]) for k in range(N)])

    return np.vectorize(p)


def ccint(f, n):
    """
    ccint(f, n)

    Perform Clenshaw-Curtis integration for the function f on n+1 nodes in [-1,1]. Return
    integral and a vector of the nodes used. Note: n must be even.
    """
    # Find Chebyshev extreme nodes.
    theta = np.linspace(0, np.pi, n + 1)
    x = -np.cos(theta)

    # Compute the C-C weights.
    c = np.zeros(n + 1)
    c[[0, n]] = 1 / (n**2 - 1)
    v = np.ones(n - 1)
    for k in range(1, int(n / 2)):
        v -= 2 * np.cos(2 * k * theta[1:-1]) / (4 * k**2 - 1)
    v -= np.cos(n * theta[1:-1]) / (n**2 - 1)
    c[1:-1] = 2 * v / n

    # Evaluate integrand and integral.
    I = np.dot(c, f(x))  # use vector inner product
    return I, x


def glint(f, n):
    """
    glint(f, n)

    Perform Gauss-Legendre integration for the function f on n nodes in (-1,1). Return
    integral and a vector of the nodes used.
    """
    # Nodes and weights are found via a tridiagonal eigenvalue problem.
    beta = 0.5 / np.sqrt(1 - (2.0 * np.arange(1, n)) ** (-2))
    T = np.diag(beta, 1) + np.diag(beta, -1)
    ev, V = eig(T)
    ev = np.real_if_close(ev)
    p = np.argsort(ev)
    x = ev[p]  # nodes
    c = 2 * V[0, p] ** 2  # weights

    # Evaluate the integrand and compute the integral.
    I = np.dot(c, f(x))  # vector inner product
    return I, x

def intinf(f, tol):
    """
    intinf(f, tol)

    Perform doubly-exponential integration of function f over (-Inf,Inf), using
    error tolerance tol. Return integral and a vector of the nodes used.
    """
    xi = lambda t: np.sinh(np.sinh(t))
    dxi_dt = lambda t: np.cosh(t) * np.cosh(np.sinh(t))
    g = lambda t: f(xi(t)) * dxi_dt(t)
    M = 3
    while (abs(g(-M)) > tol/100) or (abs(g(M)) > tol/100):
        M += 0.5
        if np.isinf(xi(M)):
            warnings.warn("Function may not decay fast enough.")
            M -= 0.5
            break

    I, t = intadapt(g,-M,M,tol)
    x = xi(t)
    return I, x

def intsing(f, tol):
    """
    intsing(f, tol)

    Adaptively integrate function f over (0,1), where f may be 
    singular at zero, with error tolerance tol. Returns the
    integral estimate and a vector of the nodes used.
    """
    xi = lambda t: 2 / (1 + np.exp( 2*np.sinh(t) ))
    dxi_dt = lambda t: np.cosh(t) / np.cosh(np.sinh(t))**2
    g = lambda t: f(xi(t)) * dxi_dt(t)
    # Find where to truncate the integration interval.
    M = 3
    while abs(g(M)) > tol/100:
        M += 0.5
        if xi(M) == 0:
            warnings.warn("Function may grow too rapidly.")
            M -= 0.5
            break

    I, t = intadapt(g, 0, M, tol)
    x = xi(t)
    return I, x
