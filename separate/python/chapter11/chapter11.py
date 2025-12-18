import numpy as np
from .chapter04 import levenberg
from .chapter10 import diffcheb
from scipy.integrate import solve_ivp

def diffper(n, xspan):
    """
    diffper(n, xspan)

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


def parabolic(phi, xspan, m, ga, gb, tspan, init):
    """
        parabolic(phi, xspan, m, ga, gb, tspan, init)

    Solve a parabolic PDE by the method of lines. The PDE is 
    ∂u/∂t = phi(t,x,u,∂u/∂x,∂^2u/∂x^2), xspan gives the space 
    domain, m gives the degree of a Chebyshev spectral discretization, ga and gb are functions of (u,∂u/∂x) at the domain ends that should be made zero, tspan is the time domain, and init is a function of x that gives the initial condition. Returns a vector x and a function of t that gives the semidiscrete solution at x. 
    """   
    x, Dx, Dxx = diffcheb(m, xspan)
    int = range(1, m)    # indexes of interior nodes

    def extend(v):
        def objective(ubc):
            u0, um = ubc
            ux = Dx @ np.hstack([u0, v, um])
            return np.array([ga(u0, ux[1]), gb(um, ux[-1])])
        ubc = levenberg(objective, np.array([0, 0]))[-1]
        return np.hstack([ubc[0], v, ubc[-1]])

    def ode(t, v):
        u = extend(v)
        ux, uxx = Dx @ u, Dxx @ u
        return phi(t, x[int], u[int], ux[int], uxx[int])

    u0 = init(x[int])
    sol = solve_ivp(ode, tspan, u0, method="BDF", dense_output=True)

    return x, lambda t: extend(sol.sol(t))
