import numpy as np
from .FNC04 import levenberg
import warnings

def euler(dudt, tspan, u0, n):
    """
    euler(dudt,tspan,u0,n)

    Apply Euler's method to solve the IVP u'=`dudt`(u,t) over the interval `tspan` with
    u(`tspan[1]`)=`u0`, using `n` subintervals/steps. Return vectors of times and solution
    values.
    """
    a, b = tspan
    h = (b - a) / n
    t = np.linspace(a, b, n+1)
    u = np.tile(np.array(u0), (n+1, 1))
    for i in range(n):
        u[i+1] = u[i] + h * dudt(t[i], u[i])

    return t, u.T

def ie2(dudt, tspan, u0, n):
    """
    ie2(dudt,tspan,u0,n)

    Apply the Improved Euler method to solve the vector-valued IVP u'=`dudt`(u,p,t) over the
    interval `tspan` with u(`tspan[1]`)=`u0`, using `n` subintervals/steps. Returns a vector
    of times and a vector of solution values/vectors.
    """
    # Time discretization.
    a, b = tspan
    h = (b - a) / n
    t = np.linspace(a, b, n + 1)

    # Initialize output.
    u = np.tile(np.array(u0), (n+1, 1))

    # Time stepping.
    for i in range(n):
        uhalf = u[i] + h / 2 * dudt(t[i], u[i])
        u[i+1] = u[i] + h * dudt(t[i] + h / 2, uhalf)

    return t, u.T

def rk4(dudt, tspan, u0, n):
    """
    rk4(dudt,tspan,u0,n)

    Apply "the" Runge-Kutta 4th order method to solve the vector-valued IVP u'=`dudt`(u,p,t)
    over the interval `tspan` with u(`tspan[1]`)=`u0`, using `n` subintervals/steps.
    Return a vector of times and a vector of solution values/vectors.
    """
    # Time discretization.
    a, b = tspan
    h = (b - a) / n
    t = np.linspace(a, b, n + 1)

    # Initialize output.
    u = np.tile(np.array(u0), (n+1, 1))

    # Time stepping.
    for i in range(n):
        k1 = h * dudt(t[i], u[i])
        k2 = h * dudt(t[i] + h / 2, u[i] + k1 / 2)
        k3 = h * dudt(t[i] + h / 2, u[i] + k2 / 2)
        k4 = h * dudt(t[i] + h, u[i] + k3)
        u[i+1] = u[i] + (k1 + 2 * (k2 + k3) + k4) / 6

    return t, u.T

def rk23(dudt, tspan, u0, tol):
    """
    rk23(dudt,tspan,u0,tol)

    Apply adaptive embedded RK formula to solve the vector-valued IVP u'=`dudt`(u,p,t)
    over the interval `tspan` with u(`tspan[1]`)=`u0`, with error tolerance `tol`.
    Return a vector of times and a vector of solution values/vectors.
    """
    # Initialize for the first time step.
    t = [tspan[0]]
    u = [np.array(u0)]
    i = 0
    h = 0.5 * tol ** (1 / 3)
    s1 = dudt(t, u[0])

    # Time stepping.
    while t[i] < tspan[-1]:
        # Detect underflow of the step size.
        if t[i] + h == t[i]:
            warnings.warn(f"Stepsize too small near t={t[i]}")
            break  # quit time stepping loop

        # New RK stages.
        s2 = dudt(t[i] + h / 2, u[i] + (h / 2) * s1)
        s3 = dudt(t[i] + 3 * h / 4, u[i] + (3 * h / 4) * s2)
        unew2 = u[i] + h * (2 * s1 + 3 * s2 + 4 * s3) / 9  # 2rd order solution
        s4 = dudt(t[i] + h, unew2)
        err = h * (-5 * s1 / 72 + s2 / 12 + s3 / 9 - s4 / 8)  # 2nd/3rd order difference
        E = np.linalg.norm(err, np.inf)  # error estimate
        maxerr = tol * (1 + np.linalg.norm(u[i], np.inf))  # relative/absolute blend

        # Accept the proposed step?
        if E < maxerr:  # yes
            t.append(t[i] + h)
            u = np.vstack([u, unew2])
            i += 1
            s1 = s4  # use FSAL property

        # Adjust step size.
        q = 0.8 * (maxerr / E) ** (1 / 3)  # conservative optimal step factor
        q = min(q, 4)  # limit stepsize growth
        h = min(q * h, tspan[-1] - t[i])  # don't step past the end

    # Convert outputs to arrays
    return np.array(t), u.T

def ab4(dudt, tspan, u0, n):
    """
    ab4(dudt,tspan,u0,n)

    Apply the Adams-Bashforth 4th order method to solve the vector-valued IVP u'=`dudt`(u,p,t)
    over the interval `tspan` with u(`tspan[1]`)=`u0`, using `n` subintervals/steps.
    """
    # Time discretization.
    a, b = tspan
    h = (b - a) / n
    t = np.linspace(a, b, n + 1)

    # Constants in the AB4 method.
    k = 4
    sigma = np.array([55, -59, 37, -9]) / 24

    # Find starting values by RK4.
    ts, us = rk4(dudt, [a, a + (k - 1) * h], u0, k - 1)
    u = np.tile(np.array(u0), (n+1, 1))
    u[:k] = us[:k].T

    # Compute history of u' values, from newest to oldest.
    f = np.array([dudt(t[k-j-2], u[k-j-2]) for j in range(k)])

    # Time stepping.
    for i in range(k-1, n):
        f = np.vstack([dudt(t[i], u[i]), f[:-1]])  # new value of du/dt
        u[i+1] = u[i] + h * np.dot(sigma, f)  # advance one step

    return t, u.T

def am2(dudt, tspan, u0, n):
    """
    am2(dudt,tspan,u0,n)

    Apply the Adams-Moulton 2nd order method to solve the vector-valued IVP u'=`dudt`(u,p,t)
    over the interval `tspan` with u(`tspan[1]`)=`u0`, using `n` subintervals/steps.
    """
    # Time discretization.
    a, b = tspan
    h = (b - a) / n
    t = np.linspace(a, b, n + 1)

    # Initialize output.
    u = np.tile(np.array(u0), (n+1, 1))

    # Time stepping.
    for i in range(n):
        # Data that does not depend on the new value.
        known = u[i] + h / 2 * dudt(t[i], u[i])
        # Find a root for the new value.
        F = lambda z: z - h / 2 * dudt(t[i+1], z) - known
        unew = levenberg(F, known)
        u[i+1] = unew[:, -1]

    return t, u.T
