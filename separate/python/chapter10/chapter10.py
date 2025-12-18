import numpy as np
from scipy.integrate import solve_ivp
from .chapter04 import levenberg

def shoot(phi, a, b, ga, gb, init):
    """
    shoot(phi, a, b, ga, gb, init)

    Use the shooting method to solve a two-point boundary value problem. 
    The ODE is u'' = phi(x, u, u') for x in (a,b). The functions 
    ga(u(a), u'(a)) and gb(u(b), u'(b)) specify the boundary conditions. 
    The value init is an initial guess for [u(a), u'(a)].

    Return vectors for the nodes, the values of u, and the values of u'.
    """

    # Tolerances for IVP solver and rootfinder.
    tol = 1e-5
    # To be solved by the IVP solver
    def shootivp(x, y):
        return np.array([y[1], phi(x, y[0], y[1])])

    # Evaluate the difference between computed and target values at x=b.
    def objective(s):
        nonlocal x, y  # change these values in outer scope

        x = np.linspace(a, b, 400)  # make decent plots on return
        sol = solve_ivp(shootivp, [a, b], s, atol=tol/10, rtol=tol/10, t_eval=x)
        x = sol.t
        y = sol.y
        residual = np.array([ga(y[0, 0], y[1, 0]), gb(y[0, -1], y[0, -1])])
        return residual

    # Find the unknown quantity at x=a by rootfinding.
    x, y = [], []    # the values will be overwritten
    s = levenberg(objective, init, tol)

    # Don't need to solve the IVP again. It was done within the
    # objective function already.
    u = y[0]        # solution
    du_dx = y[1]    # derivative

    return x, u, du_dx


def diffmat2(n, xspan):
    """
    diffmat2(n, xspan)

    Compute 2nd-order-accurate differentiation matrices on n+1 points in the
    interval xspan. Return a vector of nodes, and the matrices for the first
    and second derivatives.
    """
    a, b = xspan
    h = (b - a) / n
    x = np.linspace(a, b, n + 1)  # nodes

    # Define most of Dx by its diagonals.
    dp = 0.5 / h * np.ones(n)  # superdiagonal
    dm = -0.5 / h * np.ones(n)  # subdiagonal
    Dx = np.diag(dm, -1) + np.diag(dp, 1)

    # Fix first and last rows.
    Dx[0, :3] = np.array([-1.5, 2, -0.5]) / h
    Dx[-1, -3:] = np.array([0.5, -2, 1.5]) / h

    # Define most of Dxx by its diagonals.
    d0 = -2 / h**2 * np.ones(n + 1)  # main diagonal
    dp = np.ones(n) / h**2  # superdiagonal and subdiagonal
    Dxx = np.diag(d0, 0) + np.diag(dp, -1) + np.diag(dp, 1)

    # Fix first and last rows.
    Dxx[0, :4] = np.array([2, -5, 4, -1]) / h**2
    Dxx[-1, -4:] = np.array([-1, 4, -5, 2]) / h**2

    return x, Dx, Dxx


def diffcheb(n, xspan):
    """
    diffcheb(n, xspan)

    Compute Chebyshev differentiation matrices on n+1 points in the
    interval xspan. Return a vector of nodes, and the matrices for the first
    and second derivatives.
    """
    x = -np.cos(np.arange(n + 1) * np.pi / n)  # nodes in [-1,1]
    Dx = np.zeros([n + 1, n + 1])
    c = np.hstack([2.0, np.ones(n - 1), 2.0])  # endpoint factors

    # Off-diagonal entries
    Dx = np.zeros([n + 1, n + 1])
    for i in range(n + 1):
        for j in range(n + 1):
            if i != j:
                Dx[i, j] = (-1) ** (i + j) * c[i] / (c[j] * (x[i] - x[j]))

    # Diagonal entries by the "negative sum trick"
    for i in range(n + 1):
        Dx[i, i] = -np.sum([Dx[i, j] for j in range(n + 1) if j != i])

    # Transplant to [a,b]
    a, b = xspan
    x = a + (b - a) * (x + 1) / 2
    Dx = 2 * Dx / (b - a)

    # Second derivative
    Dxx = Dx @ Dx

    return x, Dx, Dxx


def bvplin(p, q, r, xspan, lval, rval, n):
    """
        bvplin(p, q, r, xspan, lval, rval, n)

    Use finite differences to solve a linear bopundary value problem. The ODE is
    u''+p(x)u'+q(x)u = r(x) on the interval xspan, with endpoint function
    values given as lval and rval. There will be n+1 equally spaced nodes,
    including the endpoints.

    Return vectors of the nodes and the solution values.
    """
    x, Dx, Dxx = diffmat2(n, xspan)

    P = np.diag(p(x))
    Q = np.diag(q(x))
    L = Dxx + P @ Dx + Q  # ODE expressed at the nodes

    # Replace first and last rows using boundary conditions.
    I = np.eye(n + 1)
    A = np.vstack([I[0], L[1:-1], I[-1]])
    b = np.hstack([lval, r(x[1:-1]), rval])

    # Solve the system.
    u = np.linalg.solve(A, b)

    return x, u

def bvp(phi, xspan, ga, gb, init):
    """
    bvp(phi, xspan, ga, gb, init)

    Use finite differences to solve a two-point boundary value problem. 
    The ODE is u'' = phi(x, u, u') for x in (a,b). The functions 
    ga(u(a), u'(a)) and gb(u(b), u'(b)) specify the boundary conditions. 
    The value init is an initial guess for [u(a), u'(a)].

    Return vectors for the nodes and the values of u.
    """
    n = len(init) - 1
    x, Dx, Dxx = diffmat2(n, xspan)
    h = x[1] - x[0]
    def residual(u):
        # Compute the difference between u'' and phi(x,u,u') at the
        # interior nodes and appends the error at the boundaries.
        du_dx = Dx @ u  # discrete u'
        d2u_dx2 = Dxx @ u  # discrete u''
        f = d2u_dx2 - phi(x, u, du_dx)

        # Replace first and last values by boundary conditions.
        f[0] = ga(u[0], du_dx[0]) / h
        f[n] = gb(u[n], du_dx[n]) / h
        return f

    u = levenberg(residual, init.copy())
    return x, u[-1]


def fem(c, s, f, a, b, n):
    """
    fem(c, s, f, a, b, n)

    Use a piecewise linear finite element method to solve a two-point boundary
    value problem. The ODE is (c(x)u')' + s(x)u = f(x) on the interval
    [a,b], and the boundary values are zero. The discretization uses n equal
    subintervals.

    Return vectors for the nodes and the values of u.
    """
    # Define the grid.
    h = (b - a) / n
    x = np.linspace(a, b, n + 1)

    # Templates for the subinterval matrix and vector contributions.
    Ke = np.array([[1, -1], [-1, 1]])
    Me = (1 / 6) * np.array([[2, 1], [1, 2]])
    fe = (1 / 2) * np.array([1, 1])

    # Evaluate coefficent functions and find average values.
    cval = c(x)
    cbar = (cval[:-1] + cval[1:]) / 2
    sval = s(x)
    sbar = (sval[:-1] + sval[1:]) / 2
    fval = f(x)
    fbar = (fval[:-1] + fval[1:]) / 2

    # Assemble global system, one interval at a time.
    K = np.zeros([n - 1, n - 1])
    M = np.zeros([n - 1, n - 1])
    f = np.zeros(n - 1)
    K[0, 0] = cbar[0] / h
    M[0, 0] = sbar[0] * h / 3
    f[0] = fbar[0] * h / 2
    K[-1, -1] = cbar[-1] / h
    M[-1, -1] = sbar[-1] * h / 3
    f[-1] = fbar[-1] * h / 2
    for k in range(1, n - 1):
        K[k - 1 : k + 1, k - 1 : k + 1] += (cbar[k] / h) * Ke
        M[k - 1 : k + 1, k - 1 : k + 1] += (sbar[k] * h) * Me
        f[k - 1 : k + 1] += (fbar[k] * h) * fe

    # Solve system for the interior values.
    u = np.linalg.solve(K + M, f)
    u = np.hstack([0, u, 0])  # put the boundary values into the result

    return x, u
