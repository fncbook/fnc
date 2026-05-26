---
numbering:
  enumerator: 13.4.%s
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---
```{code-cell}
:tags: [remove-cell]
from numpy import *
from scipy import linalg
from scipy.linalg import norm
from matplotlib.pyplot import *
from prettytable import PrettyTable
from timeit import default_timer as timer
import sys
sys.path.append('fncbook/')
import fncbook as FNC

# This (optional) block is for improving the display of plots.
# from IPython.display import set_matplotlib_formats
# set_matplotlib_formats("svg","pdf")
# %config InlineBackend.figure_format = 'svg'
rcParams["figure.figsize"] = [7, 4]
rcParams["lines.linewidth"] = 2
rcParams["lines.markersize"] = 4
rcParams['animation.html'] = "jshtml"  # or try "html5"
```

(section-twodim-nonlinear)=

# Nonlinear elliptic PDEs

Many nonlinear elliptic PDEs include references to the Laplacian operator.

::::{prf:example}
:label: example-mems2dmodel
Recall the micromechanical deflector modeled in a disk by {eq}`mems`. A fully two-dimensional equivalent is (see {cite}`peleskoEffectSmallaspectratio2006`)

```{math}
:label: mems2d
\Delta u - \frac{\lambda}{(u+1)^2} = 0.
```

This may be posed on any region, with $u=0$ specified everywhere on the boundary.
::::

More generally, we want to solve the nonlinear PDE

```{math}
:label: nonlinpdepde
\phi(x,y,u,u_x,u_y,u_{xx},u_{yy}) = 0
```

in the interior of a rectangle $R$, subject to the Dirichlet condition

```{math}
:label: nonlinpdebc
u(x,y) = g(x,y)
```

on the boundary of $R$.

## Implementation

```{index} quasi-Newton method
```

In order to solve for as few unknowns as possible, we use a Chebyshev discretization of the domain. The core idea is to formulate collocation equations at the grid points based on discrete approximations of {eq}`nonlinpdepde` and {eq}`nonlinpdebc`. If the PDE is nonlinear, then these equations are also nonlinear and are to be solved by a quasi-Newton iteration. @function-elliptic is our implementation.

``````{prf:algorithm} elliptic
:label: function-elliptic

```{literalinclude} chapter13.py
:filename: elliptic.py
:linenos: true
:language: python
:start-at: def elliptic
:end-at: return np.vectorize(evaluate)
```
``````

```{index} Julia; indexing arrays
```

```{note}
Take a moment to read through @function-elliptic slowly. Nearly every line is related to a topic in this book that took you some time to absorb. This is how most large problems in scientific computing are solved: by decomposing them until you reach steps that are recognizable as standard problems.
```

@function-elliptic first defines the discretization and then computes all the values of $g$ at the boundary nodes. It uses @function-levenberg as the nonlinear solver, and it translates back and forth between vector and grid shapes for the unknowns. After the discrete PDE is collocated at the grid points, the boundary terms are replaced by the boundary residual.

Lines 38–41, which produce the value returned by @function-elliptic, provide a function that evaluates the numerical solution anywhere in the domain, as is explained next.

## Off-grid evaluation

A Chebyshev grid is clustered close to the boundary of the domain. As a result, simple piecewise interpolation to evaluate off the grid, as is done by plotting routines, may be unacceptably nonsmooth and inaccurate. Instead, we should use the global polynomial interpolation that is foundational to the Chebyshev spectral method.

Let $\mathbf{U}$ be a matrix of solution values on the Chebyshev grid, defining a function $u(x,y)$, and let $(\xi,\eta)$ be a point where we want to evaluate $u(x,y)$. Column $\mathbf{u}_j$ of the grid matrix represents values spanning all the $x_i$ while $y$ is fixed at $y_j$. Therefore, we can define an interpolating polynomial $p_j(x)$ based on the values in $\mathbf{u}_j$.

Now let $v_j = p_j(\xi)$ for $j=1,\ldots,n$. The vector $\mathbf{v}$ is a discretization of $u(\xi,y)$ at the Chebyshev nodes in $y$. It defines an interpolating polynomial $q(y)$, and finally we have $u(\xi,\eta)=q(\eta)$. You can think of the total process as reducing one dimension at a time through the action of evaluating a polynomial interpolant at a point.

The function returned by @function-elliptic performs interpolation as described above, using a helper function `chebinterp` (not shown). The helper performs the evaluation of a polynomial interpolant in one variable using a modified implementation of {numref}`Function {number} <function-polyinterp>` that exploits the barycentric weights for Chebyshev nodes given in {eq}`weightcheb`.[^grideval]

[^grideval]: The interpolation algorithm in {numref}`Function {number} <function-elliptic>` is inefficient when $u$ is to be evaluated on a finer grid, as for plotting. A more careful version could re-use the same values $v_j = p_j(\xi)$ for multiple values of $\eta$.

::::{prf:example} MEMS model in 2D
:label: demo-nonlinear2d-mems

We solve the PDE

$$
\Delta u - \frac{\lambda}{(u+1)^2} = 0
$$

on the rectangle $[0,2.5] \times [0,1]$, with a zero Dirichlet condition on the boundary.

All we need to define are $\phi$ from {eq}`nonlinpdepde` for the PDE, and a trivial zero function for the boundary condition.

```{code-cell}
lamb = 1.5
phi = lambda x, y, u, ux, uxx, uy, uyy: uxx + uyy - lamb / (u + 1)**2
g = lambda x, y: 0
```

Here is the solution for $m=15$, $n=8$.

```{code-cell}
u = FNC.elliptic(phi, g, 15, [0, 2.5], 8, [0, 1])

print(f"solution at (2, 0.6) is {u(2, 0.6):.7f}")
```

```{code-cell}
:tags: hide-input
x = linspace(0, 2.5, 90)
y = linspace(0, 1, 60)
mtx, X, Y, _, _, _ = FNC.tensorgrid(x, y)
U = mtx(u)

pcolormesh(X.T, Y.T, U.T, cmap="viridis")
xlabel("$x$"),  ylabel("$y$"),  axis("equal")
colorbar()
title("Solution of the MEMS equation in 2d");
```

In the absence of an exact solution, how can we be confident that the solution is accurate? First, the Levenberg iteration converged without issuing a warning, so we should feel confident that the discrete equations were solved. We can check the boundary values easily. For example,

```{code-cell}
err = norm(u(x, 0) - g(x, 0), inf)
print(f"max error on bottom edge: {err:.2e}")
```

Assuming that we encoded the PDE correctly, the remaining source error is truncation from the discretization. We can estimate that by refining the grid a bit and seeing how much the numerical solution changes.

```{code-cell}
x_test = linspace(0, 2.5, 6)
y_test = linspace(0, 1, 6)
mtx_test, X_test, Y_test, _, _, _ = FNC.tensorgrid(x_test, y_test)

with printoptions(precision=7, suppress=True):
    print(mtx_test(u))
```

```{code-cell}
u = FNC.elliptic(phi, g, 25, [0, 2.5], 14, [0, 1])
with printoptions(precision=7, suppress=True):
    print(mtx_test(u))
```

The original solution seems to be accurate to about four digits.


::::

```{index} advection-diffusion equation
```

::::{prf:example} Steady advection-diffusion in 2D
:label: demo-nonlinear-advdiff

The steady-state limit of an advection-diffusion equation is

$$
1 - u_x - 2u_y + \epsilon \, \Delta u = 0.
$$

Here we solve it with a homogeneous Dirichlet condition on the square $[-1,1]^2$.


```{code-cell}
phi = lambda x, y, u, ux, uxx, uy, uyy: 1 - ux - 2*uy + 0.05 * (uxx + uyy)
g = lambda x, y: 0
u = FNC.elliptic(phi, g, 32, [-1, 1], 32, [-1, 1])
```

```{code-cell}
x = y = linspace(-1, 1, 70)
mtx, X, Y, _, _, _ = FNC.tensorgrid(x, y)
U = mtx(u)

pcolormesh(X.T, Y.T, U.T, cmap="viridis")
xlabel("$x$"),  ylabel("$y$"),  axis("equal")
colorbar()
title("Steady advection–diffusion");
```


::::

```{index} Allen–Cahn equation
```

::::{prf:example} Steady Allen–Cahn in 2D
:label: demo-nonlinear2d-allencahn

The stationary Allen–Cahn equation in two dimensions is

$$
u(1-u^2)+\epsilon \, \Delta u = 0.
$$


The following defines the PDE and a nontrivial Dirichlet boundary condition for the square $[0,1]^2$.

```{code-cell}
phi = lambda x, y, u, ux, uxx, uy, uyy: u * (1 - u**2) + 0.05 * (uxx + uyy)
g = lambda x, y: tanh(5 * (x + 2*y - 1))
```

We solve the PDE and then plot the result.

```{code-cell}
u = FNC.elliptic(phi, g, 36, [0, 1], 36, [0, 1])
```

```{code-cell}
:tags: hide-input
x = y = linspace(0, 1, 70)
mtx, X, Y, _, _, _ = FNC.tensorgrid(x, y)
U = mtx(u)
pcolormesh(X.T, Y.T, U.T, cmap="viridis")
xlabel("$x$"),  ylabel("$y$"),  axis("equal")
colorbar(),  title("Steady Allen–Cahn equation");
```

::::

## Exercises

``````{exercise}
:label: problem-nonlinear2d-steady
⌨ **(a)** Solve for the steady state of

$$
u_t = - u_y - x - 2 + \epsilon ( u_{xx} + u_{yy} )
$$

for $\epsilon=1$ in $[-1,1]\times[-1,1]$, subject to a homogeneous Dirichlet boundary condition. Use $m=n=30$ points and plot the solution.

**(b)** Repeat part (a) for $\epsilon=0.1$, which weakens the influence of diffusion relative to advection.
``````

```{index} soap film
```

``````{exercise}
:label: problem-nonlinear2d-soapfilm
⌨ A [soap film](wiki:Soap_film) stretched on a wire frame above the $(x,y)$ plane assumes a shape $u(x,y)$ of minimum area and is governed by

\begin{align*}
\operatorname{div} \, \left( \frac{\operatorname{grad} u}{\sqrt{1 + u_x^2 + u_y^2}} \right) &= 0 \text{ in region $R$},\\
u(x,y) &= g(x,y) \text{ on the boundary of $R$}.
\end{align*}

Solve the equation on $[-1,1]^2$ with boundary value $u(x,y)=\tanh(y-2x)$, and make a surface plot of the result. (Hints: Don't try to rewrite the PDE. Instead, modify @function-elliptic so that `ϕ` is called with arguments `(U,Dx,Dy)`, and compute the PDE in the form given. Also, since convergence is difficult in this problem, use the boundary data over the whole domain as the initial value for @function-levenberg.)
``````

``````{exercise}
:label: problem-nonlinear2d-mixedbc

⌨  Modify @function-elliptic to solve {eq}`nonlinpdepde` on $[a,b] \times [c,d]$ with the mixed boundary conditions

$$
u = 0, \text{ if } x=a \text{ or } y = d, \qquad  \frac{\partial u}{\partial n} = 0, \text{ if } x=b \text{ or } y = c,
$$

where $\frac{\partial}{\partial n}$ is the derivative in the direction of the outward normal. Either condition can be used at a corner point. (Hint: Define index vectors for each side of the boundary square.) Apply your solver to the PDE $\Delta u + \sin(3\pi x) = 0$ on $[0,1]^2$, and make a contour plot of the solution. Why do the level curves intersect two of the sides only at right angles?
``````
