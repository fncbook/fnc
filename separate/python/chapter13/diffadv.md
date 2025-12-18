---
numbering:
    enumerator: 13.2.%s
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

(section-twodim-diffadv)=

# Two-dimensional diffusion and advection

We next describe how to apply the method of lines to PDEs of the form

```{math}
:label: pde2d
  u_t = \phi(u,u_x,u_y,u_{xx},u_{xy},u_{yy}), \quad (x,y)\in [a,b]\times [c,d].
```

The PDE may be of either parabolic or hyperbolic type, with the primary difference being potential restrictions on the time step size. To keep descriptions and implementations relatively simple, we will consider only periodic conditions or Dirichlet boundary conditions.

As described in {numref}`section-twodim-tensorprod`, the rectangular domain is discretized by a grid $(x_i,y_j)$ for $i=0,\ldots,m$ and $j=0,\ldots,n$. The solution is semidiscretized as a matrix $\mathbf{U}(t)$ such that $U_{ij}$ is the approximate solution at $(x_i,y_j,t)$. Terms involving the spatial derivatives of $u$ are readily replaced by discrete counterparts, as shown in @tab-2d-derivatives.

```{table} Discrete replacements for 2D derivatives
:label: tab-2d-derivatives

| term | discrete replacement |
|:---|:---|
| $u$ | $\mathbf{U}$ |
| $u_x$ | $\mathbf{D}_x\mathbf{U}$ |
| $u_y$ | $\mathbf{U}\mathbf{D}_y^T$ |
| $u_{xx}$ | $\mathbf{D}_{xx}\mathbf{U}$ |
| $u_{yy}$ | $\mathbf{U}\mathbf{D}_{yy}^T$ |
```

## Matrix and vector shapes

```{index} ! vec,! Julia; vec operation, ! unvec, Julia; reshape
```

Our destination is an IVP that can be solved by a Runge–Kutta or multistep solver. These solvers are intended for vector problems, but our unknowns naturally have a matrix shape, which is the most convenient for the differentiation formulas {eq}`partfpartx` and {eq}`partfparty`. Fortunately, it's easy to translate back and forth between a matrix and an equivalent vector.

::::{prf:definition} vec and unvec operations
:label: definition-diffadv-vec
Let $\mathbf{A}$ be an $m\times n$ matrix. Define the **vec** function as stacking the columns of $\mathbf{A}$ into a vector, i.e.,

```{math}
:label: vecdef
\operatorname{vec}(\mathbf{A}) =
\begin{bmatrix}
A_{11} \\ \vdots \\ A_{m1}  \\ \vdots  \\ A_{1n} \\ \vdots \\ A_{m n}
\end{bmatrix}.
```

Let $\mathbf{z}$ be a vector of length $m n$. Define the **unvec** function as the inverse of vec:

```{math}
:label: unvecdef
\operatorname{unvec}(\mathbf{z}) = \begin{bmatrix}
  z_1 & z_{m+1} & \cdots & z_{m(n-1)+1} \\
  z_2 & z_{m+2} & \cdots & z_{m(n-1)+2} \\
  \vdots & \vdots & & \vdots \\
  z_m & z_{2m} & \cdots & z_{m n} \\
\end{bmatrix}.
```

::::

Suppose $\mathbf{U} = \operatorname{mtx}(u)$ is the matrix of unknowns. @tab-vec-unvec shows how to convert between $\mathbf{U}$ and $\mathbf{u} = \operatorname{vec}(\mathbf{U})$.

``````{table} vec and unvec operations
:label: tab-vec-unvec


|vec | unvec |
|:---|:---|
| `u = U.T.flatten()` | `U = np.reshape(u, (n+1, m+1)).T` |


```{dropdown} About the code
The incorporation of transposes is because NumPy uses row-major order, while MATLAB and Julia use column-major order. If performance were a concern, we would reverse our convention to avoid the transposes. We also use `flatten` rather than `ravel` to ensure making copies rather than views of the data and avoiding subtle bugs.
```


``````

::::{prf:example} Reshaping for grid functions
:label: demo-diffadv-vec


```{code-cell}
m, n = 4, 3
x = linspace(0, 2, m+1)
y = linspace(-3, 0, n+1)

f = lambda x, y: cos(0.75 * pi * x * y - 0.5 * pi * y)
mtx, X, Y, vec, unvec, _ = FNC.tensorgrid(x, y)
F = mtx(f)
print(f"function on a {m}x{n} grid:")
with printoptions(precision=4, suppress=True):
    print(F)

print("vec(F):")
with printoptions(precision=4, suppress=True):
    print(vec(F))
```

The `unvec` operation is the inverse of vec.

```{code-cell}
print("unvec(vec(F)):")
with printoptions(precision=4, suppress=True):
    print(unvec(vec(F)))
```

::::

In order to modularize our codes, we use @function-tensorgrid to define functions and values related to working with tensor-product grids. Its final output is to be discussed and used in @section-twodim-laplace.

``````{prf:algorithm} tensorgrid
:label: function-tensorgrid

```{literalinclude} chapter13.py
:filename: tensorgrid.py
:linenos: true
:language: python
:start-at: def tensorgrid
:end-at: return mtx
```
``````

## Periodic end conditions

If the boundary conditions are periodic, then the unknowns in the method of lines are the elements of the matrix $\mathbf{U}(t)$ representing grid values of the numerical solution. For the purposes of an IVP solution, this matrix is equivalent to the vector $\mathbf{u}(t)$ defined as $\mathbf{u}=\operatorname{vec}(\mathbf{U})$.

```{index} heat equation
```

::::{prf:example} Heat equation in 2D
:label: demo-diffadv-heat
We will solve a 2D heat equation, $u_t = 0.1(u_{xx} + u_{yy})$, on the square $[-1,1]\times[-1,1]$, with periodic behavior in both directions. The initial condition is

```{math}
:label: twodim-heatinit
u(x,y,0) = \sin(4\pi x) \exp[\cos(\pi y)], 
```

which is also periodic on the rectangle.


We start by defining the discretization of the rectangle.

```{code-cell}
m, n = 60, 40
x, Dx, Dxx = FNC.diffper(m, [-1, 1])
y, Dy, Dyy = FNC.diffper(n, [-1, 1])
mtx, X, Y, vec, unvec, _ = FNC.tensorgrid(x, y)
```

Here is the initial condition, evaluated on the grid.

```{code-cell}
u_init = lambda x, y: sin(4 * pi * x) * exp(cos(pi * y))
U0 = mtx(u_init)
mx = max(abs(U0))
pcolormesh(X, Y, U0, vmin=-mx, vmax=mx, cmap="RdBu", shading="gouraud")
axis("equal"),  colorbar()
xlabel("$x$"),  ylabel("$y$")
title("Initial condition");
```

The following function computes the time derivative for the unknowns, which have a vector shape. The actual calculations, however, take place using the matrix shape.

```{code-cell}
alpha = 0.1
def du_dt(t, u):
    U = unvec(u)
    Uyy = Dxx @ U
    Uxx = U @ Dyy.T
    dU_dt = alpha * (Uxx + Uyy)  # PDE
    return vec(dU_dt)
```

Since this problem is parabolic, a stiff time integrator is appropriate.

```{code-cell}
from scipy.integrate import solve_ivp
sol = solve_ivp(du_dt, (0, 0.2), vec(U0), method="BDF", dense_output=True)
U = lambda t: unvec(sol.sol(t))
```

We can use the function `U` defined above to get the solution at any time. Its output is a matrix of values on the grid. 

```{tip}
:class: dropdown
To plot the solution at any time, we use the same color scale as with the initial condition, so that the pictures are more easily compared.
```

```{code-cell}
pcolormesh(X.T, Y.T, U(0.02).T, 
    vmin=-mx, vmax=mx, cmap="RdBu", shading="gouraud")
axis("equal"),  colorbar()
xlabel("$x$"),  ylabel("$y$")
title("Heat equation, t=0.02");
```

An animation shows convergence toward a uniform value.

```{code-cell}
from matplotlib import animation
fig, ax = subplots()
obj = ax.pcolormesh(X, Y, U(0), vmin=-mx, vmax=mx, cmap="RdBu", shading="gouraud")
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
ax.set_xlabel("$x$"),  ax.set_ylabel("$y$")
ax.set_aspect("equal")
ax.set_title("Heat equation on a periodic domain")
def snapshot(t):
    global obj
    obj.remove()
    obj = ax.pcolormesh(X, Y, U(t), vmin=-mx, vmax=mx, cmap="RdBu", shading="gouraud")
    time_text.set_text(f"t = {t:.2f}")

anim = animation.FuncAnimation(fig, snapshot, frames=linspace(0, 0.2, 41))
anim.save("heat-2d.mp4", fps=30)
close()
```

![Heat equation in 2d](heat-2d.mp4)


::::

## Dirichlet conditions

```{index} boundary conditions; numerical implementation of
```

In {numref}`section-diffusion-boundaries` we coped with boundary conditions by removing the boundary values from the vector of unknowns being solved in the semidiscretized ODE. Each evaluation of the time derivative required us to extend the values to include the boundaries before applying differentiation matrices in space, then remove them from the time derivative vector.

We proceed similarly here, defining `chop` and `extend` functions that convert between the full grid and the inner grid of interior values. Mathematically speaking, the chopping operation is defined by

```{math}
:label: tensorprod-chop
\operatorname{chop}(\mathbf{U}) = \mathbf{E}_x \mathbf{U} \mathbf{E}_y^T,
```

where

$$
\mathbf{E}_x = \begin{bmatrix}
  0 & 1 & 0 & \cdots & 0 & 0 \\
  0 & 0 & 1 & \cdots & 0 & 0 \\
   &  &  & \ddots &  & \\
  0 & 0 & 0 & \cdots & 1 & 0
\end{bmatrix}
$$

is $(m-1)\times (m+1)$, and $\mathbf{E}_y$ is analogous but of size $(n-1)\times (n+1)$. The left multiplication in @tensorprod-chop deletes the first and last row of $\mathbf{U}$, and the right multiplication deletes its first and last column.

The extension operator is a bit more awkward to write out. It stars with appending rows and columns of zeros around the border of a matrix $\mathbf{W}$ of interior values:

```{math}
:label: tensorprod-extend
\tilde{\mathbf{U}} = \mathbf{E}_x^T \mathbf{W} \mathbf{E}_y
```

We can then modify the new zero values to reflect the boundary conditions, via

```{math}
:label: tensorprod-extend2
\mathbf{U} = \operatorname{extend}(\mathbf{W}) = \tilde{\mathbf{U}} + \mathbf{G},
```

Finally, we combine extending and chopping with the need to reshape between vectors and matrices. The function

$$
\operatorname{unpack}(\mathbf{w}) = \operatorname{extend}(\operatorname{unvec}(\mathbf{w}))
$$

converts a vector of unknowns (i.e., interior values) into a matrix of grid values, including the boundary values. The function

$$
\operatorname{pack}(\mathbf{U}) = \operatorname{vec}(\operatorname{chop}(\mathbf{U}))
$$

reverses the transformation, which is needed after the time derivative is computed on the full grid.

```{warning}
The `vec` and `unvec` reshaping operations in this context take place on the _interior_ grid, not the full grid. 
```

```{index} ! pack, Julia; pack operation, ! unpack, Julia; unpack operation
```
<!-- [^jacobian]: You might wonder why we use linear algebra to define the extension and deletion of boundary values rather than directly accessing row and column indices in the grid function. The algebraic expressions make it easier to express the Jacobian of a nonlinear The linear algebra approach allows `DifferentialEquations` to compute the Jacobian matrix of the implicit IVP solver quickly using *automatic differentiation* tools, greatly speeding up the solution process. Since the matrices in our expressions are sparse, multiplications by them do not affect running time perceptibly. -->

```{index} advection-diffusion equation
```

::::{prf:example} Advection-diffusion equation in 2D
:label: demo-diffadv-advdiff

We solve an advection-diffusion problem, $u_t + u_x = 1 + \epsilon(u_{xx} + u_{yy})$ on the square $[-1,1]^2$, with $u=0$ on the boundary. The outline of our approach is based on {numref}`Function {number} <function-parabolic>` for parabolic PDEs in one space dimension.


The first step is to define a discretization of the domain.

```{code-cell}
m, n = 50, 36
x, Dx, Dxx = FNC.diffcheb(m, [-1, 1])
y, Dy, Dyy = FNC.diffcheb(n, [-1, 1])
mtx, X, Y, _, _, _ = FNC.tensorgrid(x, y)
u_init = lambda x, y: (1 + y) * (1 - x)**4 * (1 + x)**2 * (1 - y**4)
```

We define functions `extend` and `chop` to deal with the Dirichlet boundary conditions. 

```{code-cell}
def chop(U):
    return U[1:-1, 1:-1]

def extend(W):
    U = zeros((m+1, n+1))
    U[1:-1, 1:-1] = W
    return U
```

Now we define the `pack` and `unpack` functions, using another call to @function-tensorgrid to get reshaping functions for the interior points. 

```{code-cell}
_, _, _, vec, unvec, _ = FNC.tensorgrid(x[1:-1], y[1:-1])
pack = lambda U: vec(chop(U))          # restrict to interior, then vectorize
unpack = lambda w: extend(unvec(w))    # reshape, then extend to boundary
```

Now we can define and solve the IVP using a stiff solver.

```{code-cell}
ep = 0.05
def dw_dt(t, w):
    U = unpack(w)
    Uyy = Dxx @ U
    Uxx = U @ Dyy.T 
    dU_dt = 1 - Dx @ U + ep * (Uxx + Uyy)
    return pack(dU_dt)

U0 = mtx(u_init)
sol = solve_ivp(dw_dt, (0, 2), pack(U0), method="BDF", dense_output=True)
```

When we evaluate the solution at a particular value of $t$, we get a vector of the interior grid values. The `unpack` converts this to a complete matrix of grid values.

```{code-cell}
U = lambda t: unpack(sol.sol(t))    # matrix-valued function of time

pcolormesh(X.T, Y.T, U(0.5).T, cmap="Blues", shading="gouraud")
colorbar()
xlabel("$x$"),  ylabel("$y$")
axis("equal"),  title("Solution at t=0.5");
```

```{code-cell}
mx = max([max(U(t)) for t in linspace(0, 2, 21)])
fig, ax = subplots()
obj = ax.pcolormesh(X.T, Y.T, U(0).T, vmin=0, vmax=2, cmap="Blues", shading="gouraud")
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
ax.set_xlabel("$x$"),  ax.set_ylabel("$y$")
ax.set_aspect("equal")
ax.set_title("Advection-diffusion in 2d")
def snapshot(t):
    global obj
    obj.remove()
    obj = ax.pcolormesh(X.T, Y.T, U(t).T, vmin=0, vmax=mx, cmap="Blues", shading="gouraud")
    time_text.set_text(f"t = {t:.2f}")

anim = animation.FuncAnimation(fig, snapshot, frames=linspace(0, 2, 81))
anim.save("advdiff-2d.mp4", fps=30)
close()
```

![Advection-diffusion in 2d](advdiff-2d.mp4)

::::

## Wave equation

```{index} wave equation
```

The wave equation introduces a little additional complexity. First, we write the 2D wave equation $u_{tt}=c^2(u_{xx}+u_{yy})$ in first-order form as

```{math}
:label: wave2dfirst
\begin{split}
    u_t &= v, \\
    v_t &= c^2(u_{xx}+u_{yy}).
\end{split}
```

Typical boundary conditions are to prescribe $u$ on the boundary and let $v$ be unspecified.

Now the grid functions are a pair of matrices $\mathbf{U}(t)$ and $\mathbf{V}(t)$. We need to chop $\mathbf{U}$ to an interior $\mathbf{W}$ and extend back using boundary data. Note that the IVP unknowns $\mathbf{W}$ and $\mathbf{V}$ have different sizes, so there are two separate reshaping operations involved. All of these details are handled within the `pack` and `unpack` functions we create.

::::{prf:example} Wave equation in 2D
:label: demo-diffadv-wave

We solve the wave equation with $c=1$ on the square $[-2,2]\times[-2,2]$, with $u=0$ on the boundary.

We start with the discretization and initial condition.

```{code-cell}
m, n = 40, 42
x, Dx, Dxx = FNC.diffcheb(m, [-2, 2])
y, Dy, Dyy = FNC.diffcheb(n, [-2, 2])
mtx, X, Y, _, _, _ = FNC.tensorgrid(x, y)

U0 = mtx(lambda x, y: (x + 0.2) * exp(-12 * (x**2 + y**2)))
V0 = zeros(U0.shape)
```

We need to define chopping and extension for the $u$ component. This looks the same as in @demo-diffadv-advdiff.

```{code-cell}
def extend(U):
    UU = zeros((m+1, n+1))
    UU[1:-1, 1:-1] = U
    return UU

def chop(U):
    return U[1:-1, 1:-1]
```

```{code-cell}
_, _, _, vec_v, unvec_v, _ = FNC.tensorgrid(x, y)
_, _, _, vec_u, unvec_u, _ = FNC.tensorgrid(x[1:-1], y[1:-1])
N = (m-1) * (n-1)    # number of interior points

def pack(U, V): 
    return hstack([vec_u(chop(U)), vec_v(V)])

def unpack(w):
    U = extend(unvec_u(w[:N]))
    V = unvec_v(w[N:])
    return U, V
```

We can now define and solve the IVP. Since this problem is hyperbolic, a nonstiff integrator is faster than a stiff one.

```{code-cell}
def dw_dt(t, w):
    U, V = unpack(w)
    dU_dt = V
    dV_dt = Dxx @ U + U @ Dyy.T
    return pack(dU_dt, dV_dt)

from scipy.integrate import solve_ivp
sol = solve_ivp(dw_dt, (0, 4), pack(U0, V0), method="RK45", dense_output=True)
U = lambda t: unpack(sol.sol(t))[0]
```

```{code-cell}
fig, ax = subplots()
obj = ax.pcolormesh(X, Y, U(0), vmin=-0.1, vmax=0.1, cmap="RdBu", shading="gouraud")
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
ax.set_xlabel("$x$"),  ax.set_ylabel("$y$")
ax.set_aspect("equal")
ax.set_title("Wave equation in 2d")

def snapshot(t):
    global obj
    obj.remove()
    obj = ax.pcolormesh(X, Y, U(t), vmin=-0.1, vmax=0.1, cmap="RdBu", shading="gouraud")
    time_text.set_text(f"t = {t:.2f}")

anim = animation.FuncAnimation(fig, snapshot, frames=linspace(0, 4, 91))
anim.save("wave-2d.mp4", fps=30);
close()
```

![Wave equation in 2d](wave-2d.mp4)

::::

## Exercises

``````{exercise}
:label: problem-diffadv-derivatives
⌨  For the given $u(x,y)$, make a plot of the given quantity on the square $[-2,2]^2$ using appropriate differentiation matrices.

**(a)** $u(x,y) = \exp(x-y^2)$; plot $u_{xx}+u_{yy}$

**(b)** $u(x,y) =\cos (\pi x)+\sin (\pi y)$; plot $u_x+u_y$

**(c)** $u(x,y) =\exp(-x^2-4y^2)$; plot $x u_y$
``````

``````{exercise}
:label: problem-diffadv-allencahn
⌨ Following @demo-diffadv-heat as a model, solve the [Allen–Cahn equation](wiki:Allen–Cahn_equation) $u_t=u(1-u^2)+0.001(u_{xx}+u_{yy})$ on the square $[-1,1]^2$ with periodic conditions, taking $u(x,y,0)=\sin(\pi x)\cos(2\pi y)$. Use $m=n=60$ to solve up to $t=4$, and make an animation of the result.
``````

``````{exercise}
:label: problem-diffadv-advdiff2d

⌨ Following @demo-diffadv-advdiff as a model, solve $u_t=y u_x-u_y+0.03(u_{xx}+u_{yy})$ on the square $[-1,1]^2$, with $u(x,y,0)=(1-x^2)(1-y^2)$ and homogeneous Dirichlet boundary conditions. Use $m=n=40$ to solve up to $t=2$, and make an animation of the result.
``````

``````{exercise}
:label: problem-diffadv-wave2d
⌨ Following @demo-diffadv-wave as a model, solve $u_{tt}=u_{xx}+u_{yy}+\cos(7t)$ on the square $[-1,1]^2$, with $u(x,y,0)=x(1-x^6)(1-y^2)$, $u_t(x,y,0)=0$, subject to homogeneous Dirichlet boundary conditions. Take $m=n=60$ to solve between $t=0$ and $t=12$, and make an animation of the result.
``````

``````{exercise}
:label: problem-diffadv-maxwell
From Maxwell's equations we can find a way to convert the wave equation to a first-order form that, unlike {eq}`wave2dfirst`, uses only first-order derivatives in space:

```{math}
:label: wave2dTM
\begin{split}
u_t &= c^2(v_y - w_x),\\
v_t &= u_y, \\
w_t &= -u_x,
\end{split}
```

subject to $u=0$ on the boundary.

**(a)** ✍ Show that a solution of {eq}`wave2dTM` satisfies $u_t=c^2(u_{xx}+u_{yy})$.

**(b)** ⌨ Solve {eq}`wave2dTM` with $c=2$ in the rectangle $x\in[-3,3]$, $y\in[-1,1]$, $u(x,y,0) = \exp(x-x^2)(9-x^2)(1-y^2)$, and $v=w=0$ at $t=0$. Use $m=50$ for $x$ and $n=25$ for $y$, solve for $0\le t \le 6$, and make an animation of the solution.
``````
