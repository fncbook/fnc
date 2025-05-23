---
kernelspec:
  display_name: Python 3
  language: python
  name: python3
numbering:
  headings: false
---
# Chapter 13

## Functions

(function-tensorgrid-python)=
``````{dropdown} Create a tensor-product grid
:open:
```{literalinclude} fncbook/fncbook/chapter13.py
:filename: tensorgrid.py
:linenos: true
:language: python
:start-at: def tensorgrid
:end-at: return mtx
```
``````

(function-poissonfd-python)=
``````{dropdown} Solution of Poisson's equation by finite differences
:open:
```{literalinclude} fncbook/fncbook/chapter13.py
:filename: poissonfd.py
:linenos: true
:language: python
:start-at: def poissonfd
:end-at: return U
```
``````

(function-elliptic-python)=
``````{dropdown} Solution of elliptic PDE by Chebyshev collocation
:open:
```{literalinclude} fncbook/fncbook/chapter13.py
:filename: elliptic.py
:linenos: true
:language: python
:start-at: def elliptic
:end-at: return np.vectorize(evaluate)
```
``````

## Examples

```{code-cell}
:tags: remove-cell
exec(open("FNC_init.py").read())
```

### 13.1 @section-twodim-tensorprod

(demo-tensorprod-gridfun-python)=
``````{dropdown} @demo-tensorprod-gridfun
:open:
Here is the grid from {numref}`Example {number} <example-tensorprod-smallgrid>`.

```{code-cell}
m = 4
x = linspace(0, 2, m+1)
n = 2
y = linspace(1, 3, n+1)
```

For a given $f(x,y)$ we can find $\operatorname{mtx}(f)$ by using a comprehension syntax.

```{code-cell}
f = lambda x, y: cos(pi * x * y - y)
F = array( [ [f(xi, yj) for yj in y] for xi in x ] )
print(F)
```

We can make a nice plot of the function by first choosing a much finer grid. However, the contour and surface plotting functions expect the *transpose* of mtx($f$).
```{tip}
:class: dropdown
To emphasize departures from a zero level, use a colormap such as `RdBu` and set the color limits to be symmetric around zero.
```

::::{warning}
The contour and surface plotting functions in `matplotlib` expect the *transpose* of the outputs of `mtx`. If you forget to do that, the $x$ and $y$ axes will be swapped.
::::

```{code-cell}
m, n = 80, 70
x = linspace(0, 2, m+1)
y = linspace(1, 3, n+1)
mtx, X, Y, _, _, _ = FNC.tensorgrid(x, y)
F = mtx(f)

pcolormesh(X.T, Y.T, F.T, cmap="RdBu", vmin=-1, vmax=1, shading="gouraud")
axis("equal"),  colorbar()
xlabel("$x$"),  ylabel("$y$");
```
``````

(demo-tensorprod-disksphere-python)=
``````{dropdown} @demo-tensorprod-disksphere
:open:
For a function given in polar form, such as $f(r,\theta)=1-r^4$, construction of a function over the unit disk is straightforward using a grid in $(r,\theta)$ space.

```{code-cell}
r = linspace(0, 1, 41)
theta = linspace(0, 2*pi, 121)
mtx, R, Theta, _, _, _ = FNC.tensorgrid(r, theta)

F = mtx(lambda r, theta: 1 - r**4)    

contourf(R.T, Theta.T, F.T, levels=20)
colorbar()
xlabel("$r$"),  ylabel("$\\theta$");
```

Of course, we are used to seeing such plots over the $(x,y)$ plane, not the $(r,\theta)$ plane. For this we create matrices for the coordinate functions $x$ and $y$.

```{code-cell}
X, Y = R * cos(Theta), R * sin(Theta)
contourf(X.T, Y.T, F.T, levels=20)
colorbar(),  axis("equal")
xlabel("$x$"),  ylabel("$y$");
```

In such functions the values along the line $r=0$ must be identical, and the values on the line $\theta=0$ should be identical to those on $\theta=2\pi$. Otherwise the interpretation of the domain as the unit disk is nonsensical. If the function is defined in terms of $x$ and $y$, then those can be defined in terms of $r$ and $\theta$ using {eq}`unitdiskparam`.

``````

(demo-tensorprod-diff-python)=
``````{dropdown} @demo-tensorprod-diff
:open:
We define a function and, for reference, its two exact partial derivatives.

```{code-cell}
u = lambda x, y: sin(pi * x * y - y)
du_dx = lambda x, y: pi * y * cos(pi * x * y - y)
du_dy = lambda x, y: (pi * x - 1) * cos(pi * x * y - y)
```

We will use an equispaced grid and second-order finite differences as implemented by `diffmat2`. First, we have a look at a plots of the exact partial derivatives.

```{code-cell}
m, n = 80, 70
x, Dx, Dxx = FNC.diffmat2(m, [0, 2])
y, Dy, Dyy = FNC.diffmat2(n, [1, 3])
mtx, X, Y, _, _, _ = FNC.tensorgrid(x, y)

U = mtx(u)
dU_dX = mtx(du_dx)
dU_dY = mtx(du_dy)

subplot(1, 2, 1)
contourf(X, Y, dU_dX)
title("$u$"),  axis("equal")
subplot(1, 2, 2)
contourf(X, Y, dU_dY)
title("$\\partial u/\\partial y$"),  axis("equal");
```

Now we compare the exact partial derivatives with their finite-difference approximations. Since these are signed errors, we use a colormap that is symmetric around zero.

```{code-cell}
subplot(1, 2, 1)
err_x = Dx @ U - dU_dX
M = max(abs(err_x))
pcolormesh(X, Y, err_x, shading="gouraud", cmap="RdBu", vmin=-M, vmax=M)
colorbar()
title("error in $\\partial u/\\partial x$"),  axis("equal")

subplot(1, 2, 2)
err_y = U @ Dy.T - dU_dY
M = max(abs(err_y))
pcolormesh(X, Y, err_y, shading="gouraud", cmap="RdBu", vmin=-M, vmax=M)
colorbar()
title("error in $\\partial u/\\partial y$"),  axis("equal");
```

Not surprisingly, the errors are largest where the derivatives themselves are largest.
``````

### 13.2 @section-twodim-diffadv

(demo-diffadv-vec-python)=
``````{dropdown} @demo-diffadv-vec
:open:

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
``````

(demo-diffadv-heat-python)=
``````{dropdown} @demo-diffadv-heat
:open:

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
anim.save("figures/heat-2d.mp4", fps=30)
close()
```

![Heat equation in 2d](figures/heat-2d.mp4)

``````

(demo-diffadv-advdiff-python)=
``````{dropdown} @demo-diffadv-advdiff
:open:

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
anim.save("figures/advdiff-2d.mp4", fps=30)
close()
```

![Advection-diffusion in 2d](figures/advdiff-2d.mp4)
``````

(demo-diffadv-wave-python)=
``````{dropdown} @demo-diffadv-wave
:open:
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
anim.save("figures/wave-2d.mp4", fps=30);
close()
```

![Wave equation in 2d](figures/wave-2d.mp4)
``````

### 13.3 @section-twodim-laplace

(demo-laplace-kron-python)=
``````{dropdown} @demo-laplace-kron
:open:

```{code-cell}
A = array([[1, 2], [-2, 0]])
B = array([[1, 10, 100], [-5, 5, 3]])
print("A:")
print(A)
print("B:")
print(B)
```

Applying the definition manually, we get

```{code-cell}
A_kron_B = vstack([ hstack([A[0, 0] * B, A[0, 1] * B]), hstack([A[1, 0] * B, A[1, 1] * B]) ])
print(A_kron_B)
```

```{index} ! Python; kron
```

But it makes more sense to use `kron` from NumPy, or the `scipy.sparse` version when sparsity is to be preserved.

```{code-cell}
kron(A, B)
```
``````

(demo-laplace-fd-python)=
``````{dropdown} @demo-laplace-fd
:open:
We make a crude discretization for illustrative purposes.

```{code-cell}
m, n = 5, 6
x, Dx, Dxx = FNC.diffmat2(m, [0, 3])
y, Dy, Dyy = FNC.diffmat2(n, [-1, 1])
mtx, X, Y, vec, unvec, is_boundary = FNC.tensorgrid(x, y)
```

Here is a look at the matrix we called $\mathbf{L}$ (the discrete Laplacian), before any modifications are made for the boundary conditions. The combination of Kronecker products and finite differences produces a characteristic sparsity pattern.

```{code-cell}
import scipy.sparse as sp
Dxx = sp.lil_matrix(Dxx)
Dyy = sp.lil_matrix(Dyy)
Ix = sp.eye(m+1)
Iy = sp.eye(n+1)
A = sp.kron(Iy, Dxx) + sp.kron(Dyy, Ix)

spy(A)
title("Matrix before imposing BC");
```

The number of equations is equal to $(m+1)(n+1)$, which is the total number of points on the grid.

```{code-cell}
N = (m+1) * (n+1)
```

We now use the final output from @function-tensorgrid, which is a Boolean array indicating where the boundary points lie in the grid.

```{code-cell}
spy(is_boundary)
title("Boundary points");
```

In order to impose Dirichlet boundary conditions, we use the boundary indicator to index into the rows of the system.

```{tip}
:class: dropdown
Changing rows of a sparse array requires that the operands be in a particular sparse representation called `lil`. The conversion isn't done automatically because it can be slow and you are encouraged to avoid it when possible. We're just trying to keep things conceptually simple here.
```

```{code-cell}
I = sp.eye(N, format="lil")
idx = vec(is_boundary)
A = A.tolil()
A[idx, :] = I[idx, :];    # Dirichlet conditions

spy(A)
title("Matrix with Dirichlet BC imposed");
```

Next, we evaluate $\phi$ on the grid to get the forcing vector of the linear system, and then modify the boundary rows to hold the boundary values—in this case, zero.

```{code-cell}
f = lambda x, y: x**2 - y + 2
b = vec(mtx(f))
b[idx] = 0
```

Now we can solve for $\mathbf{u}$ and reinterpret it as the matrix-shaped $\mathbf{U}$, the solution on our grid.

```{code-cell}
from scipy.sparse.linalg import spsolve
u = spsolve(A.tocsr(), b)
U = unvec(u)
with printoptions(precision=4, suppress=True):
    print(U)
```
``````

(demo-laplace-poisson-python)=
``````{dropdown} @demo-laplace-poisson
:open:

First we define the problem on $[0,1]\times[0,2]$.

```{code-cell}
f = lambda x, y: -sin(3 * x * y - 4 * y) * (9 * y**2 + (3 * x - 4) ** 2)
g = lambda x, y: sin(3 * x * y - 4 * y)
xspan = [0, 1]
yspan = [0, 2]
```

Here is the finite-difference solution.

```{code-cell}
U, X, Y = FNC.poissonfd(f, g, 50, xspan, 80, yspan)
```

```{code-cell}
:tags: hide-input
pcolormesh(X.T, Y.T, U.T, cmap="viridis")
xlabel("$x$"),  ylabel("$y$"),  axis("equal")
colorbar(), title("Solution of Poisson's equation");
```

Since this is an artificial problem with a known solution, we can plot the error, which is a smooth function of $x$ and $y$. It must be zero on the boundary; otherwise, we have implemented boundary conditions incorrectly.

```{code-cell}
error = g(X, Y) - U    # because we set up g as the exact solution
M = max(abs(error))

pcolormesh(X.T, Y.T, error.T, vmin=-M, vmax=M, cmap="RdBu")
xlabel("$x$"),  ylabel("$y$"),  axis("equal")
colorbar(),  title("Error");
```
``````

### 13.4 @section-twodim-nonlinear

(demo-nonlinear2d-mems-python)=
``````{dropdown} @demo-nonlinear2d-mems
:open:
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

``````

(demo-nonlinear-advdiff-python)=
``````{dropdown} @demo-nonlinear-advdiff
:open:

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

``````

(demo-nonlinear2d-allencahn-python)=
``````{dropdown} @demo-nonlinear2d-allencahn
:open:

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
``````