---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
# Chapter 13

## Functions

(function-tensorgrid-julia)=
``````{dropdown} Create a tensor-product grid
:open:
```{literalinclude} FNCFunctions/src/chapter13.jl
:filename: tensorgrid.jl
:linenos: true
:language: julia
:start-after: # begin tensorgrid
:end-before: # end tensorgrid
```
``````

(function-poissonfd-julia)=
``````{dropdown} Solution of Poisson's equation by finite differences
:open:
```{literalinclude} FNCFunctions/src/chapter13.jl
:filename: poissonfd.jl
:linenos: true
:language: julia
:start-after: # begin poissonfd
:end-before: # end poissonfd
```
``````

(function-elliptic-julia)=
``````{dropdown} Solution of elliptic PDE by Chebyshev collocation
:open:
```{literalinclude} FNCFunctions/src/chapter13.jl
:filename: elliptic.jl
:linenos: true
:language: julia
:start-after: # begin elliptic
:end-before: # end elliptic
```
::::{admonition} About the code
:class: dropdown
The boundary values are accessed using Boolean indexing. One advantage of this style, though it is not exploited here, is that the complementary points can also be accessed via the Boolean NOT operator `!`. Note that any indexing array either has to be the same size as the object of the indexing, or a vector with the same number of elements. In this function, for example, `X[idx]`, `X[isboundary]`, and `u[idx]` would all be valid, but `u[isboundary]` would not be.
::::
``````

## Examples

```{code-cell}
:tags: remove-output
include("FNC_init.jl")
```

### 13.1 @section-twodim-tensorprod

(demo-tensorprod-gridfun-julia)=
``````{dropdown} @demo-tensorprod-gridfun
:open:
Here is the grid from {numref}`Example {number} <example-tensorprod-smallgrid>`.

```{code-cell}
m = 4
x = range(0, 2, m+1)
n = 2
y = range(1, 3, n+1);
```

For a given $f(x,y)$ we can find $\operatorname{mtx}(f)$ by using a comprehension syntax.

```{code-cell}
f = (x, y) -> cos(π * x * y - y)
F = [f(x, y) for x in x, y in y]
```

We can make a nice plot of the function by first choosing a much finer grid. However, the contour and surface plotting functions expect the *transpose* of mtx($f$).
```{tip}
:class: dropdown
To emphasize departures from a zero level, use a colormap such as `redsblues`, and use `clims` to set balanced color differences.
```

```{code-cell}
using Plots
m, n = 80, 60
x = range(0, 2, m+1);
y = range(1, 3, n+1);
F = [f(x, y) for x in x, y in y]
contour(x, y, F';
    levels=21,  aspect_ratio=1,
    color=:redsblues,  clims=(-1, 1),
    xlabel="x",  ylabel="y" )
``` 

``````

(demo-tensorprod-disksphere-julia)=
``````{dropdown} @demo-tensorprod-disksphere
:open:
For a function given in polar form, such as $f(r,\theta)=1-r^4$, construction of a function over the unit disk is straightforward using a grid in $(r,\theta)$ space.

```{code-cell}
r = range(0, 1, 41)
θ = range(0, 2π, 81)
F = [1 - r^4 for r in r, θ in θ]
plot(r, θ, F';
    legend=:none, 
    color=:viridis,  fill=true,
    xlabel="r",  ylabel="θ", 
    title="A polar function")
```

Of course, we are used to seeing such plots over the $(x,y)$ plane, not the $(r,\theta)$ plane. 

In such functions the values along the line $r=0$ must be identical, and the values on the line $\theta=0$ should be identical to those on $\theta=2\pi$. Otherwise the interpretation of the domain as the unit disk is nonsensical. If the function is defined in terms of $x$ and $y$, then those can be defined in terms of $r$ and $\theta$ using {eq}`unitdiskparam`.

``````

(demo-tensorprod-diff-julia)=
``````{dropdown} @demo-tensorprod-diff
:open:
We define a function and, for reference, its two exact partial derivatives.

```{code-cell}
u = (x, y) -> sin(π * x * y - y);
∂u_∂x = (x, y) -> π * y * cos(πx * y - y);
∂u_∂y = (x, y) -> (π * x - 1) * cos(π * x * y - y);
```

We use an equispaced grid and second-order finite differences as implemented by `diffmat2`.

```{code-cell}
m = 80;
x, Dx, _ = FNC.diffmat2(m, [0, 2]);
n = 60;
y, Dy, _ = FNC.diffmat2(n, [1, 3]);
mtx = (f, x, y) -> [f(x, y) for x in x, y in y]
U = mtx(u, x, y)
∂xU = Dx * U
∂yU = U * Dy';
```

Now we compare the exact $\frac{\partial u}{\partial y}$ with its finite-difference approximation.

```{code-cell}
M = maximum(abs, ∂yU)    # find the range of the result
plot(layout=(1, 2), 
    aspect_ratio=1,   clims=(-M, M), 
    xlabel="x", ylabel="y")
contour!(x, y, mtx(∂u_∂y, x, y)';
    levels=15,  subplot=1,
    color=:redsblues,
    title="∂u/∂y")
contour!(x, y, ∂yU';
    levels=15,  subplot=2,
    color=:redsblues, 
    title="approximation")
```

To the eye there is little difference to be seen, though the results have no more than a few correct digits at these discretization sizes:

```{code-cell}
exact = mtx(∂u_∂y, x, y)
# Relative difference in Frobenius norm:
norm(exact - ∂yU) / norm(exact)
```
``````

### 13.2 @section-twodim-diffadv

(demo-diffadv-vec-julia)=
``````{dropdown} @demo-diffadv-vec
:open:

```{code-cell}
m = 2;
n = 3;
V = rand(1:9, m, n);
v = vec(V)
```

The `unvec` operation is the inverse of vec.

```{code-cell}
unvec = z -> reshape(z, m, n)
unvec(v)
```
``````

(demo-diffadv-heat-julia)=
``````{dropdown} @demo-diffadv-heat
:open:

```{code-cell}
m, n = (60, 25)
x, Dx, Dxx = FNC.diffper(m, [-1, 1])
y, Dy, Dyy = FNC.diffper(n, [-1, 1])
mtx = f -> [f(x, y) for x in x, y in y]
unvec = z -> reshape(z, m, n);
```

Note that the initial condition should also be periodic on the domain.

```{code-cell}
using Plots
u_init = (x, y) -> sin(4 * π * x) * exp(cos(π * y))
U₀ = mtx(u_init)
M = maximum(abs, U₀)
contour(x, y, U₀';
    color=:redsblues,  clims=(-M, M), 
    aspect_ratio=1,
    xaxis=("x", (-1, 1)),  yaxis=("y", (-1, 1)), 
    title="Initial condition" )
```

This function computes the time derivative for the unknowns. The actual calculations take place using the matrix shape.

```{code-cell}
function du_dt(u, α, t)
    U = unvec(u)
    Uxx = Dxx * U
    Uyy = U * Dyy'            # 2nd partials
    du_dt = α * (Uxx + Uyy)    # PDE
    return vec(du_dt)
end;
```

Since this problem is parabolic, a stiff integrator is appropriate.

```{code-cell}
using OrdinaryDiffEq
IVP = ODEProblem(du_dt, vec(U₀), (0, 0.2), 0.1)
sol = solve(IVP, Rodas4P());
```

Here is an animation of the solution.
```{tip}
:class: dropdown
Here `clims` are set so that colors remain at fixed values throughout the animation.
```

```{code-cell}
:tags: remove-output, hide-input
anim = @animate for t in range(0, 0.2, 81)
    surface(x, y, unvec(sol(t))';
        color=:redsblues,  clims=(-M, M),
        xaxis=(L"x", (-1, 1)), 
        yaxis=(L"y", (-1, 1)), 
        zlims=(-M, M),
        title=@sprintf("Heat equation, t=%.3f", t),
        dpi=150, colorbar=:none)
end
mp4(anim, "figures/2d-heat.mp4");
```

![Heat equation in 2d](figures/2d-heat.mp4)

``````

(demo-diffadv-advdiff-julia)=
``````{dropdown} @demo-diffadv-advdiff
:open:

The first step is to define a discretization of the domain. 

```{code-cell}
m, n = 50, 36
x, Dx, Dxx = FNC.diffcheb(m, [-1, 1])
y, Dy, Dyy = FNC.diffcheb(n, [-1, 1])
mtx, X, Y, _ = FNC.tensorgrid(x, y)
U₀ = mtx( (x, y) -> (1 + y) * (1 - x)^4 * (1 + x)^2 * (1 - y^4) );
```

There are really two grids now: the full grid and the subset grid of interior points. Since the IVP unknowns are on the interior grid, that is the one we need to change shapes on. We also need the functions `extend` and `chop` to add and remove boundary values.

```{code-cell}
chop = U -> U[2:m, 2:n]
extend = U -> [zeros(m+1) [zeros(1, n-1); U; zeros(1, n-1)] zeros(m+1)]
unvec = u -> reshape(u, m-1, n-1)
pack = U -> vec(chop(U))
unpack = w -> extend(unvec(w))
```

Now we can define and solve the IVP using a stiff solver.

```{code-cell}
function dw_dt(w, ϵ, t)
    U = unpack(w)
    Ux, Uxx = Dx * U, Dxx * U
    Uyy = U * Dyy'
    du_dt = @. 1 - Ux + ϵ * (Uxx + Uyy)
    return pack(du_dt)
end

IVP = ODEProblem(dw_dt, pack(U₀), (0.0, 2), 0.05)
w = solve(IVP, Rodas4P());
```

When we evaluate the solution at a particular value of $t$, we get a vector of the interior grid values. The same `unpack` function above converts this to a complete matrix of grid values.

```{code-cell}
U = t -> unpack(w(t))
contour(x, y, U(0.5)';
    fill=true,  color=:blues,  levels=20, l=0,
    aspect_ratio=1,  xlabel=L"x",  ylabel=L"y",
    title="Solution at t = 0.5")
```

```{code-cell}
:tags: remove-output, hide-input
anim = @animate for t in 0:0.02:2
    U = unpack(w(t))
    surface(x, y, U';
        layout=(1, 2),  size=(640, 320),
        xlabel=L"x",  ylabel=L"y",  zaxis=((0, 2), L"u(x,y)"),
        color=:blues,  alpha=0.66,  clims=(0, 2), colorbar=:none,
        title="Advection-diffusion",  dpi=150)
    contour!(x, y, U'; 
        levels=24, 
        aspect_ratio=1,  subplot=2, 
        xlabel=L"x",  ylabel=L"y",
        color=:blues,  clims=(0, 2),  colorbar=:none,
        title=@sprintf("t = %.2f", t))
end
mp4(anim, "figures/2d-advdiff.mp4");
```

![Advection-diffusion in 2d](figures/2d-advdiff.mp4)
``````

(demo-diffadv-wave-julia)=
``````{dropdown} @demo-diffadv-wave
:open:
We start with the discretization and initial condition.

```{code-cell}
m, n = 40, 40
x, Dx, Dxx = FNC.diffcheb(m, [-2, 2])
y, Dy, Dyy = FNC.diffcheb(n, [-2, 2])
mtx, X, Y, _ = FNC.tensorgrid(x, y)
U₀ = mtx( (x, y) -> (x + 0.2) * exp(-12 * (x^2 + y^2)) )
V₀ = zeros(size(U₀));
```

Note that because $u$ is known on the boundary, while $v$ is unknown over the full grid, there are two different sizes of vec/unvec operations. We also need to define functions to pack grid unknowns into a vector and to unpack them. When the unknowns for $u$ are packed, the boundary values are chopped off, and these are restored when unpacking.

```{code-cell}
_, _, _, unvec_v, _ = FNC.tensorgrid(x, y)
_, _, _, unvec_u, _ = FNC.tensorgrid(x[2:m], y[2:n])
chop = U -> U[2:m, 2:n]
extend = U -> [zeros(m+1) [zeros(1, n-1); U; zeros(1, n-1)] zeros(m+1)]
pack = (U, V) -> [vec(chop(U)); vec(V)]
N = (m-1) * (n-1)    # number of interior unknowns
unpack = w -> ( extend(unvec_u(w[1:N])), unvec_v(w[N+1:end]) )
```

We can now define and solve the IVP. Since this problem is hyperbolic, not parabolic, a nonstiff integrator is faster than a stiff one.

```{code-cell}
function dw_dt(w, c, t)
    U, V = unpack(w)
    du_dt = V
    dv_dt = c^2 * (Dxx * U + U * Dyy')
    return pack(du_dt, dv_dt)
end

IVP = ODEProblem(dw_dt, pack(U₀, V₀), (0, 4.0), 1)
sol = solve(IVP, Tsit5())
U = t -> unpack(sol(t))[1]
```

```{code-cell}
:tags: remove-output, hide-input
anim = @animate for t in 0:4/100:4
    Ut = U(t)
    surface(x, y, Ut';
        layout=(1, 2), size=(640, 320),
        xlabel=L"x",  ylabel=L"y",  zaxis=((-0.1, 0.1), L"u(x,y)"),
        color=:redsblues,  alpha=0.66,  clims=(-0.1, 0.1), colorbar=:none,
        title="Wave equation",  dpi=150)
    contour!(x, y, Ut'; 
        levels=24,  subplot=2, 
        aspect_ratio=1,
        xlabel=L"x",  ylabel=L"y",
        color=:redsblues,  clims=(-0.1, 0.1), 
        colorbar=:none,  title=@sprintf("t = %.2f", t))
end
mp4(anim, "figures/2d-wave.mp4");
```

![Wave equation in 2d](figures/2d-wave.mp4)
``````

### 13.3 @section-twodim-laplace

(demo-laplace-kron-julia)=
``````{dropdown} @demo-laplace-kron
:open:

```{code-cell}
A = [1 2; -2 0]
```

```{code-cell}
B = [1 10 100; -5 5 3]
```

Applying the definition manually, we get

```{code-cell}
A_kron_B = [
    A[1, 1]*B A[1, 2]*B;
    A[2, 1]*B A[2, 2]*B
    ]
```

```{index} ! Julia; kron
```

That result should be the same as the following.

```{code-cell}
kron(A, B)
```
``````

(demo-laplace-fd-julia)=
``````{dropdown} @demo-laplace-fd
:open:
We make a crude discretization for illustrative purposes.

```{code-cell}
m, n = 6, 5
x, Dx, Dxx = FNC.diffmat2(m, [0, 3])
y, Dy, Dyy = FNC.diffmat2(n, [-1, 1])
mtx, X, Y, unvec, is_boundary = FNC.tensorgrid(x, y)
```

Next, we evaluate $\phi$ on the grid to get the forcing vector of the linear system.

```{code-cell}
ϕ = (x, y) -> x^2 - y + 2
b = vec(mtx(ϕ));
```

Here are the coefficients for the PDE collocation, before any modifications are made for the boundary conditions. The combination of Kronecker products and finite differences produces a characteristic sparsity pattern.

```{code-cell}
using SparseArrays, Plots
A = kron(I(n+1), sparse(Dxx)) + kron(sparse(Dyy), I(m+1))
spy(A;
    color=:blues,  m=3,
    title="System matrix before boundary conditions")
```

The number of equations is equal to $(m+1)(n+1)$, which is the total number of points on the grid.

```{code-cell}
N = length(b)
```

The combination of Kronecker products and finite differences produces a characteristic sparsity pattern.

```{code-cell}

```

We now use the Boolean array that indicates where the boundary points lie in the grid.

```{code-cell}
spy(sparse(is_boundary);
    m=3,  color=:darkblue, 
    legend=:none,  title="Boundary points",
    xaxis=("column index", [0, n+2]), 
    yaxis=("row index", [0, m+2]))
```

In order to impose Dirichlet boundary conditions, we replace the boundary rows of the system by rows of the identity.

```{code-cell}
I_N = I(N)
idx = vec(is_boundary)
A[idx, :] .= I_N[idx, :];     # Dirichlet conditions
```

```{code-cell}
:tags: hide-input
spy(A;
    color=:blues,  m=3,
    title="System matrix with boundary conditions")
```

Finally, we must replace the rows in the vector $\mathbf{b}$ by the boundary values being assigned to the boundary points. Here, we let the boundary values be zero everywhere.

```{code-cell}
b[idx] .= 0;                 # Dirichlet values
```

Now we can solve for $\mathbf{u}$ and reinterpret it as the matrix-shaped $\mathbf{U}$, the solution on our grid.

```{code-cell}
u = A \ b
U = unvec(u)
```
``````

(demo-laplace-poisson-julia)=
``````{dropdown} @demo-laplace-poisson
:open:

First we define the problem on $[0,1]\times[0,2]$.

```{code-cell}
f = (x, y) -> -sin(3x * y - 4y) * (9y^2 + (3x - 4)^2);
g = (x, y) -> sin(3x * y - 4y);
xspan = [0, 1];
yspan = [0, 2];
```

Here is the finite-difference solution.

```{code-cell}
x, y, U = FNC.poissonfd(f, g, 40, xspan, 60, yspan);
```

```{code-cell}
surface(x, y, U';
    color=:viridis,
    title="Solution of Poisson's equation",
    xaxis=(L"x"),  yaxis=(L"y"),  zaxis=(L"u(x,y)"),
    right_margin=3Plots.mm,  camera=(70, 50))
```

The error is a smooth function of $x$ and $y$. It must be zero on the boundary; otherwise, we have implemented boundary conditions incorrectly.

```{code-cell}
error = [g(x, y) for x in x, y in y] - U;
M = maximum(abs, error)
contour(x, y, error';
    levels=17, 
    clims=(-M, M), color=:redsblues, 
    colorbar=:bottom,  aspect_ratio=1,
    title="Error", 
    xaxis=(L"x"),  yaxis=(L"y"),
    right_margin=7Plots.mm)
plot!([0, 1, 1, 0, 0], [0, 0, 2, 2, 0], l=(2, :black))
```
``````

### 13.4 @section-twodim-nonlinear

(demo-nonlinear2d-mems-julia)=
``````{dropdown} @demo-nonlinear2d-mems
:open:
All we need to define are $\phi$ from {eq}`nonlinpdepde` for the PDE, and a trivial zero function for the boundary condition.

```{code-cell}
λ = 1.5
ϕ = (X, Y, U, Ux, Uxx, Uy, Uyy) -> @. Uxx + Uyy - λ / (U + 1)^2;
g = (x, y) -> 0;
```

Here is the solution for $m=15$, $n=8$.

```{code-cell}
u = FNC.elliptic(ϕ, g, 15, [0, 2.5], 8, [0, 1]);
```

```{code-cell}
using Plots
x = range(0, 2.5, 100)
y = range(0, 1, 50)
U = [u(x, y) for x in x, y in y]
contourf(x, y, U';
    color=:blues,  l=0,
    aspect_ratio=1,
    xlabel=L"x",  ylabel=L"y",  zlabel=L"u(x,y)",
    title="Deflection of a MEMS membrane",
    right_margin=3Plots.mm)
```

In the absence of an exact solution, how can we be confident that the solution is accurate? First, the Levenberg iteration converged without issuing a warning, so we should feel confident that the discrete equations were solved. We can check the boundary values easily. For example,

```{code-cell}
x_test = range(0, 2.5, 100)
norm([u(x, 0) - g(x, 0) for x in x_test], Inf)
```

Assuming that we encoded the PDE correctly, the remaining source error is truncation from the discretization. We can estimate that by refining the grid a bit and seeing how much the numerical solution changes.

```{code-cell}
x_test = range(0, 2.5, 6)
y_test = range(0, 1, 6)
mtx_test, _ = FNC.tensorgrid(x_test, y_test)
mtx_test(u)
```

```{code-cell}
u = FNC.elliptic(ϕ, g, 25, [0, 2.5], 14, [0, 1]);
mtx_test(u)
```

The original solution seems to be accurate to about four digits.

``````

(demo-nonlinear-advdiff-julia)=
``````{dropdown} @demo-nonlinear-advdiff
:open:

```{code-cell}
ϕ = (X, Y, U, Ux, Uxx, Uy, Uyy) -> @. 1 - Ux - 2Uy + 0.05 * (Uxx + Uyy)
g = (x, y) -> 0
u = FNC.elliptic(ϕ, g, 32, [-1, 1], 32, [-1, 1]);
```

```{code-cell}
x = y = range(-1, 1, 80)
U = [u(x, y) for x in x, y in y]
contourf(x, y, U';
    color=:viridis, 
    aspect_ratio=1,
    xlabel=L"x",  ylabel=L"y",  zlabel=L"u(x,y)",
    title="Steady advection–diffusion")
```

``````

(demo-nonlinear2d-allencahn-julia)=
``````{dropdown} @demo-nonlinear2d-allencahn
:open:

The following defines the PDE and a nontrivial Dirichlet boundary condition for the square $[0,1]^2$.

```{code-cell}
ϕ = (X, Y, U, Ux, Uxx, Uy, Uyy) -> @. U * (1 - U^2) + 0.05 * (Uxx + Uyy)
g = (x, y) -> tanh(5 * (x + 2y - 1));
```

We solve the PDE and then plot the result.

```{code-cell}
u = FNC.elliptic(ϕ, g, 36, [0, 1], 36, [0, 1]);
```

```{code-cell}
x = y = range(0, 1, 80)
U = [u(x, y) for x in x, y in y]
contourf(x, y, U';
    color=:viridis, 
    aspect_ratio=1,
    xlabel=L"x",  ylabel=L"y",  zlabel=L"u(x,y)", 
    title="Steady Allen-Cahn",
    right_margin=3Plots.mm)
```
``````