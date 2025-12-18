---
numbering:
    enumerator: 13.2.%s
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.12
---
```{code-cell}
:tags: [remove-cell]
import Pkg; Pkg.activate("/Users/driscoll/Documents/GitHub/fnc")

using FNCFunctions

using Plots
default(
    titlefont=(11,"Helvetica"),
    guidefont=(11,"Helvetica"),
    linewidth = 2,
    markersize = 3,
    msa = 0,
    size=(500,320),
    label="",
    html_output_format = "svg"
)

using PrettyTables, LaTeXStrings, Printf
using LinearAlgebra

@ptconf backend = Val(:html) tf = tf_html_simple
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
| `u = vec(U)` | `U = reshape(u, m+1, n+1)` |


``````

::::{prf:example} Reshaping for grid functions
:label: demo-diffadv-vec


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

::::

In order to modularize our codes, we use @function-tensorgrid to define functions and values related to working with tensor-product grids. Its final output is to be discussed and used in @section-twodim-laplace.

``````{prf:algorithm} tensorgrid
:label: function-tensorgrid

```{literalinclude} chapter13.jl
:filename: tensorgrid.jl
:linenos: true
:language: julia
:start-after: # begin tensorgrid
:end-before: # end tensorgrid
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
m, n = (60, 25)
x, Dx, Dxx = FNC.diffper(m, [-1, 1])
y, Dy, Dyy = FNC.diffper(n, [-1, 1])
mtx, X, Y, unvec, _ = FNC.tensorgrid(x, y);
```

Here is the initial condition, evaluated on the grid.

```{code-cell}
u_init = (x, y) -> sin(4 * π * x) * exp(cos(π * y))
U₀ = mtx(u_init)
M = maximum(abs, U₀)
contour(x, y, U₀';
    color=:redsblues,  clims=(-M, M), 
    aspect_ratio=1,
    xaxis=("x", (-1, 1)),  yaxis=("y", (-1, 1)), 
    title="Initial condition" )
```

The following function computes the time derivative for the unknowns, which have a vector shape. The actual calculations, however, take place using the matrix shape.

```{code-cell}
function du_dt(u, α, t)
    U = unvec(u)
    Uxx = Dxx * U
    Uyy = U * Dyy'            # 2nd partials
    du_dt = α * (Uxx + Uyy)    # PDE
    return vec(du_dt)
end;
```

Since this problem is parabolic, a stiff time integrator is appropriate.

```{code-cell}
using OrdinaryDiffEq
IVP = ODEProblem(du_dt, vec(U₀), (0, 0.2), 0.1)
sol = solve(IVP, Rodas4P());
```

We can use the function `sol` defined above to get the solution at any time. Its output is the matrix $\mathbf{U}$ of values on the grid. An animation shows convergence of the solution toward a uniform value.

```{tip}
:class: dropdown
To plot the solution at any time, we use the same color scale as with the initial condition, so that the pictures are more easily compared.
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
mtx, X, Y, _ = FNC.tensorgrid(x, y)
U₀ = mtx( (x, y) -> (1 + y) * (1 - x)^4 * (1 + x)^2 * (1 - y^4) );
```

We define functions `extend` and `chop` to deal with the Dirichlet boundary conditions. 

```{code-cell}
chop = U -> U[2:m, 2:n]
z = zeros(1, n-1)
extend = W -> [zeros(m+1) [z; W; z] zeros(m+1)];
```

Next, we define the `pack` and `unpack` functions, using another call to @function-tensorgrid to get reshaping functions for the interior points. 

```{code-cell}
_, _, _, unvec, _ = FNC.tensorgrid(x[2:m], y[2:n])
pack = U -> vec(chop(U))
unpack = w -> extend(unvec(w));
```

Now we can define and solve the IVP using a stiff solver.

```{code-cell}
function dw_dt(w, ϵ, t)
    U = unpack(w)
    Ux, Uxx = Dx * U, Dxx * U
    Uyy = U * Dyy'
    dU_dt = @. 1 - Ux + ϵ * (Uxx + Uyy)
    return pack(dU_dt)
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
        xlabel=L"x",  ylabel=L"y",  zaxis=((0, 2.3), L"u(x,y)"),
        color=:blues,  alpha=0.66,  clims=(0, 2.3), colorbar=:none,
        title="Advection-diffusion",  dpi=150)
    contour!(x, y, U'; 
        levels=24, 
        aspect_ratio=1,  subplot=2, 
        xlabel=L"x",  ylabel=L"y",
        color=:blues,  clims=(0, 2.3),  colorbar=:none,
        title=@sprintf("t = %.2f", t))
end
mp4(anim, "figures/2d-advdiff.mp4");
```

![Advection-diffusion in 2d](figures/2d-advdiff.mp4)

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
m, n = 40, 40
x, Dx, Dxx = FNC.diffcheb(m, [-2, 2])
y, Dy, Dyy = FNC.diffcheb(n, [-2, 2])
mtx, X, Y, _ = FNC.tensorgrid(x, y)
U₀ = mtx( (x, y) -> (x + 0.2) * exp(-12 * (x^2 + y^2)) )
V₀ = zeros(size(U₀));
```

We need to define chopping and extension for the $u$ component. This looks the same as in @demo-diffadv-advdiff.

```{code-cell}
chop = U -> U[2:m, 2:n]
extend = U -> [zeros(m+1) [zeros(1, n-1); U; zeros(1, n-1)] zeros(m+1)];
```

While `vec` is the same for both the interior and full grids, the `unvec` operation is defined differently for them. 

```{code-cell}
_, _, _, unvec_v, _ = FNC.tensorgrid(x, y)
_, _, _, unvec_u, _ = FNC.tensorgrid(x[2:m], y[2:n])
N = (m-1) * (n-1)    # number of interior unknowns

pack = (U, V) -> [vec(chop(U)); vec(V)]
unpack = w -> (extend(unvec_u(w[1:N])), unvec_v(w[N+1:end]));
```

We can now define and solve the IVP. Since this problem is hyperbolic, a nonstiff integrator is faster than a stiff one.

```{code-cell}
function dw_dt(w, c, t)
    U, V = unpack(w)
    du_dt = V
    dv_dt = c^2 * (Dxx * U + U * Dyy')
    return pack(du_dt, dv_dt)
end

IVP = ODEProblem(dw_dt, pack(U₀, V₀), (0, 4.0), 1)
sol = solve(IVP, Tsit5())
U = t -> unpack(sol(t))[1];
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
