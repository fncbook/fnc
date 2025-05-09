---
numbering:
    enumerator: 13.2.%s
---
(section-twodim-diffadv)=
# Two-dimensional diffusion and advection

We next describe how to apply the method of lines to PDEs of the form

:::{math}
:label: pde2d
  u_t = \phi(u,u_x,u_y,u_{xx},u_{xy},u_{yy}), \quad (x,y)\in [a,b]\times [c,d].
:::

The PDE may be of either parabolic or hyperbolic type, with the primary difference being potential restrictions on the time step size. To keep descriptions and implementations relatively simple, we will consider only periodic conditions or Dirichlet boundary conditions.

As described in {numref}`section-twodim-tensorprod`, the rectangular domain is discretized by a grid $(x_i,y_j)$ for $i=0,\ldots,m$ and $j=0,\ldots,n$. The solution is semidiscretized as a matrix $\mathbf{U}(t)$ such that $U_{ij}$ is the approximate solution at $(x_i,y_j,t)$. Terms involving the spatial derivatives of $u$ are readily replaced by discrete counterparts, as shown in @tab-2s-derivatives.

```{table} Discrete replacements for 2D derivatives
:label: tab-2d-derivatives
:align: center

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

(definition-diffadv-vec)=
::::{prf:definition} vec and unvec operations
Let $\mathbf{A}$ be an $m\times n$ matrix. Define the **vec** function as stacking the columns of $\mathbf{A}$ into a vector, i.e.,

:::{math}
:label: vecdef
\operatorname{vec}(\mathbf{A}) =
\begin{bmatrix}
A_{11} \\ \vdots \\ A_{m1}  \\ \vdots  \\ A_{1n} \\ \vdots \\ A_{m n}
\end{bmatrix}.
:::

Let $\mathbf{z}$ be a vector of length $m n$. Define the **unvec** function as the inverse of vec:

:::{math}
:label: unvecdef
\operatorname{unvec}(\mathbf{z}) = \begin{bmatrix}
  z_1 & z_{m+1} & \cdots & z_{m(n-1)+1} \\
  z_2 & z_{m+2} & \cdots & z_{m(n-1)+2} \\
  \vdots & \vdots & & \vdots \\
  z_m & z_{2m} & \cdots & z_{m n} \\
\end{bmatrix}.
:::
::::

Suppose $\mathbf{U} = \operatorname{mtx}(u)$ is the matrix of unknowns. @tab-vec-unvec shows how to convert between $\mathbf{U}$ and $\mathbf{u} = \operatorname{vec}(\mathbf{U})$.

``````{table} vec and unvec operations
:label: tab-vec-unvec
:align: center

`````{tab-set}
````{tab-item} Julia
:sync: julia

|vec | unvec |
|:---|:---|
| `u = vec(U)` | `U = reshape(u, m+1, n+1)` |

````

````{tab-item} MATLAB
:sync: matlab

|vec | unvec |
|:---|:---|
| `u = U(:)` | `U = reshape(u, m+1, n+1)` |

````

````{tab-item} Python
:sync: python

|vec | unvec |
|:---|:---|
| `u = U.T.flatten()` | `U = np.reshape(u, (n+1, m+1)).T` |


```{dropdown} About the code
The incorporation of transposes is because NumPy uses row-major order, while MATLAB and Julia use column-major order. If performance were a concern, we would reverse our convention to avoid the transposes. We also use `flatten` rather than `ravel` to ensure making copies rather than views of the data and avoiding subtle bugs.
```

`````

``````

(demo-diffadv-vec)=
::::{prf:example} Reshaping for grid functions

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-diffadv-vec-julia
:::
````

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-diffadv-vec-matlab
:::
````

````{tab-item} Python
:sync: python
:::{embed} #demo-diffadv-vec-python
:::
````
`````

::::

In order to modularize our codes, we use @function-tensorgrid to define functions and values related to working with tensor-product grids. Its final output is to be discussed and used in @section-twodim-laplace.

(function-tensorgrid)=
``````{prf:algorithm} tensorgrid
`````{tab-set} 
````{tab-item} Julia
:sync: julia
:::{embed} #function-tensorgrid-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #function-tensorgrid-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #function-tensorgrid-python
:::
````
`````
``````

## Periodic end conditions

If the boundary conditions are periodic, then the unknowns in the method of lines are the elements of the matrix $\mathbf{U}(t)$ representing grid values of the numerical solution. For the purposes of an IVP solution, this matrix is equivalent to the vector $\mathbf{u}(t)$ defined as $\mathbf{u}=\operatorname{vec}(\mathbf{U})$.

```{index} heat equation
```

(demo-diffadv-heat)=
::::{prf:example} Heat equation in 2D
We will solve a 2D heat equation, $u_t = 0.1(u_{xx} + u_{yy})$, on the square $[-1,1]\times[-1,1]$, with periodic behavior in both directions. The initial condition is 

```{math}
:label: twodim-heatinit
u(x,y,0) = \sin(4\pi x) \exp[\cos(\pi y)], 
```

which is also periodic on the rectangle.

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-diffadv-heat-julia
:::
````

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-diffadv-heat-matlab
:::
````

````{tab-item} Python
:sync: python
:::{embed} #demo-diffadv-heat-python
:::
````
`````
::::

## Dirichlet conditions

```{index} boundary conditions; numerical implementation of
```

In {numref}`section-diffusion-boundaries` we coped with boundary conditions by removing the boundary values from the vector of unknowns being solved in the semidiscretized ODE. Each evaluation of the time derivative required us to extend the values to include the boundaries before applying differentiation matrices in space, then remove them from the time derivative vector. 

We proceed similarly here, defining `chop` and `extend` functions that convert between the full grid and the inner grid of interior values. Mathematically speaking, the chopping operation is defined by

:::{math}
:label: tensorprod-chop
\operatorname{chop}(\mathbf{U}) = \mathbf{E}_x \mathbf{U} \mathbf{E}_y^T,
:::

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

The extension operator is a bit more awkward to write out. It stars with appending rows and columns of zeros around the border of a matrix $\mathbf{V}$ of interior values:

```{math}
:label: tensorprod-extend
\tilde{\mathbf{U}} = \mathbf{E}_x^T \mathbf{V} \mathbf{E}_y
```

We can then modify the new zero values to reflect the boundary conditions, via

```{math}
:label: tensorprod-extend2
\mathbf{U} = \operatorname{extend}(\mathbf{V}) = \tilde{\mathbf{U}} + \mathbf{G},
```

except that we have changes in the shape as well as boundary conditions to consider. Our strategy is to define two functions, `pack` and `unpack`, that convert between the vector of unknown values only and the full matrix of grid values.

Suppose we are given a matrix $\mathbf{U}$ that represents the solution on an $(m+1)\times (n+1)$ grid, including boundary values. Then we define


 All that remains, then, are the interior values, which are converted into a vector by the vec operator.

For the inverse transformation, suppose we are given a vector $\mathbf{w}$ of interior solution values. Then we define

$$
\operatorname{unpack}(\mathbf{w}) = .
$$

This operator reshapes the vector to a grid of interior values, then appends one extra zero row and column on each side of the grid.[^jacobian]

[^jacobian]: You might wonder why we use linear algebra to define the extension and deletion of boundary values rather than directly accessing row and column indices in the grid function. The algebraic expressions make it easier to express the Jacobian of a nonlinear The linear algebra approach allows `DifferentialEquations` to compute the Jacobian matrix of the implicit IVP solver quickly using *automatic differentiation* tools, greatly speeding up the solution process. Since the matrices in our expressions are sparse, multiplications by them do not affect running time perceptibly.

Now suppose the ODE unknowns for the interior solution values are in the vector $\mathbf{w}(t)$. When we form $\operatorname{unpack}(\mathbf{w})$, we reinterpret the values on the tensor-product grid and then extend these values to zero around the boundary. If the boundary values are given as $g(x,y)$, then $g$ has to be evaluated at the boundary nodes of the grid and inserted into the grid function matrix. Then the grid values are used to compute partial derivatives in $x$ and $y$, the discrete form of the PDE is evaluated, and we pack the result as the computation of $\mathbf{w}'$.

```{index} advection-diffusion equation
```

(demo-diffadv-advdiff)=
::::{prf:example} Advection-diffusion equation in 2D

We will solve an advection-diffusion problem, $u_t + u_x = 1 + \epsilon(u_{xx} + u_{yy})$, where $u=0$ on the boundary of the square $[-1,1]^2$. The outline of our approach is based on {numref}`Function {number} <function-parabolic>` for parabolic PDEs in one space dimension.

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-diffadv-advdiff-julia
:::
````

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-diffadv-advdiff-matlab
:::
````

````{tab-item} Python
:sync: python
:::{embed} #demo-diffadv-advdiff-python
:::
````
`````
::::


```{index} wave equation
```

The wave equation introduces a little additional complexity. First, we write the 2D wave equation $u_{tt}=c^2(u_{xx}+u_{yy})$ in first-order form as

:::{math}
:label: wave2dfirst
\begin{split}
    u_t &= v, \\
    v_t &= c^2(u_{xx}+u_{yy}).
\end{split}
:::

Now the grid unknowns are a pair of matrices $\mathbf{U}(t)$ and $\mathbf{V}(t)$. Typical boundary conditions would prescribe $u$ on all of the boundary and let $v$ be unspecified. Since the boundary values of $\mathbf{U}$ are prescribed, those values are omitted from the semidiscretization IVP, while all of $\mathbf{V}$ is included. All of these unknowns need to be packed into and unpacked from a single vector $\mathbf{w}(t)$ for the IVP solver.

(demo-diffadv-wave)=
::::{prf:example} Wave equation in 2D

We solve the wave equation with $c=1$ on the square $[-2,2]\times[-2,2]$, where $u=0$ on the boundary. 

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-diffadv-wave-julia
:::
````

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-diffadv-wave-matlab
:::
````

````{tab-item} Python
:sync: python
:::{embed} #demo-diffadv-wave-python
:::
````
`````

::::

## Exercises

``````{exercise}
⌨  For the given $u(x,y)$, make a plot of the given quantity on the square $[-2,2]^2$ using appropriate differentiation matrices.

**(a)** $u(x,y) = \exp(x-y^2)$; plot $u_{xx}+u_{yy}$

**(b)** $u(x,y) =\cos (\pi x)+\sin (\pi y)$; plot $u_x+u_y$

**(c)** $u(x,y) =\exp(-x^2-4y^2)$; plot $x u_y$
``````

``````{exercise}
⌨ Following {numref}`Demo %s <demo-diffadv-heat>` as a model, solve the Allen–Cahn equation $u_t=u(1-u^2)+0.001(u_{xx}+u_{yy})$ on the square $[-1,1]^2$ with periodic conditions, taking $u(x,y,0)=\sin(\pi x)\cos(2\pi y)$. Use $m=n=60$ to solve up to $t=4$, and make an animation of the result.
``````

``````{exercise}
⌨ Following {numref}`Demo %s <demo-diffadv-advdiff>` as a model, solve $u_t=y u_x-u_y+0.03(u_{xx}+u_{yy})$ on the square $[-1,1]^2$, with $u(x,y,0)=(1-x^2)(1-y^2)$ and homogeneous Dirichlet boundary conditions. Use $m=n=40$ to solve up to $t=2$, and make an animation of the result.
``````

``````{exercise}
⌨ Following {numref}`Demo %s <demo-diffadv-wave>` as a model, solve $u_{tt}=u_{xx}+u_{yy}+\cos(7t)$ on the square $[-1,1]^2$, with $u(x,y,0)=x(1-x^6)(1-y^2)$, $u_t(x,y,0)=0$, subject to homogeneous Dirichlet boundary conditions. Take $m=n=60$ to solve between $t=0$ and $t=12$, and make an animation of the result.
``````

``````{exercise}
From Maxwell's equations we can find a way to convert the wave equation to a first-order form that, unlike {eq}`wave2dfirst`, uses only first-order derivatives in space:

:::{math}
:label: wave2dTM
\begin{split}
u_t &= c^2(v_y - w_x),\\
v_t &= u_y, \\
w_t &= -u_x,
\end{split}
:::

subject to $u=0$ on the boundary.

**(a)** ✍ Show that a solution of {eq}`wave2dTM` satisfies $u_t=c^2(u_{xx}+u_{yy})$.

**(b)** ⌨ Solve {eq}`wave2dTM` with $c=2$ in the rectangle $x\in[-3,3]$, $y\in[-1,1]$, $u(x,y,0) = \exp(x-x^2)(9-x^2)(1-y^2)$, and $v=w=0$ at $t=0$. Use $m=50$ for $x$ and $n=25$ for $y$, solve for $0\le t \le 6$, and make an animation of the solution.
``````
