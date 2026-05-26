---
numbering:
  enumerator: 13.1.%s
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

(section-twodim-tensorprod)=

# Tensor-product discretizations

As you learned when starting double integration in vector calculus, the simplest extension of an interval to two dimensions is a rectangle. We will use a particular notation for rectangles.

```{index} ! tensor-product domain
```

```{prf:definition} Tensor-product domain
:label: definition-tensordomain
A {term}`tensor-product domain` is two dimensions is the rectangle

```{math}
:label: rectangleTP
  [a,b] \times [c,d] = \bigl\{ (x,y)\in\mathbb{R}^2 : a\le x \le b,\; c\le y \le d \bigr\}.
```

```{note}
The $\times$ notation in @rectangleTP is called a **tensor product**. 
```

The idea of the tensor product is that each variable independently varies over a fixed set. When the interval is the same in each dimension, we may write $[a,b]^2$.

The discretization of a two-dimensional tensor-product domain is straightforward.

````{prf:definition} Tensor-product grid
:label: definition-tensorgrid
Given discretizations of two intervals,

```{math}
:label: twointervals
a= x_0< x_1 < \cdots < x_m = b,  \qquad c = y_0 < y_1 < \cdots < y_n = d,
```

then a **tensor-product grid** on $[a,b]\times[c,d]$ is the set

```{math}
:label: rectangledisc
  \bigl\{ (x_i, y_j): i=0,\ldots,m,\; j=0,\ldots,n \bigr\}.
```
````

@fig-tensor-grid shows a tensor-product grid on a rectangle. The grid is constructed from discretizations of $x$ with $m=3$ and $y$ with $n=5$. Each grid point $(x_i, y_j)$ lies at the intersection of a vertical line through $x_i$ and a horizontal line through $y_j$.

```{figure} ../_static/tensor-grid.svg
:label: fig-tensor-grid
:alt: Tensor-product grid
:align: center
A tensor-product grid (dark gray) on a rectangle constructed from discretizations of $x$ with $m=3$ (blue) and $y$ with $n=5$ (red).
```

```{index} ! tensor-product grid
```

## Functions on grids

The double indexing of the grid set {eq}`rectangledisc` implies an irresistible connection to matrices. Corresponding to any function $f(x,y)$ defined on the rectangle is an $(m+1)\times(n+1)$ matrix $\mathbf{F}$ defined by collecting the values of $f$ at the points in the grid. This transformation of a function to a matrix is so important that we give it a formal name

````{prf:definition} Function-to-matrix map
:label: definition-mtx
The function mtx maps a function $f(x,y)$ to an $(m+1)\times (n+1)$  matrix $\mathbf{F}$ by

```{math}
:label: fun2mtx
\mathbf{F} = \mtx(f) = \Bigl[f(x_i,y_j)\Bigr]_{\substack{i=0,\ldots,m\\j=0,\ldots,n}},
```

where the evaluations are on the grid {eq}`rectangledisc`.
````

```{warning}
Traditionally, the first dimension of a matrix varies in the _vertical_ direction, while the first space coordinate $x$ varies _horizontally_. This clash has long been a source of confusion when coding, and each language or package has its own conventions around it. In this book, the first dimension of a matrix always corresponds to the first coordinate of the function. Sometimes, this forces us to use a transpose for making plots.
```

::::{prf:example}
:label: example-tensorprod-smallgrid
Let the interval $[0,2]$ be divided into $m=4$ equally sized pieces, and let $[1,3]$ be discretized in $n=2$ equal pieces. Then the grid in the rectangle $[0,2]\times[1,3]$ is given by all points $(i/2,1+j)$ for all choices $i=0,1,2,3,4$ and $j=0,1,2$. If $f(x,y)=\sin(\pi xy)$, then

\begin{equation*}
  \mtx(f) =
    \begin{bmatrix}
    \sin(\pi\cdot 0\cdot 1) & \sin(\pi\cdot0\cdot 2) & \sin(\pi\cdot0\cdot 3) \\[1mm]
    \sin\left(\pi\cdot\tfrac{1}{2} \cdot 1 \right) & \sin\left(\pi\cdot\tfrac{1}{2} \cdot 2 \right) & \sin\left(\pi\cdot\tfrac{1}{2} \cdot 3 \right) \\[1mm]
    \sin\left(\pi \cdot 1 \cdot 1 \right) & \sin\left(\pi \cdot 1 \cdot 2 \right) & \sin\left(\pi \cdot 1 \cdot 3 \right) \\[1mm]
    \sin\left(\pi\cdot \tfrac{3}{2} \cdot 1 \right) & \sin\left(\pi\cdot\tfrac{3}{2} \cdot 2 \right) & \sin\left(\pi\cdot\tfrac{3}{2} \cdot 3 \right) \\[1mm]
    \sin\left(\pi \cdot 2 \cdot 1 \right) & \sin\left(\pi \cdot 2 \cdot 2 \right) & \sin\left(\pi \cdot 2 \cdot 3 \right)
    \end{bmatrix}
    = \begin{bmatrix}
    0 & 0 & 0 \\ 1 & 0 & -1 \\ 0 & 0 & 0 \\ -1 & 0 & 1 \\ 0 & 0 & 0
    \end{bmatrix}.
\end{equation*}

::::

::::{prf:example} Functions on 2D grids
:label: demo-tensorprod-gridfun

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

::::

## Parameterized surfaces

We are not limited to rectangles by tensor products. Many regions and surfaces may be parameterized by means of $x(u,v)$, $y(u,v)$, and $z(u,v)$, where $u$ and $v$ lie in a rectangle. Such "logically rectangular" surfaces include the unit disk,

```{math}
:label: unitdiskparam
\left\{
\begin{aligned}
x &= u \cos v, \\
y &= u \sin v, \\
\end{aligned}
\right.
\qquad \qquad
\left.
\begin{aligned}
0 & \le u < 1, \\
0 & \le v \le 2\pi,
\end{aligned}
\right.
```

and the unit sphere,

```{math}
:label: spheredomain
\left\{
\begin{aligned}
x &= \cos u \sin v,\\
y &= \sin u \sin v,\\
z &= \cos v,
\end{aligned}
\right.
\qquad \qquad
  \left.
\begin{aligned}
0 & \le u < 2\pi, \\
0 &\le v \le \pi.
\end{aligned}
\right.
```

::::{prf:example} Functions on parameterized surfaces
:label: demo-tensorprod-disksphere

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


::::

## Partial derivatives

In order to solve boundary-value problems in one dimension by collocation, we replaced an unknown function $u(x)$ by a vector of its values at selected nodes and discretized the derivatives in the equation using differentiation matrices. We use the same ideas in the 2D case: we represent a function by its values on a grid, and multiplication by differentiation matrices to construct discrete analogs of the partial derivatives $\frac{\partial u}{\partial x}$ and $\frac{\partial u}{\partial y}$.

```{index} differentiation matrix
```

Consider first $\frac{\partial u}{\partial x}$. In the definition of this partial derivative, the independent variable $y$ is held constant. Note that $y$ is constant within each column of $\mathbf{U} = \mtx(u)$. Thus, we may regard a single column $\mathbf{u}_j$ as a discretized function of $x$ and, as usual, left-multiply by a differentiation matrix $\mathbf{D}_x$ such as {eq}`diffmat12b`. We need to do this for each column of $\mathbf{U}$ by $\mathbf{D}_x$, which is accomplished by $\mathbf{D}_x \mathbf{U}$. Altogether,

```{math}
:label: partfpartx
\mtx\left( \frac{\partial u}{\partial x} \right) \approx \mathbf{D}_x \, \mtx(u).
```

This relation is not an equality, because the left-hand side is a discretization of the exact partial derivative, while the right-hand side is a
finite-difference approximation. Yet it is a natural analog for partial differentiation when we are given not $u(x,y)$ but only the grid value matrix $\mathbf{U}.$

Now we tackle $\frac{\partial u}{\partial y}$. Here the inactive coordinate $x$ is held fixed within each *row* of $\mathbf{U}$. However, if we transpose $\mathbf{U}$, then the roles of rows and columns are swapped, and now $y$ varies independently down each column. This is analogous to the situation for the $x$-derivative, so we left-multiply by a finite-difference matrix $\mathbf{D}_y$, and then transpose the entire result to restore the roles of $x$ and $y$ in the grid. Fortunately, linear algebra allows us to express the sequence transpose–left-multiply–transpose more compactly:

```{math}
:label: partfparty
\mtx\left( \frac{\partial u}{\partial y} \right) \approx \Bigl(\mathbf{D}_y \mathbf{U}^T\Bigr)^T = \mtx(u)\, \mathbf{D}_y^T.
```

Keep in mind that the differentiation matrix $\mathbf{D}_x$ is based on the discretization $x_0,\ldots,x_m$, and as such it must be $(m+1)\times (m+1)$. On the other hand, $\mathbf{D}_y$ is based on $y_0,\ldots,y_n$ and is $(n+1)\times (n+1)$. This is exactly what is needed dimensionally to make the products in {eq}`partfpartx` and {eq}`partfparty` consistent. More subtly, if the differentiation is based on equispaced grids in each variable, the value of $h$ in a formula such as {eq}`centerFD12` will be different for $\mathbf{D}_x$ and $\mathbf{D}_y$.

::::{prf:example} Numerical partial derivatives
:label: demo-tensorprod-diff

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

::::

## Exercises

``````{exercise}
:label: problem-tensorprod-surfcontourplot
⌨ In each part, make side-by-side surface and contour plots of the given function over the given domain.

**(a)** $f(x,y) = 2y + e^{x-y}$, $\quad[0,2]\times[-1,1]$

**(b)** $f(x,y) = \tanh[5(x+xy-y^3)]$, $\quad[-2,2]\times[-1,1]$

**(c)** $f(x,y) = \exp \bigl[-6(x^2+y^2-1)^2 \bigr]$, $\quad[-2,2]\times[-2,2]$

``````

``````{exercise}
:label: problem-tensorprod-derivatives

⌨ For each function in @problem-tensorprod-surfcontourplot, make side-by-side surface plots of $f_x$ and $f_y$ using Chebyshev spectral differentiation.
``````

``````{exercise}
:label: problem-tensorprod-mixed

⌨ For each function in @problem-tensorprod-surfcontourplot, make a contour plot of the mixed derivative $f_{xy}$ using Chebyshev spectral differentiation.
``````

``````{exercise}
:label: problem-tensorprod-polar

⌨ In each case, make a plot of the function given in polar or Cartesian coordinates over the unit disk.

**(a)** $f(r,\theta) = r^2 - 2r\cos \theta$

**(b)** $f(r,\theta) = e^{-10r^2}$

**(c)** $f(x,y) = xy - 2 \sin (x)$
``````

``````{exercise}
:label: problem-tensorprod-sphere

⌨ Plot $f(x,y,z)=x y - x z - y z$ as a function on the unit sphere.
%(Use `aspect_ratio=1` in a plot call to get equal aspect ratios for the axes.)
``````

``````{exercise}
:label: problem-tensorprod-cylinder

⌨ Plot $f(x,y,z)=x y - x z - y z$ as a function on the cylinder $r=1$ for $-1\le z \le 2$.
%(Use `aspect_ratio=1` in a plot call to get equal aspect ratios for the axes.)
``````
