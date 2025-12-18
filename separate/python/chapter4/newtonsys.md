---
numbering:
  enumerator: 4.5.%s
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

(section-nonlineqn-newtonsys)=

# Newton for nonlinear systems

The rootfinding problem becomes much more difficult when multiple variables and equations are involved.

```{index} ! rootfinding problem; multidimensional
```

````{prf:definition} Multidimensional rootfinding problem
:label: definition-vectorroot
Given a continuous vector-valued function $\mathbf{f}$ mapping from $\mathbb{R}^n$ into $\mathbb{R}^n$, find a vector $\mathbf{r}$ such that

```{math}
:label: rootvector
\begin{split}
  f_1(r_1,\dots,r_n) &= 0,\\
  f_2(r_1,\dots,r_n) &= 0,\\
  &\vdots\\
  f_n(r_1,\dots,r_n) &= 0.
\end{split}
```
````

Particular problems are often posed using scalar variables and equations.

::::{prf:example}
:label: example-newtonsys-predprey
The steady state of interactions between the population $w(t)$ of a predator species and the population $h(t)$ of a prey species might be modeled as

$$
\begin{split}
ah - b h w &= 0, \\
-cw + d w h &= 0,
\end{split}
$$

for positive parameters $a,b,c,d$. To cast this in the form of {eq}`rootvector`, we could define $\mathbf{x}=[h,w]$, $f_1(x_1,x_2) = ax_1 - bx_1x_2$, and $f_2(x_1,x_2)= -c x_2 + d x_1 x_2$.
::::

While the equations of {numref}`Example {number} <example-newtonsys-predprey>` are easy to solve by hand, in practice even establishing the existence and uniqueness of solutions for any particular system is typically quite difficult.

## Linear model

To extend rootfinding methods to systems, we will keep to the basic philosophy of constructing easily managed models of the exact function. As usual, the starting point is a linear model. We first need to define what it means to take a derivative of a vector-valued function of a vector variable.

```{index} ! derivative; vector-valued function
```

```{index} ! Jacobian matrix
```

::::{prf:definition} Jacobian matrix
:label: definition-jacobian
The {term}`Jacobian matrix` of $\mathbf{f}(\mathbf{x})$, where $\mathbf{f}$ is $m$-dimensional and $\mathbf{x}$ is $n$-dimensional, is the $m\times n$ matrix

```{math}
:label: jacobian
\mathbf{J}(\mathbf{x}) =
  \begin{bmatrix}
    \rule[2mm]{0pt}{1em}\frac{\partial f_1}{\partial x_1} & \frac{\partial f_1}{\partial x_2} & \cdots & \frac{\partial f_1}{\partial x_n}\\[2mm]
    \frac{\partial f_2}{\partial x_1} & \frac{\partial f_2}{\partial x_2} & \cdots & \frac{\partial f_2}{\partial x_n}\\[1mm]
    \vdots & \vdots & & \vdots\\[1mm]
    \rule[-3mm]{0pt}{1em} \frac{\partial f_n}{\partial x_1} & \frac{\partial f_n}{\partial x_2} & \cdots & \frac{\partial f_n}{\partial x_n}
  \end{bmatrix} = \left[ \frac{\partial f_i}{\partial x_j} \right]_{\,i=1,\ldots,m,\, j=1,\ldots,n.}
```

::::

(example-nonlinsystem)=

````{prf:example}
Let
  
```{math}
\begin{split}
    f_1(x_1,x_2,x_3) &= -x_1\cos(x_2) - 1\\
    f_2(x_1,x_2,x_3) &= x_1x_2 + x_3\\
    f_3(x_1,x_2,x_3) &= e^{-x_3}\sin(x_1+x_2) + x_1^2 - x_2^2.
\end{split}
```

Then
  
```{math}
    \mathbf{J}(x) =
    \begin{bmatrix}
       -\cos(x_2) & x_1 \sin(x_2) & 0\\
      x_2 & x_1 & 1\\
       e^{-x_3}\cos(x_1+x_2)+2x_1 & e^{-x_3}\cos(x_1+x_2)-2x_2 &
       -e^{-x_3}\sin(x_1+x_2)
    \end{bmatrix}.
```

If we were to start writing out the terms in {eq}`multitaylor`, we would begin with
  
```{math}
\begin{split}
    f_1(x_1+h_1,x_2+h_2,x_3+h_3) &= -x_1\cos(x_2)-1 -\cos(x_2)h_1 +
    x_1\sin(x_2)h_2 + O\bigl(\| \mathbf{h} \|^2\bigr) \\
    f_2(x_1+h_1,x_2+h_2,x_3+h_3) &= x_1x_2 + x_3 + x_2h_1 +x_1h_2 +
    h_3 + O\bigl(\| \mathbf{h} \|^2\bigr),
  \end{split}
```

and so on.
````

A multidimensional Taylor series begins with the linear approximation

```{math}
:label: multitaylor
\mathbf{f}(\mathbf{x}+\mathbf{h}) = \mathbf{f}(\mathbf{x}) + \mathbf{J}(\mathbf{x})\mathbf{h} + O(\| \mathbf{h} \|^2),
```

where $\mathbf{J}$ is the Jacobian matrix of $\mathbf{f}$ at $\mathbf{x}$, and $\mathbf{h}$ is a small perturbation from $\mathbf{x}$.  

The terms $\mathbf{f}(\mathbf{x})+\mathbf{J}(\mathbf{x})\mathbf{h}$ in {eq}`multitaylor` represent the linear part of $\mathbf{f}$ near $\mathbf{x}$, while the $O(\| \mathbf{h} \|^2)$ term represents the omitted higher-order terms in the Taylor series. If $\mathbf{f}$ is actually linear, i.e., $\mathbf{f}(\mathbf{x})=\mathbf{A}\mathbf{x}-\mathbf{b}$, then the Jacobian matrix is the constant matrix $\mathbf{A}$ and the higher-order terms in {eq}`multitaylor` disappear.

:::{note}
Because of the Jacobian's role in {eq}`multitaylor`, we may write $\mathbf{J}(\mathbf{x})$ as $\mathbf{f}{\,}'(\mathbf{x})$. Like any derivative, it is a function of the independent variable $\mathbf{x}$.
:::

## The multidimensional Newton iteration

```{index} Newton's method
```

With a method in hand for constructing a linear model for the vector system $\mathbf{f}(\mathbf{x})$, we can generalize Newton's method. Specifically, at a root estimate $\mathbf{x}_k$, we set $\mathbf{h} = \mathbf{x}-\mathbf{x}_k$ in {eq}`multitaylor` and get

```{math}
\mathbf{f}(\mathbf{x}) \approx \mathbf{q}(\mathbf{x})  = \mathbf{f}(\mathbf{x}_k) + \mathbf{J}(\mathbf{x}_k)(\mathbf{x}-\mathbf{x}_k).
```

We define the next iteration value $\mathbf{x}_{k+1}$ by requiring $\mathbf{q}(\mathbf{x}_{k+1})=\boldsymbol{0}$,

```{math}
\begin{split}
  \boldsymbol{0} &=  \mathbf{f}(\mathbf{x}_k) + \mathbf{J}(\mathbf{x}_k)(\mathbf{x}_{k+1}-\mathbf{x}_k),\\
\end{split}
```

which can be rearranged into

```{math}
:label: newtonsys
\mathbf{x}_{k+1} = \mathbf{x}_k - \bigl[\mathbf{J}(\mathbf{x}_k)\bigr]^{-1} \mathbf{f}(\mathbf{x}_k).
```

Note that $\mathbf{J}^{-1}\mathbf{f}$ now plays the role that $f/f'$ had in the scalar case; in fact, the two are the same in one dimension. In computational practice, however, we don't compute matrix inverses.

```{index} ! Newton's method; multidimensional
```

::::{prf:algorithm} Multidimensional Newton's method
:label: definition-newtonsmethodsys
Given $\mathbf{f}$ and a starting value $\mathbf{x}_1$, for each $k=1,2,3,\ldots$

1. Compute $\mathbf{y}_k = \mathbf{f}(\mathbf{x}_k)$ and $\mathbf{A}_k=\mathbf{f\,}'(\mathbf{x}_k)$.
2. Solve the linear system $\mathbf{A}_k\mathbf{s}_k = -\mathbf{y}_k$ for the **Newton step** $\mathbf{s}_k$.
3. Let $\mathbf{x}_{k+1} = \mathbf{x}_k + \mathbf{s}_k$.
::::

An extension of our series analysis of the scalar Newton's method shows that the vector version is also quadratically convergent in any vector norm, under suitable circumstances and when the iteration converges at all.

## Implementation

```{index} Julia; \\
```

An implementation of Newton's method for systems is given in {numref}`Function {number} <function-newtonsys>`. Other than computing the Newton step using backslash and taking vector magnitudes with `norm`, {numref}`Function {number} <function-newtonsys>` is virtually identical to the scalar version {numref}`Function {number} <function-newton>` presented earlier.

``````{prf:algorithm} newtonsys
:label: function-newtonsys

```{literalinclude} chapter04.py
:filename: newtonsys.py
:start-line: 68
:end-line: 96
:language: python
:linenos: true
```
````{admonition} About the code
:class: dropdown
The output of {numref}`Function {number} <function-newtonsys>` is a vector of vectors representing the entire history of root estimates. Since these should be in floating point, the starting value is converted with `float` before the iteration starts.
````
``````

::::{prf:example} Convergence of Newton's method for systems
:label: demo-newtonsys-converge

A system of nonlinear equations is defined by its residual and Jacobian.

```{code-cell}
def func(x):
    return array([
        exp(x[1] - x[0]) - 2, 
        x[0] * x[1] + x[2], 
        x[1] * x[2] + x[0]**2 - x[1]
    ])

def jac(x):
    return array([
            [-exp(x[1] - x[0]), exp(x[1] - x[0]), 0],
            [x[1], x[0], 1],
            [2 * x[0], x[2] - 1, x[1]],
    ])
```

Our initial guess at a root is the origin. 

```{code-cell}
x1 = zeros(3)
x = FNC.newtonsys(func, jac, x1)
print(x)
```

The output has one row per iteration, so the last row contains the final Newton estimate. Let's compute its residual.

```{code-cell}
r = x[-1]
f = func(r)
print("final residual:", f)
```

Let's check the convergence rate:

```{code-cell}
logerr = [log(norm(x[k] - r)) for k in range(x.shape[0] - 1)]
for k in range(len(logerr) - 1):
    print(logerr[k+1] / logerr[k])
```

The ratio is apparently converging toward 2, as expected for quadratic convergence.

::::

## Exercises

``````{exercise}
:label: problem-newtonsys-byhand
✍ Suppose that

```{math}
:numbered: false
\mathbf{f}(\mathbf{x}) =
\begin{bmatrix}
x_1x_2 + x_2^2 - 1 \\[1mm] x_1 x_2^3 + x_1^2 x_2^2 + 1
\end{bmatrix}.
```

Let $\mathbf{x}_1=[-2,1]^T$. Use Newton's method to find $\mathbf{x}_2$.
``````

``````{exercise}
:label: problem-newtonsys-linear
✍ Suppose that $\mathbf{f}(\mathbf{x}) = \mathbf{A}\mathbf{x} - \mathbf{b}$ for a constant $n\times n$ matrix $\mathbf{A}$ and constant $n\times 1$ vector $\mathbf{b}$. Show that Newton's method converges to the exact root in one iteration.
``````

``````{exercise}
:label: problem-newtonsys-spherepotential
Two curves in the $(u,v)$ plane are defined implicitly by the equations $u\log u + v \log v = -0.3$ and $u^4 + v^2 = 1$.

**(a)** ✍ Write the intersection of these curves in the form $\mathbf{f}(\mathbf{x}) = \boldsymbol{0}$ for two-dimensional $\mathbf{f}$ and $\mathbf{x}$.

**(b)** ✍ Find the Jacobian matrix of $\mathbf{f}$.

**(c)** ⌨ Use {numref}`Function {number} <function-newtonsys>` to find an intersection point starting from $u=1$, $v=0.1$.

**(d)** ⌨ Use {numref}`Function {number} <function-newtonsys>` to find an intersection point starting from $u=0.1$, $v=1$.

``````

``````{exercise}
:label: problem-newtonsys-orbitintersect
Two elliptical orbits $(x_1(t),y_1(t))$ and $(x_2(t),y_2(t))$ are described by the equations

```{math}
:numbered: false
\begin{bmatrix}
x_1(t) \\ y_1(t)
\end{bmatrix}
=
\begin{bmatrix}
-5+10\cos(t) \\ 6\sin(t)
\end{bmatrix}, \qquad
\begin{bmatrix}
x_2(t)\\y_2(t)
\end{bmatrix} =
\begin{bmatrix}
8\cos(t) \\ 3+12\sin(t)
\end{bmatrix},
```

where $t$ represents time.

**(a)** ⌨ Make a plot of the two orbits with the following code:

``` python
import numpy as np
import matplotlib.pyplot as plt
t = np.linspace(0, 2*np.pi, 100)
x1 = -5 + 10*np.cos(t);   y1 = 6*np.sin(t)
plt.plot(x1, y1, label='Orbit 1')
x2 = 8*np.cos(t);   y2 = 3 + 12*np.sin(t)
plt.plot(x2, y2, label='Orbit 2')
plt.axis('equal'); plt.grid(True); plt.legend()
plt.show()
```

**(b)** ✍ Write out a $2\times 2$ nonlinear system of equations that describes an intersection of these orbits. (Note: An intersection is not the same as a collision—they don't have to occupy the same point at the same time.)

**(c)** ✍ Write out the Jacobian matrix of this nonlinear system.

**(d)** ⌨ Use {numref}`Function {number} <function-newtonsys>` to find all of the unique intersections. Add them to the plot from part (a).

``````

``````{exercise}
:label: problem-newtonsys-ellipsemin
⌨  Suppose one wants to find the points on the ellipsoid $x^2/25 + y^2/16 + z^2/9 = 1$ that are closest to and farthest from the point $(5,4,3)$. The method of Lagrange multipliers implies that any such point satisfies

```{math}
:numbered: false
\begin{split}
x-5 &= \frac{\lambda x}{25}, \\[1mm]
y-4 &= \frac{\lambda y}{16}, \\[1mm]
z-3 &= \frac{\lambda z}{9}, \\[1mm]
1 &=  \frac{1}{25}x^2 + \frac{1}{16}y^2 + \frac{1}{9}z^2
\end{split}
```

for an unknown value of $\lambda$.

**(a)** Write out this system in the form $\mathbf{f}(\mathbf{u}) = \boldsymbol{0}$. (Note that the system has four variables to go with the four equations.)

**(b)** Write out the Jacobian matrix of this system.

**(c)** Use {numref}`Function {number} <function-newtonsys>` with different initial guesses to find the two roots of this system. Which is the closest point to $(5,4,3)$, and which is the farthest?
``````

``````{exercise}
:label: problem-newtonsys-circlefit
⌨  Any three noncollinear points in the plane determine a unique circle. Suppose the points are given as $(x_i,y_i)$ for $i=1,2,3$. We can define the circle in terms of its center $(a,b)$ and radius $r$. Then 

$$f_i(a,b,r) = (a-x_i)^2 + (b-y_i)^2 - r^2$$ 

should be made zero for all $i=1,2,3$. This defines a nonlinear system $\mathbf{f}(\mathbf{v})=\boldsymbol{0}$ for $\mathbf{v}=[a,b,r]$. 

Use {numref}`Function {number} <function-newtonsys>` on this system to find the circle passing through $(-5,0)$, $(1,-3)$, and $(4,2)$. Make a plot that shows you found the correct circle.
``````
