---
numbering:
  enumerator: 2.8.%s
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

(section-linsys-condition-number)=

# Conditioning of linear systems

```{index} condition number; of linear system
```

We are ready to consider the conditioning of solving the square linear system $\mathbf{A}\mathbf{x}=\mathbf{b}$. In this problem, the data are $\mathbf{A}$ and $\mathbf{b}$, and the solution is $\mathbf{x}$. Both data and result are multidimensional, so we will use norms to measure their magnitudes.

The motivation for the definition of relative condition number in Chapter 1 was to quantify the response of the result to perturbations of the data. For simplicity, we start by allowing perturbations to $\mathbf{b}$ only while $\mathbf{A}$ remains fixed.

Let $\mathbf{A}\mathbf{x}=\mathbf{b}$ be perturbed to

```{math}
  \mathbf{A}(\mathbf{x}+\mathbf{h}) = \mathbf{b}+\mathbf{d}.
```

The condition number should be the relative change in the solution divided by relative change in the data,

```{math}
  \frac{\quad\dfrac{\| \mathbf{h} \| }{\| \mathbf{x} \| }\quad}{\dfrac{\| \mathbf{d} \| }{\| \mathbf{b} \|}} 
  = \frac{\| \mathbf{h} \|\;\| \mathbf{b} \| }{\| \mathbf{d} \|\; \| \mathbf{x} \| }.
```

We can bound $\| \mathbf{h} \|$ in terms of $\| \mathbf{d} \|$:

```{math}
\begin{split}
  \mathbf{A}\mathbf{x} +  \mathbf{A} \mathbf{h} &= \mathbf{b} + \mathbf{d}, \\
  \mathbf{A} \mathbf{h} &= \mathbf{d},\\
  \mathbf{h} &= \mathbf{A}^{-1} \mathbf{d},\\
  \| \mathbf{h} \| &\le \| \mathbf{A}^{-1}\| \,\| \mathbf{d} \|,
\end{split}
```

where we have applied $\mathbf{A}\mathbf{x}=\mathbf{b}$ and {eq}`normineq1`.
Since also $\mathbf{b}=\mathbf{A}\mathbf{x}$ implies $\| \mathbf{b} \|\le
\| \mathbf{A} \|\, \| \mathbf{x} \|$, we derive

```{math}
   \frac{\| \mathbf{h} \|\; \| \mathbf{b} \|}{\| \mathbf{d} \|\; \| \mathbf{x} \|} 
   \le \frac{\bigl(\| \mathbf{A}^{-1} \|\, \| \mathbf{d} \|\bigr)
    \bigl(\| \mathbf{A} \|\,\| \mathbf{x} \|\bigr)}{\| \mathbf{d} \|\,\| \mathbf{x} \|} 
    = \| \mathbf{A}^{-1}\| \, \| \mathbf{A} \|.
```

```{index} condition number; of a matrix
```

It is possible to show that this bound is tight, in the sense that the inequalities are in fact equalities for some choices of $\mathbf{b}$ and $\mathbf{d}$. This result motivates a new definition.

::::{prf:definition} Matrix condition number
:label: definition-matrixcond
The {term}`matrix condition number` of an invertible square matrix $\mathbf{A}$ is

```{math}
:label: mxcond
\kappa(\mathbf{A}) = \| \mathbf{A}^{-1}\| \, \| \mathbf{A} \|.
```

This value depends on the choice of norm; a subscript on $\kappa$ such as $1$, $2$, or $\infty$ is used if clarification is needed. If $\mathbf{A}$ is singular, we define $\kappa(\mathbf{A}) = \infty$.
::::

## Main result

The matrix condition number {eq}`mxcond` is equal to the condition number of solving a linear system of equations. Although we derived this fact above only for perturbations of $\mathbf{b}$, a similar statement holds when $\mathbf{A}$ is perturbed.

Using a traditional $\Delta$ notation for the perturbation in a quantity, we can write the following.

````{prf:theorem} Conditioning of linear systems
If $\mathbf{A}(\mathbf{x} + \Delta \mathbf{x}) = \mathbf{b} + \Delta \mathbf{b}$, then

```{math}
:label: linsyscondb
\frac{\| \Delta \mathbf{x} \|}{\| \mathbf{x} \|} \le \kappa(\mathbf{A}) \frac{\| \Delta \mathbf{b} \|}{\| \mathbf{b} \|}.
```

If $(\mathbf{A}+\Delta \mathbf{A}) (\mathbf{x} + \Delta \mathbf{x}) = \mathbf{b}$, then

```{math}
:label: linsyscondA
\frac{\| \Delta \mathbf{x} \|}{\| \mathbf{x} \|} \le \kappa(\mathbf{A}) \frac{\| \Delta \mathbf{A} \|}{\| \mathbf{A} \|},
```

in the limit $\| \Delta \mathbf{A} \| \to 0$.
````

Note that for any induced matrix norm,

```{math}
  1 = \| \mathbf{I} \| = \| \mathbf{A} \mathbf{A}^{-1} \| \le \| \mathbf{A} \|\, \| \mathbf{A}^{-1} \| = \kappa(\mathbf{A}).
```

A condition number of 1 is the best we can hope for—in that case, the relative perturbation of the solution has the same size as that of the data.  A condition number of size $10^t$ indicates that in floating-point arithmetic, roughly $t$ digits are lost (i.e., become incorrect) in computing the solution $\mathbf{x}$. And if $\kappa(\mathbf{A}) > \epsilon_\text{mach}^{-1}$, then for computational purposes the matrix is effectively singular.

::::{prf:example} Matrix condition number
:label: demo-condition-bound


```{index} ! Python; cond
```

The function `cond` from `scipy.linalg` is used to computes matrix condition numbers. By default, the 2-norm is used. As an example, the family of *Hilbert matrices* is famously badly conditioned. Here is the $6\times 6$  case. 

```{code-cell} 
A = array([ 
    [1/(i + j + 2) for j in range(6)] 
    for i in range(6) 
    ])
print(A)
```

```{code-cell} 
from numpy.linalg import cond
from scipy import linalg
kappa = cond(A)
print(f"kappa is {kappa:.3e}")
```

Next we engineer a linear system problem to which we know the exact answer.

```{code-cell} 
x_exact = 1.0 + arange(6)
b = A @ x_exact
```

Now we perturb the data randomly with a vector of norm $10^{-12}$. 

```{code-cell} 
dA = random.randn(6, 6)
dA = 1e-12 * (dA / linalg.norm(dA, 2))
db = random.randn(6)
db = 1e-12 * (db / linalg.norm(db, 2))
```

We solve the perturbed problem using built-in pivoted LU and see how the solution was changed.

```{code-cell} 
x = linalg.solve(A + dA, b + db) 
dx = x - x_exact
```

Here is the relative error in the solution.

```{code-cell} 
print(f"relative error is {linalg.norm(dx) / linalg.norm(x_exact):.2e}")
```

And here are upper bounds predicted using the condition number of the original matrix. 

```{code-cell} 
print(f"b_bound: {kappa * 1e-12 / linalg.norm(b):.2e}")
print(f"A_bound: {kappa * 1e-12 / linalg.norm(A, 2):.2e}")
```

Even if we don't make any manual perturbations to the data, machine epsilon does when we solve the linear system numerically.

```{code-cell}
x = linalg.solve(A, b)
print(f"relative error: {linalg.norm(x - x_exact) / linalg.norm(x_exact):.2e}")
print(f"rounding bound: {kappa / 2**52:.2e}")

```

Because $\kappa\approx 10^8$, it's possible to lose 8 digits of accuracy in the process of passing from $A$ and $b$ to $x$. That's independent of the algorithm; it's inevitable once the data are expressed in double precision. 

Larger Hilbert matrices are even more poorly conditioned.

```{code-cell} 
A = array([ [1/(i+j+2) for j in range(14)] for i in range(14) ])
kappa = cond(A)
print(f"kappa is {kappa:.3e}")
```

Before we compute the solution, note that $\kappa$ exceeds `1/eps`. In principle we therefore might end up with an answer that is completely wrong (i.e., a relative error greater than 100%).

```{code-cell} 
print(f"rounding bound: {kappa / 2**52:.2e}")
```

```{code-cell} 
x_exact = 1.0 + arange(14)
b = A @ x_exact  
x = linalg.solve(A, b)
```

We got an answer. But in fact, the error does exceed 100%:

```{code-cell} 
print(f"relative error: {linalg.norm(x - x_exact) / linalg.norm(x_exact):.2e}")
```

::::

## Residual and backward error

Suppose that $\mathbf{A}\mathbf{x}=\mathbf{b}$ and $\tilde{\mathbf{x}}$ is a computed estimate of the solution $\mathbf{x}$. The most natural quantity to study is the error, $\mathbf{x}-\tilde{\mathbf{x}}$. Normally we can't compute it because we don't know the exact solution. However, we can compute something related.

```{index} ! residual; of a linear system
```

::::{prf:definition} Residual of a linear system
:label: definition-linresidual
For the problem $\mathbf{A}\mathbf{x}=\mathbf{b}$, the **residual** at a solution estimate $\tilde{\mathbf{x}}$ is

```{math}
:label: residual
  \mathbf{r} = \mathbf{b} - \mathbf{A}\tilde{\mathbf{x}}.
```

::::

```{index} backward error; in a linear system
```

Obviously, a zero residual means that $\tilde{\mathbf{x}}=\mathbf{x}$, and we have the exact solution. What happens more generally? Note that $\mathbf{A}\tilde{\mathbf{x}}=\mathbf{b}-\mathbf{r}$. That is, $\tilde{\mathbf{x}}$ solves the linear system problem for a right-hand side that is changed by $-\mathbf{r}$. This is precisely what is meant by backward error.

Hence residual and backward error are the same thing for a linear system. What is the connection to the (forward) error? We can reconnect with {eq}`linsyscondb` by the definition $\mathbf{h} = \tilde{\mathbf{x}}-\mathbf{x}$, in which case

$$\mathbf{d} = \mathbf{A}(\mathbf{x}+\mathbf{h})-\mathbf{b}=\mathbf{A}\mathbf{h} = -\mathbf{r}.$$

Thus, {eq}`linsyscondb` is equivalent to

```{math}
:label: residualcond
  \frac{\| \mathbf{x}-\tilde{\mathbf{x}} \|}{\| \mathbf{x} \|} \le
  \kappa(\mathbf{A}) \frac{\| \mathbf{r} \|}{\| \mathbf{b} \|}.
```

Equation {eq}`residualcond` says that the gap between relative error and the relative residual is a multiplication by the matrix condition number.

```{prf:observation}
When solving a linear system, all that can be expected is that the backward error, not the error, is small.
```

## Exercises

``````{exercise}
:label: problem-cond-hilbert

⌨ Refer to @demo-condition-bound for the definition of a Hilbert matrix. Make a table of the values of $\kappa(\mathbf{H}_n)$ in the 2-norm for $n=2,3,\ldots,16$. Speculate as to why the growth of $\kappa$ appears to slow down at $n=13$.
``````

``````{exercise}
:label: problem-cond-verify

⌨ The purpose of this problem is to verify, like in @demo-condition-bound, the error bound

```{math}
:numbered: false
\frac{\| \mathbf{x}-\tilde{\mathbf{x} \|}}{\| \mathbf{x} \|} \le \kappa(\mathbf{A})
\frac{\| \mathbf{h} \|}{\| \mathbf{b} \|}.
```

Here $\tilde{\mathbf{x}}$ is a numerical approximation to the exact solution $\mathbf{x}$, and $\mathbf{h}$ is an unknown perturbation caused by machine roundoff. We will assume that $\| \mathbf{h} \|/\| \mathbf{b} \|$ is roughly `eps()`.

You will need to import the `rogues` package from PyPi. For each $n=10,20,\ldots,70$, let `A = rogues.prolate(n, 0.4)` and let $\mathbf{x}$ have components $x_k=k/n$ for $k=1,\ldots,n$. Define `b = A@x` and let `xtilde` be the solution produced numerically by `numpy.linalg.solve`.

Make a table including columns for $n$, the condition number of $\mathbf{A}$, the observed relative error in $\tilde{\mathbf{x}}$, and the right-hand side of the inequality above. You should find that the inequality holds in every case.
``````

``````{exercise}
:label: problem-cond-triangular

⌨ @problem-linearsystems-triangillcond suggests that the solutions of linear systems

```{math}
:numbered: false
\mathbf{A} = \begin{bmatrix} 1 & -1 & 0 & \alpha-\beta & \beta \\ 0 & 1 & -1 &
0 & 0 \\ 0 & 0 & 1 & -1 & 0 \\ 0 & 0 & 0 & 1 & -1  \\ 0 & 0 & 0 & 0 & 1
\end{bmatrix}, \quad
\mathbf{b} = \begin{bmatrix} \alpha \\ 0 \\ 0 \\ 0 \\ 1 \end{bmatrix}
```

become less accurate as $\beta$ increases. Using $\alpha=0.1$ and $\beta=10,10^2,\ldots,10^{12}$, make a table with columns for $\beta$, $|x_1-1|$, and the condition number of the matrix.
``````

``````{exercise}
:label: problem-cond-vandermonde

⌨ Let $\mathbf{A}_n$ denote the $(n+1)\times(n+1)$ version of the Vandermonde matrix in @vandersystem based on the equally spaced interpolation nodes $t_i=i/n$ for $i=0,\ldots,n$. Using the 1-norm, plot $\kappa(\mathbf{A}_n)$ as a function of $n$ for $n=4,5,6,\ldots,20$, using a log scale on the $y$-axis. (The graph is nearly a straight line.)

``````

``````{exercise}
:label: problem-cond-unpivoted

⌨ The matrix $\mathbf{A}$ in {eq}`plu-stab-A` has unpivoted LU factors given in {eq}`plu-stab-LU` as a function of parameter $\epsilon$. For $\epsilon = 10^{-2},10^{-4},\ldots,10^{-10}$, make a table with columns for $\epsilon$, $\kappa(\mathbf{A})$, $\kappa(\mathbf{L})$, and $\kappa(\mathbf{U})$. (This shows that solution via unpivoted LU factorization is arbitrarily unstable.)
``````

``````{exercise}
:label: problem-cond-inverse

✍  Define $\mathbf{A}_n$ as the $n\times n$ matrix $\displaystyle\begin{bmatrix}
1 & -2 & & &\\
& 1 & -2 & & \\
& & \ddots & \ddots & \\
& & & 1 & -2 \\
& & & & 1
\end{bmatrix}.$

**(a)** Write out $\mathbf{A}_2^{-1}$ and $\mathbf{A}_3^{-1}$.

**(b)** Write out $\mathbf{A}_n^{-1}$ in the general case $n>1$. (If necessary, look at a few more cases on the computer until you are certain of the pattern.) Make a clear argument why it is correct.

**(c)** Using the $\infty$-norm, find $\kappa(\mathbf{A}_n)$.
``````

``````{exercise}
:label: problem-cond-product

✍ **(a)** Prove that for $n\times n$ nonsingular matrices $\mathbf{A}$ and $\mathbf{B}$, $\kappa(\mathbf{A}\mathbf{B})\le \kappa(\mathbf{A})\kappa(\mathbf{B})$.

**(b)** Show by means of an example that the result of part (a) cannot be an equality in general.
``````

``````{exercise}
:label: problem-cond-diagonal
✍  Let $\mathbf{D}$ be a diagonal $n\times n$ matrix, not necessarily invertible. Prove that in the 1-norm,

```{math}
:numbered: false
\kappa(\mathbf{D}) = \frac{\max_i |D_{ii}|}{\min_i |D_{ii}|}.
```

(Hint: See @problem-norms-diagnorm.)
``````
