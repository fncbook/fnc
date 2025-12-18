---
numbering:
  enumerator: 10.3.%s
kernelspec:
  display_name: MATLAB
  language: matlab
  name: jupyter_matlab_kernel
---
```{code-cell}
:tags: [remove-cell]
clear all
format short
set(0, 'defaultaxesfontsize', 12)
set(0, 'defaultlinelinewidth', 1.5)
set(0, 'defaultFunctionLinelinewidth', 1.5)
set(0, 'defaultscattermarkerfacecolor', 'flat')
gcf;
set(gcf, 'Position', [0 0 600 350])
addpath FNC-matlab
```

(section-bvp-diffmats)=

# Differentiation matrices

```{index} finite differences; matrix
```

In {numref}`section-localapprox-finitediffs` we used finite differences to turn a discrete collection of function values into an estimate of the derivative of the function at a point. Just as with differentiation in elementary calculus, we can generalize differences at a point into an operation that maps discretized functions to discretized derivatives.

```{index} ! differentiation matrix
```

:::{prf:definition} Differentiation matrix
:label: definition-diffmat
A {term}`differentiation matrix` is a square matrix whose associated linear transformation maps a vector of function values at nodes to a vector of approximate derivative values at the same nodes.
:::

The derivative of a function is a unique result. The same is *not* true for our finite-dimensional approximation of the derivative, though.

:::{important}
There are many ways to approximate the derivative operator by a differentiation matrix.
:::

## Matrices for finite differences

:::{note}
In [Chapter 5](#chapter-localapprox) and [Chapter 9](#chapter-globalapprox) we used the notation $t_k$ for nodes to more clearly distinguish them from the continuous variable $x$. Since we will soon have problems that discretize both $x$ and $t$ variables, though, we now switch to using the same letter for nodes as the continuous variable.
:::

We first discretize the interval $x\in[a,b]$ into equal pieces of length $h=(b-a)/n$, leading to the nodes

```{math}
:label: bvp-fdnodes
x_i=a + i h, \quad  i=0,\ldots,n, \qquad h=\frac{b-a}{n}.
```

Note again that our node indexing scheme starts at zero. Our goal is to find a vector $\mathbf{g}$ such that $g_i \approx f'(x_i)$ for $i=0,\ldots,n$. Our first try is the forward difference formula {eq}`forwardFD11`,

$$
g_i = \frac{f_{i+1}-f_i}{h}, \qquad t=0,\ldots,n-1.
$$

However, this leaves $g_n$ undefined, because the formula would refer to the unavailable value $f_{n+1}$. For $g_n$ we could resort to the backward difference

$$
g_n = \frac{f_{n}-f_{n-1}}{h}.
$$

We can summarize the entire set of formulas by defining

$$
\mathbf{f} =
\begin{bmatrix}
f(x_0) \\[1mm] f(x_1) \\[1mm] \vdots \\[1mm] f(x_{n-1}) \\[1mm] f(x_n)
\end{bmatrix}, \quad
$$

and then the vector equation

```{math}
:label: diffmat11
\begin{bmatrix}
f'(x_0) \\[1mm] f'(x_1) \\[1mm] \vdots \\[1mm] f'(x_{n-1}) \\[1mm] f'(x_n)
\end{bmatrix}
\approx
\mathbf{D}_x \mathbf{f}, \qquad
\mathbf{D}_x
= \frac{1}{h}
\begin{bmatrix}
-1 & 1 & & & \\[1mm]
& -1 & 1 & & \\[1mm]
& & \ddots & \ddots & \\[1mm]
& & & -1 & 1 \\[1mm]
& & & -1 & 1
\end{bmatrix}.
```

The matrix $\mathbf{D}_x$ in @diffmat11 is a finite-difference differentiation matrix.  Here as elsewhere, elements of $\mathbf{D}_x$ that are not shown are zero. Each row of $\mathbf{D}_x$ gives the weights of a finite-difference formula being used at one of the nodes.

```{index} order of accuracy; of a finite-difference formula
```

The differentiation matrix $\mathbf{D}_x$ in {eq}`diffmat11` is not a unique choice. In fact, it's about the least accurate choice possible, as explained in {numref}`section-localapprox-fd-converge`. We are theoretically free to use whatever finite-difference formulas we want in each row, such as those in @table-fdcenter and @table-fdforward, although it makes sense to choose rows that are as similar as possible. For instance, using second-order centered differences where possible and second-order one-sided formulas at the boundary points leads to

```{index} finite-difference formula; first derivative
```

```{index} differentiation matrix
```

```{prf:definition} Finite-difference differentiation matrix for a first derivative
:label: definition-diffmatfd1

A second-order accurate differentiation matrix for the first derivative using finite differences is 

```{math}
:label: diffmat12b
\mathbf{D}_x = \frac{1}{h}
\begin{bmatrix}
-\frac{3}{2} & 2    & -\frac{1}{2}   &        &        &  \\[1mm]
-\frac{1}{2} & 0    & \frac{1}{2}    &        &        &  \\[1mm]
& -\frac{1}{2} & 0      & \frac{1}{2}    &        &  \\
&      & \ddots & \ddots & \ddots &  \\
&      &        & -\frac{1}{2}   & 0      & \frac{1}{2} \\[1mm]
&      &        & \frac{1}{2}    & -2     & \frac{3}{2}
\end{bmatrix},
```

to be applied to a vector of function values at the nodes @bvp-fdnodes.

```{index} banded matrix, sparse matrix
```

:::{note}
The differentiation matrix @diffmat12b is sparse and banded, i.e., all the nonzero values are along diagonals close to the main diagonal.
:::

::::{prf:example} Differentiation matrix for finite differences
:label: example-diffmatfd1

Let's use @diffmat12b to compute the derivative of $f(x) = \sin(\pi x)$ over the interval $[-1, 0]$ with $n=4$. This gives $h = 1 / 4$ and the nodes $x_i = -1 + i/4$ for $i=0,\ldots,4$. The matrix $\mathbf{D}_x$ is

```{math}
\mathbf{D}_x = 4
\begin{bmatrix}
-3/2 & 2    & -1/2   &        &    \\  
-1/2 & 0    & 1/2    &        &   \\ 
 & -1/2 & 0      & 1/2    &      \\
 & & -1/2 & 0      & 1/2      \\
& & 1/2 & -2    & 3/2  
\end{bmatrix}.
```

The discrete approximation to $f'(x) = \pi \cos(\pi x)$ is

```{math}
& 
\begin{bmatrix}
-6 & 8    & -2   &        &   \\  
-2 & 0    & 2    &        &    \\
 & -2 & 0      & 2    &      \\
 & & -2 & 0      & 2      \\
& & 2 & -8    & 6  
\end{bmatrix}
\begin{bmatrix}
\sin(-\pi) \\[1mm] \sin(-3\pi/4) \\[1mm] \sin(-\pi/2) \\[1mm] \sin(-\pi/4) \\[1mm] \sin(0)
\end{bmatrix} \\ 
& =
\begin{bmatrix}
-6 & 8    & -2   &        &   \\   
-2 & 0    & 2    &        &  \\  
 & -2 & 0      & 2    &      \\
 & & -2 & 0      & 2      \\
& & 2 & -8  & 6  
\end{bmatrix}
\begin{bmatrix}
0 \\ -1/\sqrt{2} \\ -1 \\ -1/\sqrt{2} \\ 0
\end{bmatrix} \\ 
& =
\begin{bmatrix}
2 - 8\sqrt{2} + 2\\ -2 \\ 0 \\ 2 \\ -2 + 8/\sqrt{2}
\end{bmatrix} \\ 
& \approx 
\begin{bmatrix}
-3.657 \\ -2 \\ 0 \\ 2 \\ 3.657
\end{bmatrix}.
```

This result is not especially close to the exact values of $f'$ at the nodes, which are $\bigl[ -3.142, -2.221, 0, 2.221, 3.142 \bigr]$, but we should not expect it to be for such a small $n$.
::::

## Second derivative

In a TPBVP, we will need to take the second derivative of the unknown solution. One option is to apply a first derivative twice, that is, to multiply the value vector $\mathbf{f}$ on the left by $\mathbf{D}_x^2$. This is usually not the best option, however (see {numref}`section-localapprox-finitediffs`). Instead, the following is usually a better choice:

```{index} finite-difference formula; second derivative
```

```{index} differentiation matrix
```

````{prf:definition} Finite-difference differentiation matrix for a second derivative
:label: definition-diffmatfd2

To second order,

$$
\begin{bmatrix}
f''(x_0) \\[1mm] f''(x_1) \\[1mm] f''(x_2) \\[1mm] \vdots \\[1mm] f''(x_{n-1}) \\[1mm] f''(x_n)
\end{bmatrix}
\approx \mathbf{D}_{xx} \begin{bmatrix}
f(x_0) \\[1mm] f(x_1) \\[1mm] f(x_2) \\[1mm] \vdots \\[1mm] f(x_{n-1}) \\[1mm] f(x_n)
\end{bmatrix}, 
$$

where

```{math}
:label: diffmat22
\mathbf{D}_{xx} = 
\frac{1}{h^2}
\begin{bmatrix}
2 & -5 & 4      & -1     &        &     \\[1mm]
1 & -2 & 1      &        &        &    \\[1mm]
& 1  & -2     & 1      &        &    \\[1mm]
&    & \ddots & \ddots & \ddots &    \\[1mm]
&        &        & 1      & -2 & 1 \\[1mm]
&        &     -1  & 4      & -5 & 2
\end{bmatrix},
```

and the nodes are given in @bvp-fdnodes.
````

## Implementation

Together, the matrices {eq}`diffmat12b` and {eq}`diffmat22` give second-order approximations of the first and second derivatives at all nodes. These matrices, as well as the nodes $x_0,\ldots,x_n$, are returned by {numref}`Function {number} <function-diffmat2>`.

``````{prf:algorithm} diffmat2
:label: function-diffmat2
```{literalinclude} ../FNC_matlab/diffmat2.m
:language: matlab
:linenos: true
```
``````

::::{prf:example} Differentiation matrices
:label: demo-diffmats-2nd

We test first-order and second-order differentiation matrices for the function $x + \exp(\sin 4x)$ over $[-1,1]$.

```{code-cell}
f = @(x) x + exp(sin(4*x));
```

For reference, here are the exact first and second derivatives.

```{code-cell}
df_dx = @(x) 1 + 4 * exp(sin(4*x)) .* cos(4*x);
d2f_dx2 = @(x) 4 * exp(sin(4*x)) .* (4*cos(4*x).^2 - 4*sin(4*x));
```

We discretize on equally spaced nodes and evaluate $f$ at the nodes.

```{code-cell}
[t, Dx, Dxx] = diffmat2(12, [-1 1]);
y = f(t);
```

Then the first two derivatives of $f$ each require one matrix-vector multiplication.

```{code-cell}
yx = Dx * y;
yxx = Dxx * y;
```

The results show poor accuracy for this small value of $n$.

```{code-cell}
clf,  subplot(2, 1, 1)
fplot(df_dx, [-1, 1]),  hold on
plot(t, yx, 'ko')
xlabel('x'),  ylabel('f''(x)')
subplot(2, 1, 2)
fplot(d2f_dx2, [-1, 1]),  hold on
plot(t, yxx, 'ko')
xlabel('x'),  ylabel('f''''(x)')
```

An convergence experiment confirms the order of accuracy. Because we expect an algebraic convergence rate, we use a log-log plot of the errors.

```{code-cell}
n = round( 2 .^ (4:.5:11)' );
err = zeros(length(n), 2);
for k = 1:length(n)
    [t, Dx, Dxx] = diffmat2(n(k), [-1, 1]);
    y = f(t);
    err(k, 1) = norm(df_dx(t) - Dx * y, Inf);
    err(k, 2) = norm(d2f_dx2(t) - Dxx * y, Inf);
end
clf
loglog(n, err, 'o-'), hold on
loglog(n, 100 * n.^(-2), 'k--')
legend("f'", "f''", '2nd order')
xlabel('n'),  ylabel('max error')
title('Convergence of finite differences')
```

::::

## Spectral differentiation

Recall that finite-difference formulas are derived in three steps:

1. Choose a node index set $S$ near node $i$.
2. Interpolate with a polynomial using the nodes in $S$.
3. Differentiate the interpolant and evaluate at node $i$.

```{index} Chebyshev points; second kind
```

We can modify this process by using a global interpolant, either polynomial or trigonometric, as in [Chapter 9](../globalapprox/overview). Rather than choosing a different index set for each node, we use all the nodes each time.

In a nonperiodic setting, we use Chebyshev second-kind points for stability:

```{math}
:label: bvp-chebnodes
x_k = -\cos\left(\frac{k \pi}{n}\right), \qquad k=0,\ldots,n.
```

```{index} differentiation matrix
```

````{prf:definition} Chebyshev differentiation matrix
:label: definition-diffmatcheb

Using a vector $\mathbf{f}$ of function values at the Chebyshev nodes @bvp-chebnodes, the vector $\mathbf{D}_x \mathbf{f}$ will have approximate first derivative values at those nodes, where the entries of $\mathbf{D}_x$ are given by 

```{math}
:label: chebdiffmat
  \begin{gathered}
    D_{00} = \dfrac{2n^2+1}{6}, \qquad D_{n n} = -\dfrac{2n^2+1}{6}, \\
    D_{ij} =
    \begin{cases}
      -\dfrac{x_i}{2(1-x_i^2)}, & i=j, \\[4mm]
      \dfrac{c_i}{c_j}\, \dfrac{(-1)^{i+j}}{x_i-x_j}, & i\neq j,
    \end{cases}
  \end{gathered}
```

where $c_0=c_n=2$ and $c_i=1$ for $i=1,\ldots,n-1$.

For the second derivative, $\mathbf{D}_{xx} = \mathbf{D}_x^2$.
````

:::{note}
The Chebyshev differentiation matrix is not sparse. There are compact formulas available elsewhere for the entries of $\mathbf{D}_{xx}$ ({cite}`trefethenSpectralMethods2000`).
:::

{numref}`Function {number} <function-diffcheb>` returns the matrices from @definition-diffmatcheb. This function uses a change of variable to transplant the standard $[-1,1]$ for Chebyshev nodes to any $[a,b]$. It also takes a different approach to computing the diagonal elements of $\mathbf{D}_x$ than the formulas in {eq}`chebdiffmat` (see @problem-diffmats-negsumtrick).

``````{prf:algorithm} diffcheb
:label: function-diffcheb

```{literalinclude} ../FNC_matlab/diffcheb.m
:language: matlab
:linenos: true
```
``````

::::{prf:example} Chebyshev differentiation matrices
:label: demo-diffmats-cheb

Here is a $4\times 4$ Chebyshev differentiation matrix.

```{code-cell}
[t, Dx] = diffcheb(3, [-1, 1]);
format rat
Dx
```

We again test the convergence rate.

```{code-cell}
f = @(x) x + exp(sin(4*x));
df_dx = @(x) 1 + 4 * exp(sin(4*x)) .* cos(4*x);
d2f_dx2 = @(x) 4 * exp(sin(4*x)) .* (4*cos(4*x).^2 - 4*sin(4*x));

```

```{code-cell}
n = 5:5:70;
err = zeros(length(n), 2);
for k = 1:length(n)
    [t, Dx, Dxx] = diffcheb(n(k), [-1, 1]);
    y = f(t);
    err(k, 1) = norm(df_dx(t) - Dx * y, Inf);
    err(k, 2) = norm(d2f_dx2(t) - Dxx * y, Inf);
end
```

Since we expect a spectral convergence rate, we use a semi-log plot for the error.

```{code-cell}
clf,  format
semilogy(n, err, 'o-'), hold on
legend("f'", "f''")
xlabel('n'),  ylabel('max error')
title('Convergence of finite differences')
```

::::

According to @theorem-spectral, the convergence of polynomial interpolation to $f$ using Chebyshev nodes is spectral if $f$ is analytic (at least having infinitely many derivatives) on the interval. The derivatives of $f$ are also approximated with spectral accuracy.

## Exercises

``````{exercise}
:label: problem-diffmats-square
**(a)** ✍  Using {eq}`diffmat11` to define $\mathbf{D}_x$, calculate $\mathbf{D}_x^2$ in the general case.

**(b)** ⌨ Repeat the convergence experiment in the second part of @demo-diffmats-2nd, but using this version of $\mathbf{D}_x^2$ in place of $\mathbf{D}_{xx}$ to estimate $f''$. Why does it fail to converge as $n\to \infty$?
``````

``````{exercise}
:label: problem-diffmats-jump
**(a)** ✍ Find the derivative of $f(x) =\operatorname{sign}(x)x^2$ on the interval $[-1,1]$. (If this gives you trouble, use an equivalent piecewise definition of $f$.) What is special about this function at $x=0$? 

**(b)** ⌨ Adapt @demo-diffmats-2nd to operate on the function from part (a), computing only the first derivative. What is the observed order of accuracy?

**(c)** ✍ Show that for even values of $n$, there is only one node at which the error for computing $f'$ in part (b) is nonzero.
``````

``````{exercise}
:label: problem-diffmats-4thorder
⌨ To get a fourth-order accurate version of $\mathbf{D}_x$, five points per row are needed, including two special rows at each boundary. For a fourth-order $\mathbf{D}_{xx}$, five symmetric points per row are needed for interior rows and six points are needed for the rows near a boundary.

**(a)** Modify {numref}`Function {number} <function-diffmat2>` to a function `diffmat4`, which outputs fourth-order accurate differentiation matrices. You may want to use {numref}`Function {number} <function-fdweights>`.

**(b)** Repeat the experiment of @function-diffmat2, and compare observed errors to fourth-order accuracy.
``````

``````{exercise}
:label: problem-diffmats-interval
✍ Explain in detail how lines 23-24 in {numref}`Function {number} <function-diffcheb>` correctly change the interval from $[-1,1]$ to $[a,b]$.
``````

``````{exercise}
:label: problem-diffmats-negsumtrick

**(a)** ✍ What is the derivative of a constant function?

**(b)** ✍ Explain why for any reasonable differentiation matrix $\mathbf{D}$, we should find $\displaystyle \sum_{j=0}^nD_{ij}=0$ for all $i$.

**(c)** ✍ What does this have to do with {numref}`Function {number} <function-diffcheb>`? Refer to specific line(s) in the function for your answer.
``````

``````{exercise}
:label: problem-diffmats-integration

Define the $(n+1)\times (n+1)$ matrix $\mathbf{T} = \displaystyle \begin{bmatrix}
1 & & & \\ 1 & 1 & & \\ \vdots & & \ddots & \\ 1 & 1 & \cdots & 1
\end{bmatrix}$.

**(a)** ✍ Write out $\mathbf{T}\mathbf{u}$ for a generic vector $\mathbf{u}$ (start with a zero index). How is this like integration?

**(b)** ✍ Find the inverse of $\mathbf{T}$ for any $n$.  (You can use a computer to find the pattern, but show that the result is correct in general.) What does this have to do with the inverse of integration?
``````
