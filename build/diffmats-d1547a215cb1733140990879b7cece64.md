---
numbering:
  enumerator: 10.3.%s
---
(section-bvp-diffmats)=
# Differentiation matrices

```{index} finite differences; matrix
```
In {numref}`section-localapprox-finitediffs` we used finite differences to turn a discrete collection of function values into an estimate of the derivative of the function at a point. Just as with differentiation in elementary calculus, we can generalize differences at a point into an operation that maps discretized functions to discretized derivatives.


## Matrices for finite differences

We first discretize the interval $x\in[a,b]$ into equal pieces of length $h=(b-a)/n$, leading to the nodes

$$
x_i=a+i h, \qquad  i=0,\ldots,n.
$$

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

:::{math}
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
:::

```{index} ! differentiation matrix
```
Here as elsewhere, elements of $\mathbf{D}_x$ that are not shown are zero. We call $\mathbf{D}_x$ a **differentiation matrix**. Each row of $\mathbf{D}_x$ gives the weights of the finite-difference formula being used at one of the nodes.

```{index} order of accuracy; of a finite-difference formula
```

The differentiation matrix $\mathbf{D}_x$ in {eq}`diffmat11` is not a unique choice. We are free to use whatever finite-difference formulas we like in each row. However, it makes sense to choose rows that are as similar as possible. Using second-order centered differences where possible and second-order one-sided formulas (see {numref}`table-FDforward`) at the boundary points leads to

:::{math}
:label: diffmat12b
\mathbf{D}_x = \frac{1}{h}
\begin{bmatrix}
-\frac{3}{2} & 2    & -\frac{1}{2}   &        &        &  \\[1mm]
-\frac{1}{2} & 0    & \frac{1}{2}    &        &        &  \\[1mm]
& -\frac{1}{2} & 0      & \frac{1}{2}    &        &  \\
&      & \ddots & \ddots & \ddots &  \\
&      &        & -\frac{1}{2}   & 0      & \frac{1}{2} \\[1mm]
&      &        & \frac{1}{2}    & -2     & \frac{3}{2}
\end{bmatrix}.
:::

```{index} banded matrix, sparse matrix
```

The differentiation matrices so far are banded matrices, i.e., all the nonzero values are along diagonals close to the main diagonal.[^sparseconv]

[^sparseconv]: In order to exploit this structure efficiently in Julia, these matrices first need to be constructed as or converted to sparse or tridiagonal form.

## Second derivative

Similarly, we can define differentiation matrices for second derivatives. For example,

:::{math}
:label: diffmat22
\begin{bmatrix}
f''(x_0) \\[1mm] f''(x_1) \\[1mm] f''(x_2) \\[1mm] \vdots \\[1mm] f''(x_{n-1}) \\[1mm] f''(x_n)
\end{bmatrix}
\approx
\frac{1}{h^2}
\begin{bmatrix}
2 & -5 & 4      & -1     &        &     \\[1mm]
1 & -2 & 1      &        &        &    \\[1mm]
& 1  & -2     & 1      &        &    \\[1mm]
&    & \ddots & \ddots & \ddots &    \\[1mm]
&        &        & 1      & -2 & 1 \\[1mm]
&        &     -1  & 4      & -5 & 2
\end{bmatrix}
\begin{bmatrix}
f(x_0) \\[1mm] f(x_1) \\[1mm] f(x_2) \\[1mm] \vdots \\[1mm] f(x_{n-1}) \\[1mm] f(x_n)
\end{bmatrix} = \mathbf{D}_{xx} \mathbf{f}.
:::

We have multiple choices again for $\mathbf{D}_{xx}$, and it need not be the square of any particular $\mathbf{D}_x$. As pointed out in {numref}`section-localapprox-finitediffs`, squaring the first derivative is a valid approach but would place entries in $\mathbf{D}_{xx}$ farther from the diagonal than is necessary.

(function-diffmat2)=
``````{prf:algorithm} diffmat2
`````{tab-set} 
````{tab-item} Julia
:sync: julia
:::{embed} #function-diffmat2-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #function-diffmat2-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #function-diffmat2-python
:::
````
`````
``````

Together the matrices {eq}`diffmat12b` and {eq}`diffmat22` give second-order approximations of the first and second derivatives at all nodes. These matrices, as well as the nodes $x_0,\ldots,x_n$, are returned by {numref}`Function {number} <function-diffmat2>`.

(demo-diffmats-2nd)=
::::{prf:example} Differentiation matrices
`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-diffmats-2nd-julia
:::
````

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-diffmats-2nd-matlab
:::
````

````{tab-item} Python
:sync: python
:::{embed} #demo-diffmats-2nd-python
:::
````
`````
::::


## Spectral differentiation

Recall that finite-difference formulas are derived in three steps:

1. Choose a node index set $S$ near node $i$.
2. Interpolate with a polynomial using the nodes in $S$.
3. Differentiate the interpolant and evaluate at node $i$.

```{index} Chebyshev points; second kind
```
We can modify this process by using a global interpolant, either polynomial or trigonometric, as in [Chapter 9](../globalapprox/overview). Rather than choosing a different index set for each node, we use all of the nodes each time. 

In a nonperiodic setting we use Chebyshev second-kind points for stability:

$$
x_k = -\cos\left(\frac{k \pi}{n}\right), \qquad k=0,\ldots,n;
$$

```{index} differentiation matrix
```
then the resulting **Chebyshev differentiation matrix** has entries 

:::{math}
:label: chebdiffmat
  \begin{gathered}
    D_{00} = \dfrac{2n^2+1}{6}, \qquad D_{n n} = -\dfrac{2n^2+1}{6}, \\
    D_{ij} =
    \begin{cases}
      -\dfrac{x_i}{2(1-x_i^2)}, & i=j, \\[4mm]
      \dfrac{c_i}{c_j}\, \dfrac{(-1)^{i+j}}{x_i-x_j}, & i\neq j,
    \end{cases}
  \end{gathered}
:::

where $c_0=c_n=2$ and $c_i=1$ for $i=1,\ldots,n-1$. Note that this matrix is dense. The simplest way to compute a second derivative is by squaring $\mathbf{D}_x$, as there is no longer any concern about the bandwidth of the result.

{numref}`Function {number} <function-diffcheb>` returns these two matrices. The function uses a change of variable to transplant the standard $[-1,1]$ for Chebyshev nodes to any $[a,b]$. It also takes a different approach to computing the diagonal elements of $\mathbf{D}_x$ than the formulas in {eq}`chebdiffmat` (see {ref}`Exercise 5 <problem-diffmats-negsumtrick>`).

(function-diffcheb)=
``````{prf:algorithm} diffcheb
`````{tab-set} 
````{tab-item} Julia
:sync: julia
:::{embed} #function-diffcheb-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #function-diffcheb-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #function-diffcheb-python
:::
````
`````
``````

(demo-diffmats-cheb)=
::::{prf:example} Chebyshev differentiation matrices
`````{tab-set}
````{tab-item} Julia
`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-diffmats-cheb-julia
:::
````

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-diffmats-cheb-matlab
:::
````

````{tab-item} Python
:sync: python
:::{embed} #demo-diffmats-cheb-python
:::
````
`````
::::

According to {numref}`Theorem %s <theorem-spectral>`, the convergence of polynomial interpolation to $f$ using Chebyshev nodes is spectral if $f$ is analytic (at least having infinitely many derivatives) on the interval. The derivatives of $f$ are also approximated with spectral accuracy.

## Exercises

1. 
    **(a)** ✍ Calculate $\mathbf{D}_x^2$ using {eq}`diffmat11` to define $\mathbf{D}_x$.
    
    **(b)** ⌨ Repeat the experiment of {numref}`Demo %s <demo-diffmats-2nd>`, but using this version of $\mathbf{D}_x^2$ to estimate $f''$. What is the apparent order of accuracy in $f''$?

2. 
    **(a)** ✍ Find the derivative of $f(x) =\operatorname{sign}(x)x^2$ on the interval $[-1,1]$. (If this gives you trouble, use an equivalent piecewise definition of $f$.) What is special about this function at $x=0$? 

    **(b)** ⌨ Adapt {numref}`Demo %s <demo-diffmats-2nd>` to operate on the function from part (a), computing only the first derivative. What is the observed order of accuracy?
  
    **(c)** ✍ Show that for even values of $n$, there is only one node at which the error for computing $f'$ in part (b) is nonzero.

3. ⌨ To get a fourth-order accurate version of $\mathbf{D}_x$, five points per row are needed, including two special rows at each boundary. For a fourth-order $\mathbf{D}_{xx}$, five symmetric points per row are needed for interior rows and six points are needed for the rows near a boundary.

    **(a)** Modify {numref}`Function {number} <function-diffmat2>` to a function `diffmat4`, which outputs fourth-order accurate differentiation matrices. You may want to use {numref}`Function {number} <function-fdweights>`.
    
    **(b)** Repeat the experiment of {numref}`Demo %s <demo-diffmats-2nd>` using `diffmat4` in place of {numref}`Function {number} <function-diffmat2>`, and compare observed errors to fourth-order accuracy.

4. ✍ Explain in detail how lines 23-24 in {numref}`Function {number} <function-diffcheb>` correctly change the interval from $[-1,1]$ to $[a,b]$.

(problem-diffmats-negsumtrick)=
5. 
    **(a)** ✍ What is the derivative of a constant function?
    
    **(b)** ✍ Explain why for any reasonable differentiation matrix $\mathbf{D}$, we should find $\displaystyle \sum_{j=0}^nD_{ij}=0$ for all $i$.
    
    **(c)** ✍ What does this have to do with {numref}`Function {number} <function-diffcheb>`? Refer to specific line(s) in the function for your answer.

6. Define the $(n+1)\times (n+1)$ matrix $\mathbf{T} = \displaystyle \begin{bmatrix}
    1 & & & \\ 1 & 1 & & \\ \vdots & & \ddots & \\ 1 & 1 & \cdots & 1
    \end{bmatrix}$.

    **(a)** ✍ Write out $\mathbf{T}\mathbf{u}$ for a generic vector $\mathbf{u}$ (start with a zero index). How is this like integration?
    
    **(b)** ✍ Find the inverse of $\mathbf{T}$ for any $n$.  (You can use Julia to find the pattern, but show that the result is correct in general.) What does this have to do with the inverse of integration?
