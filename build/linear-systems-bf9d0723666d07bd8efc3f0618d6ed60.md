---
numbering:
  enumerator: 2.3.%s
---
(section-linsys-linear-systems)=

# Linear systems

We now attend to the central problem of this chapter: Given a square, $n\times n$ matrix $\mathbf{A}$ and an $n$-vector $\mathbf{b}$, find an $n$-vector $\mathbf{x}$ such that $\mathbf{A}\mathbf{x}=\mathbf{b}$. Writing out these equations, we obtain

```{math}
\begin{split}
  A_{11}x_1 + A_{12}x_2 + \cdots + A_{1n}x_n &= b_1, \\
  A_{21}x_1 + A_{22}x_2 + \cdots + A_{2n}x_n &= b_2, \\
  \vdots  \\
  A_{n1}x_1 + A_{n2}x_2 + \cdots + A_{nn}x_n &= b_n.
\end{split}
```

```{index} matrix inverse
```

If $\mathbf{A}$ is invertible, then the mathematical expression of the solution is $\mathbf{x}=\mathbf{A}^{-1}\mathbf{b}$ because

```{math}
\begin{split}
  \mathbf{A}^{-1}\mathbf{b} = \mathbf{A}^{-1} (\mathbf{A} \mathbf{x}) = (\mathbf{A}^{-1}\mathbf{A}) \mathbf{x} = \mathbf{I} \mathbf{x}
  = \mathbf{x}.
\end{split}
```

When $\mathbf{A}$ is singular, then $\mathbf{A}\mathbf{x}=\mathbf{b}$ may have no solution or
infinitely many solutions.

````{prf:example}
:label: example-singmatrix
If we define

```{math}
\mathbf{S} =  \begin{bmatrix}
 0 & 1\\0 & 0
\end{bmatrix},
```

then it is easy to check that for any real value of $\alpha$ we have

```{math}
\mathbf{S}
\begin{bmatrix}
 \alpha \\ 1
\end{bmatrix}
=
\begin{bmatrix}
 1 \\ 0
\end{bmatrix}.
```

Hence the linear system $\mathbf{S}\mathbf{x}=\mathbf{b}$ with $\mathbf{b}=\begin{bmatrix} 1\\0\end{bmatrix}$ has infinitely many solutions. For most other choices of $\mathbf{b}$, the system has no solution.
````

## Don't use the inverse

Matrix inverses are indispensable for mathematical discussion and derivations. However, as you may remember from a linear algebra course, they are not trivial to compute from the entries of the original matrix. You might be surprised to learn that matrix inverses play almost no role in scientific computing.

:::{important}
Computing the inverse of a matrix is not a good way to solve a linear system of equations.
:::

In fact, when we encounter an expression such as $\mathbf{x} = \mathbf{A}^{-1} \mathbf{b}$ in computing, we interpret it as "solve the linear system $\mathbf{A} \mathbf{x} = \mathbf{b}$" and apply whatever algorithm is most expedient based on what we know about $\mathbf{A}$.

```{index} Julia; \\
```

`````{tab-set}
````{tab-item} Julia
:sync: julia
As demonstrated in @demo-interp-vander, the backslash (the `\` symbol, not to be confused with the slash `/` used in web addresses) invokes a linear system solution. 
```` 

````{tab-item} MATLAB
:sync: matlab
In MATLAB, the backslash operator `\` is used to solve linear systems. 
```` 

````{tab-item} Python
:sync: python
In Python, the `numpy.linalg.solve` function is used to solve linear systems. 
```` 
`````

::::{prf:example} Solving linear systems
:label: demo-systems-backslash

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-systems-backslash-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-systems-backslash-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-systems-backslash-python
:::
```` 
`````

::::

(section-linsys-triangular)=

## Triangular systems

```{index} triangular matrix
```

The solution process is especially easy to demonstrate for a system with a **triangular matrix**. For example, consider the lower triangular system

```{math}
  \begin{bmatrix}
    4 & 0 & 0 & 0 \\
    3 & -1 & 0 & 0 \\
    -1 & 0 & 3 & 0 \\
    1 & -1 & -1 & 2
  \end{bmatrix} \mathbf{x} =
  \begin{bmatrix}
    8 \\ 5 \\ 0 \\ 1
  \end{bmatrix}.
```

The first row of this system states simply that $4x_1=8$, which is easily solved as $x_1=8/4=2$. Now, the second row states that $3x_1-x_2=5$. As $x_1$ is already known, it can be replaced to find that $x_2 = -(5-3\cdot 2)=1$. Similarly, the third row gives $x_3=(0+1\cdot 2)/3 = 2/3$, and the last row yields $x_4=(1-1\cdot 2 + 1\cdot 1 + 1\cdot 2/3)/2 = 1/3$. Hence, the solution is

```{math}
  \mathbf{x} =
  \begin{bmatrix} 2 \\ 1 \\ 2/3 \\ 1/3
  \end{bmatrix}.
```

```{index} forward substitution
```

The process just described is called {term}`forward substitution`. In the $4\times 4$ lower triangular case of $\mathbf{L}\mathbf{x}=\mathbf{b}$ it leads to the formulas

```{math}
:label: forwardsub
\begin{split}
  x_1 &= \frac{b_1}{L_{11}}, \\
  x_2 &= \frac{b_2 - L_{21}x_1}{L_{22}}, \\
  x_3 &= \frac{b_3 - L_{31}x_1 - L_{32}x_2}{L_{33}}, \\
  x_4 &= \frac{b_4 - L_{41}x_1 - L_{42}x_2 - L_{43}x_3}{L_{44}}.
\end{split}
```

```{index} backward substitution
```

For upper triangular systems $\mathbf{U}\mathbf{x}=\mathbf{b}$ an analogous process of {term}`backward substitution` begins by solving for the last component $x_n=b_n/U_{nn}$ and working backward. For the $4\times 4$ case we have

```{math}
  \begin{bmatrix}
    U_{11} & U_{12} & U_{13} & U_{14} \\
    0 & U_{22} & U_{23} & U_{24} \\
    0 & 0 & U_{33} & U_{34} \\
    0 & 0 & 0 & U_{44}
  \end{bmatrix} \mathbf{x} =
  \begin{bmatrix}
    b_1 \\ b_2 \\ b_3 \\ b_4
  \end{bmatrix}.
```

Solving the system backward, starting with $x_4$ first and then proceeding in descending order, gives

```{math}
:label: eq-backsub
\begin{split}
  x_4 &= \frac{b_4}{U_{44}}, \\
  x_3 &= \frac{b_3 - U_{34}x_4}{U_{33}}, \\
  x_2 &= \frac{b_2 - U_{23}x_3 - U_{24}x_4}{U_{22}}, \\
  x_1 &= \frac{b_1 - U_{12}x_2 - U_{13}x_3 - U_{14}x_4}{U_{11}}.
\end{split}
```

It should be clear that forward or backward substitution fails if and only if one of the diagonal entries of the system matrix is zero. We have essentially proved the following theorem.

```{prf:theorem} Triangular singularity
:label: theorem-triangle-invert

A triangular matrix is singular if and only if at least one of its diagonal elements is zero.
```

## Implementation

Consider how to implement the sequential process implied by Equation {eq}`forwardsub`. It seems clear that we want to loop through the elements of $\mathbf{x}$ in order. Within each iteration of that loop, we have an expression whose length depends on the iteration number. This leads to a nested loop structure.

``````{prf:algorithm} forwardsub
:label: function-forwardsub

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #function-forwardsub-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #function-forwardsub-julia
:::
```` 


````{tab-item} Python
:sync: python
:::{embed} #function-forwardsub-python
:::
````
`````
``````

The implementation of backward substitution is much like forward substitution and is given in {numref}`Function {number} <function-backsub>`.

``````{prf:algorithm} backsub
:label: function-backsub

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #function-backsub-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #function-backsub-julia
:::
```` 


````{tab-item} Python
:sync: python
:::{embed} #function-backsub-python
:::
````
`````
``````

::::{prf:example} Triangular systems of equations
:label: demo-systems-triangular

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-systems-triangular-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-systems-triangular-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-systems-triangular-python
:::
```` 
`````

::::

The example in @demo-systems-triangular is our first clue that linear system problems may have large condition numbers, making inaccurate solutions inevitable in floating-point arithmetic. We will learn how to spot such problems in {numref}`section-linsys-condition-number`. Before reaching that point, however, we need to discuss how to solve general linear systems, not just triangular ones.

## Exercises

``````{exercise}
:label: problem-linearsystems-singular
✍ Find a vector $\mathbf{b}$ such that the system $\begin{bmatrix} 0&1\\0&0 \end{bmatrix} \mathbf{x}=\mathbf{b}$ has no solution.
``````

``````{exercise}
:label: problem-linearsystems-triangular
✍ Solve the following triangular systems by hand.

**(a)** $\displaystyle \begin{aligned}
-2x_1  &= -4 \\
x_1  - x_2        &= 2 \\
3x_1 + 2x_2  + x_3 &= 1
\end{aligned} \quad$
**(b)** $\displaystyle \begin{bmatrix}
4 & 0 & 0 & 0 \\
1 & -2 & 0 & 0 \\
-1 & 4 & 4 & 0 \\
2 & -5 & 5 & 1
\end{bmatrix} \mathbf{x} = \begin{bmatrix}
-4 \\ 1 \\ -3 \\ 5
\end{bmatrix}\quad$
**(c)** $\displaystyle \begin{aligned}
3x_1 +  2x_2  +  x_3      &= 1 \\
x_2   -  x_3      &= 2 \\
2 x_3    &= -4
\end{aligned}$
``````

``````{exercise}
:label: problem-linearsystems-triangularcomp
⌨ Use {numref}`Function {number} <function-forwardsub>` or {numref}`Function {number} <function-backsub>` to solve each system for $\mathbf{x}$ from the preceding exercise. Verify that the solution is correct by computing $\mathbf{b} - \mathbf{A}\mathbf{x}$.
``````

``````{exercise}
:label: problem-linearsystems-backsub
⌨  Use {numref}`Function {number} <function-backsub>` to solve the following systems.  Verify that the solution is correct by computing $\mathbf{b} - \mathbf{U}\mathbf{x}$.

**(a)** $\;\begin{bmatrix}
3 & 1 & 0  \\
0 & -1 & -2  \\
0 & 0 & 3  \\
\end{bmatrix} \mathbf{x} = \begin{bmatrix}
1 \\ 1 \\ 6
\end{bmatrix}\qquad$
**(b)** $\;\begin{bmatrix}
3 & 1 & 0 & 6 \\
0 & -1 & -2 & 7 \\
0 & 0 & 3 & 4 \\
0 & 0 & 0 & 5
\end{bmatrix} \mathbf{x} = \begin{bmatrix}
4 \\ 1 \\ 1 \\ 5
\end{bmatrix}$
``````

``````{exercise}
:label: problem-linearsystems-lumpstring
Suppose a string is stretched with tension $\tau$ horizontally between two anchors at $x=0$ and $x=1$. At each of the $n-1$ equally spaced positions $x_k=k/n$, $k=1,\ldots,n-1$, we attach a little mass $m_i$ and allow the string to come to equilibrium. This causes vertical displacement of the string. Let $q_k$ be the amount of displacement at $x_k$. If the displacements are not too large, then an approximate force balance equation is

```{math}
:numbered: false
n \tau (q_k - q_{k-1}) + n\tau (q_k - q_{k+1}) =
m_k g, \qquad k=1,\ldots,n-1,
```

where $g=-9.8$ m/s$^2$ is the acceleration due to gravity, and we define $q_0=0$ and $q_n=0$ due to the anchors. This defines a linear system for $q_1,\ldots,q_{n-1}$.

**(a)** ✍ Show that the force balance equations can be written as a linear system $\mathbf{A}\mathbf{q}=\mathbf{f}$, where $\mathbf{q}$ is a vector of the unknown displacements and $\mathbf{A}$ is a tridiagonal matrix (i.e., $A_{ij}=0$ if $|i-j|>1$) of size $(n-1)\times(n-1)$.

**(b)** ⌨  Let $\tau=10$ N, and $m_k=(1/10n)$ kg for every $k$. Using backslash, find the displacements when $n=8$ and $n=40$, and superimpose plots of $\mathbf{q}$ over $0\le x \le 1$ for the two cases. (Be sure to include the zero values at $x=0$ and $x=1$ in your plots.)

**(c)** ⌨  Repeat (b) for the case $m_k = (k/5n^2)$ kg.
``````

``````{exercise}
:label: problem-linearsystems-inverse
⌨  If $\mathbf{B}\in\mathbb{R}^{n \times p}$ has columns $\mathbf{b}_1,\ldots,\mathbf{b}_p$, then we can pose $p$ linear systems at once by writing $\mathbf{A} \mathbf{X} = \mathbf{B}$, where $\mathbf{X}$ is $n\times p$. Specifically, this equation implies $\mathbf{A} \mathbf{x}_j = \mathbf{b}_j$ for $j=1,\ldots,p$.

**(a)** Modify {numref}`Function {number} <function-forwardsub>` and {numref}`Function {number} <function-backsub>` so that they solve the case where the second input is $n\times p$ for $p\ge 1$.

**(b)** If $\mathbf{A} \mathbf{X}=\mathbf{I}$, then $\mathbf{X}=\mathbf{A}^{-1}$. Use this fact to write a function `ltinverse` that uses your modified **forwardsub** to compute the inverse of a lower triangular matrix. Test your function on at least two nontrivial matrices. (We remind you here that this is just an exercise; matrix inverses are rarely a good idea in numerical practice!)
``````

``````{exercise}
:label: problem-linearsystems-triangillcond
⌨ @demo-systems-triangular showed solutions of $\mathbf{A}\mathbf{x}=\mathbf{b}$, where

```{math}
:numbered: false
\mathbf{A} = \begin{bmatrix} 1 & -1 & 0 & \alpha-\beta & \beta \\ 0 & 1 & -1 &
0 & 0 \\ 0 & 0 & 1 & -1 & 0 \\ 0 & 0 & 0 & 1 & -1  \\ 0 & 0 & 0 & 0 & 1
\end{bmatrix}, \quad
\mathbf{b} = \begin{bmatrix} \alpha \\ 0 \\ 0 \\ 0 \\ 1 \end{bmatrix}.
```

Use {numref}`Function {number} <function-backsub>` to solve with $\alpha=0.1$ and $\beta=10,100,10^3,\ldots,10^{12}$, tabulating the values of $\beta$ and $|x_1-1|$. (This kind of behavior is explained in {numref}`section-linsys-condition-number`.)
``````
