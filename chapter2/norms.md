---
numbering:
  enumerator: 2.7.%s
---
(section-linsys-norms)=

# Vector and matrix norms

The manipulations on matrices and vectors so far in this chapter have been algebraic, much like those in an introductory linear algebra course. In order to progress to the analysis of the algorithms we have introduced, we need a way to measure the size of vectors and matrices—size in the sense of magnitude or distance, not the number of rows and columns.

## Vector norms

```{index} ! norm; vector
```

For vectors, we use a **norm** $\| \cdot \|$, which is a function from $\real^n$ to $\real$ with the following properties for all $n$-vectors $\mathbf{x},\mathbf{y}$ and scalars $\alpha$:[^complexnorm]

[^complexnorm]: The same statements work for vectors with complex entries, with complex modulus in place of absolute values.
  
```{math}
\begin{align*}
\label{norm-properties}
\norm{\mathbf{x}} &\ge 0, \\
\norm{\mathbf{x}} &=0 \;\Leftrightarrow \; \mathbf{x}=\boldsymbol{0}, \\
\norm{\alpha \mathbf{x}} &= \abs{\alpha}\, \norm{\mathbf{x}}, \\
\norm{\mathbf{x}+\mathbf{y}} & \le \norm{\mathbf{x}} + \norm{\mathbf{y}}.
\end{align*}
```

The last of these properties is known as the **triangle inequality**. It is natural to interpret $\| \mathbf{x} \|=\| \mathbf{x}-\boldsymbol{0} \|$ as the distance from $\mathbf{x}$ to the origin and $\| \mathbf{x}-\mathbf{y} \|$ as the distance from $\mathbf{x}$ to $\mathbf{y}$. We will be using only the three most important vector norms, defined as follows.

:::{prf:definition} Common vector norms
:label: definition-vectornorms
**2-norm:** $\quad \twonorm{\mathbf{x}} = \left( \displaystyle \sum_{i=1}^n |x_i|^2 \right)^{\frac{1}{2}} = \sqrt{\rule[1mm]{0pt}{0.75em}\mathbf{x}^T \mathbf{x}}$

**$\infty$-norm** or **max-norm:** $\quad \infnorm{\mathbf{x}} = \displaystyle \max_{i=1,\dots,n} |x_i|$

**1-norm:** $\quad \onenorm{\mathbf{x}} = \displaystyle \sum_{i=1}^n |x_i|$
:::

The 2-norm corresponds to ordinary Euclidean distance.

```{index} ! unit vector
```

```{prf:definition} Unit vector
:label: definition-unitvector
A vector $\mathbf{x}$ is a {term}`unit vector` if $\| \mathbf{x} \|=1$.
```

For any nonzero vector $\mathbf{v}$, we can find a unit vector through the normalization $\mathbf{x}=\mathbf{v}/\|\mathbf{v}\|$. Thus, we can interpret

```{math}
:label: magdir
  \mathbf{v} = \| \mathbf{v} \| \,\cdot\, \frac{\mathbf{v}}{\| \mathbf{v} \|}
```

as writing a nonzero vector $\mathbf{v}$ in magnitude–direction form.

::::{prf:example} Vector norms
:label: demo-norms-vector
Given the vector $\mathbf{x}= \bigl[ 2 ,\, -3 ,\, 1 ,\, -1 \bigr]^T$, we have
\begin{align*}
    \| \mathbf{x} \|_2 &= \sqrt{ 4 + 9 + 1 + 1 } = \sqrt{15}, \\[1ex]
    \| \mathbf{x} \|_\infty &= \max\{ 2,3,1,1 \} = 3,\\[1ex]
    \| \mathbf{x} \|_1 &= 2 + 3 + 1 + 1 = 7.
\end{align*}

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-norms-vector-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-norms-vector-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-norms-vector-python
:::
```` 
`````

::::

:::{note}
Most of the time, when just $\| \mathbf{x} \|$ is written, the 2-norm is implied. However, in this section we use it to mean a generic, unspecified vector norm.
:::

We say that a sequence of vectors $\mathbf{x}_1,\mathbf{x}_2,\ldots$ **converges** to $\mathbf{x}$ if

```{math}
:label: vectorconverge
  \lim_{k\rightarrow\infty} \norm{\mathbf{x}_k - \mathbf{x}} = 0.
```

By definition, a sequence is convergent in the infinity norm if and only if it converges componentwise. The same is true for a convergent sequence in _any_ norm.

```{prf:theorem} Norm equivalence
:label: theorem-normequivalence
  In a finite-dimensional space, convergence in any norm implies convergence in all norms.
```

## Matrix norms

```{index} ! norm; Frobenius
```

Although we view matrices as two-dimensional, we can also interpret them as vectors: simply stack the columns on top of one another.[^colmajor] Hence, we can define a matrix norm via the vector 2-norm.

::::{prf:definition} Frobenius norm
:label: definition-frobenius
The {term}`Frobenius norm` of matrix $\mathbf{A}$ is defined as

```{math}
:label: frobenius
\| \mathbf{A} \|_F = \left( \sum_{i,j} |A_{ij}|^2 \right)^{1/2}.
```

::::

```{index} Julia; norm
```

<!-- This is the norm computed by the `norm` function in Julia. -->

[^colmajor]: Column stacking is actually how matrices are stored in memory within Julia and is known as **column-major order**. MATLAB and FORTRAN also use column-major order, while C and Python use row-major order, in which the rows are stacked.

However, it often proves to be more useful to define matrix norms differently.

```{index} ! norm; matrix
```

::::{prf:definition} Induced matrix norm
:label: definition-inducednorm
Given a vector norm $\| \cdot \|_p$, we define an {term}`induced matrix norm` for any $m\times n$ matrix $\mathbf{A}$ as

```{math}
:label: matrixnorm
\| \mathbf{A} \|_{p} = \max_{\| \mathbf{x} \|_p=1} \| \mathbf{A}\mathbf{x} \|_p =
\max_{\mathbf{x}\neq \boldsymbol{0}} \frac{\| \mathbf{A}\mathbf{x} \|_p}{\| \mathbf{x} \|_p}.
```

::::

The last equality above follows from linearity (as shown in @problem-norms-linearity).  It is derived from the interpretation of a matrix as a linear operator between $\real^n$ and $\real^m$. Thus in the 2-norm, for instance,

```{math}
\| \mathbf{A} \|_2 = \max_{\| \mathbf{x} \|_2=1} \| \mathbf{A}\mathbf{x} \|_2.
```

````{prf:example} Norm of the identity matrix
:label: example-norms-identity

In any induced matrix norm, the identity matrix $\mathbf{I}$ has norm $1$. This is because $\| \mathbf{I}\mathbf{x} \| = \| \mathbf{x} \|$ for all vectors $\mathbf{x}$, so the max of $\| \mathbf{I}\mathbf{x} \|$ over unit vectors is $1$. 
````

One can interpret the definition of an induced norm geometrically.  Each vector $\mathbf{x}$ on the unit "sphere" (as defined by the chosen vector norm) is mapped to its image $\mathbf{A}\mathbf{x}$, and the norm of $\mathbf{A}$ is the radius of the smallest "sphere" that encloses all such images.

The definition of an induced matrix norm may seem oddly complicated. However, there are some key properties that follow directly from the definition.

````{prf:theorem} Norm inequalities
:label: theorem-norms-inequalities
Let $\| \cdot \|$ designate a matrix norm and the vector norm that induced it. Then for all matrices and vectors of compatible sizes,

```{math}
:label: normineq1
\| \mathbf{A}\mathbf{x} \| \le \| \mathbf{A} \|\cdot \| \mathbf{x} \|.
```

For all matrices of compatible sizes,

```{math}
:label: normineq2
\| \mathbf{A}\mathbf{B} \| \le \| \mathbf{A} \|\cdot\| \mathbf{B} \|.
```

For a square matrix $\mathbf{A}$,

```{math}
:label: normineq3
\| \mathbf{A}^k \| \le \| \mathbf{A} \|^k \text{ for any integer $k\ge 0$.}
```
````

````{prf:proof}
:enumerated: false

The first result is trivial if $\mathbf{x}=\boldsymbol{0}$; otherwise,

```{math}
:enumerated: false
\frac{ \| \mathbf{A}\mathbf{x} \| }{\| \mathbf{x} \|} \le
\max_{\mathbf{x}\neq \boldsymbol{0}}  \frac{\| \mathbf{A}\mathbf{x} \|}{\| \mathbf{x} \|} = \| \mathbf{A} \|.
```

Inequality {eq}`normineq2` then follows because

```{math}
:enumerated: false
\| \mathbf{A}\mathbf{B}\mathbf{x} \| =\| \mathbf{A}(\mathbf{B}\mathbf{x}) \|\le \| \mathbf{A} \|\cdot\| \mathbf{B}\mathbf{x} \| \le
\| \mathbf{A} \|\cdot\| \mathbf{B} \|\cdot\| \mathbf{x} \|,
```

and then

```{math}
:enumerated: false
\| \mathbf{A}\mathbf{B} \| = \max_{\mathbf{x}\neq \boldsymbol{0}} \frac{\| \mathbf{A}\mathbf{B}\mathbf{x} \|}{\| \mathbf{x} \|} \le
\max_{\mathbf{x}\neq \boldsymbol{0}} \| \mathbf{A} \|\cdot\| \mathbf{B} \| = \| \mathbf{A} \|\cdot\| \mathbf{B} \|.
```

Finally,  {eq}`normineq3` results from repeated application of {eq}`normineq2`.
````

Two of the vector norms we have encountered induce matrix norms that are easy to compute from the matrix elements.

::::{prf:theorem} Matrix $\infty$-norm and 1-norm

```{math}
:label: mxinfnorm
\| \mathbf{A} \|_\infty = \max_{1\le \,i \,\le n} \sum_{j=1}^n |A_{ij}|,
```

```{math}
:label: mxonenorm
\| \mathbf{A} \|_1 = \max_{1\le \,j\, \le n} \sum_{i=1}^n |A_{ij}|.
```

::::

````{tip}
:open: false
A mnemonic for @mxinfnorm and @mxonenorm is that the $\infty$ symbol extends horizontally while the $1$ character extends vertically, each indicating the direction of the summation in its formula. Also, both formulas give the same result for $m\times 1$ matrices as the vector norm. In both cases you must take absolute values of the matrix elements before summing.
````

````{note}
There is no general simple formula like {eq}`mxinfnorm` or {eq}`mxonenorm` for the matrix 2-norm. The computation of the matrix 2-norm is discussed further in Chapter 7.
````

Despite the lack of a simple formula for it, the 2-norm is the default choice for many applications and theorems.

````{important}
The usual default matrix norm is the induced 2-norm.
````

::::{prf:example} Matrix norms
:label: demo-norms-matrix

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-norms-matrix-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-norms-matrix-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-norms-matrix-python
:::
```` 
`````

::::

## Exercises

``````{exercise}
:label: problem-norms-taxicab
✍ Why is the vector 1-norm also called the *taxicab norm*?
``````

``````{exercise}
:label: problem-norms-circle

✍ **(a)** Draw the unit "circle" in the $\infty$-norm, i.e., the set of all vectors $\mathbf{x}\in\real^2$ such that $\| \mathbf{x} \|_\infty=1$.

**(b)** Draw the unit "circle" in the 1-norm.
``````

``````{exercise}
:label: problem-norms-inequalities

✍ Prove that for all vectors $\mathbf{x}\in\real^n,$

**(a)** $\| \mathbf{x} \|_\infty \le \| \mathbf{x} \|_2; \qquad$ 
**(b)** $\| \mathbf{x} \|_2 \le \| \mathbf{x} \|_1$.
``````

``````{exercise}
:label: problem-norms-csinequality
✍ Prove that for any vectors $\mathbf{x}$, $\mathbf{y}$ in $\real^n$, $|\mathbf{x}^T\mathbf{y}| \le \norm{\mathbf{x}}_1\, \norm{\mathbf{y}}_\infty.$
``````

``````{exercise}
:label: problem-norms-linearity
✍ Prove using @definition-inducednorm that for any induced matrix norm, matrix $\mathbf{A}$, and scalar $c$, $\| c\mathbf{A} \| = |c|\cdot \| \mathbf{A} \|.$
``````

``````{exercise}
:label: problem-norms-extremal

✍ Let $\mathbf{A} =
\displaystyle \begin{bmatrix}
-1 & 1 \\ 2 & 2
\end{bmatrix}.$

**(a)** Find all vectors satisfying $\|\mathbf{x}\|_\infty=1$ and $\| \mathbf{A}\mathbf{x} \|_\infty=\| \mathbf{A} \|_\infty$.

**(b)** Find a vector satisfying $\|\mathbf{x}\|_1=1$ and $\| \mathbf{A}\mathbf{x} \|_1=\| \mathbf{A} \|_1$.

**(c)** Find a vector satisfying $\|\mathbf{x}\|_2=1$ such that $\| \mathbf{A}\mathbf{x} \|_2=\| \mathbf{A} \|_2$. (Hint: Use the definition of $\|\mathbf{A}\|_2$ as the maximum of $\|\mathbf{A}\mathbf{x}\|_2$ over unit vectors, and parameterize the unit vectors as $(\cos(\theta),\sin(\theta))$. Then $\twonorm{\mathbf{A}\mathbf{x}}^2$ is a function of $\theta$ that is easy to maximize.)
``````

``````{exercise}
:label: problem-norms-equivalence
✍ Prove the equivalence of the two formulas for a matrix norm in @matrixnorm.
``````

``````{exercise}
:label: problem-norms-inverse
✍ Prove that for any induced matrix norm and nonsingular matrix $\mathbf{A},$ $\| \mathbf{A}^{-1} \| \ge (\| \mathbf{A} \|)^{-1}.$ (Hint: Apply @theorem-norms-inequalities.)
``````

``````{exercise}
:label: problem-norms-maxelement

✍ **(a)** Prove that for any $\mathbf{v}\in \real^n,$

```{math}
\| \mathbf{v} \|_p \ge \max_{i=1,\ldots,n} |v_i|,
```

where $p=1,$ $2,$ or $\infty.$

**(b)** Prove that for any $\mathbf{A}\in\real^{n \times n},$

```{math}
\| \mathbf{A} \|_p \ge \max_{i,j=1,\ldots,n} |A_{ij}|,
```

where $p=1,$ $2,$ or $\infty.$ (Hint: For $p=2$, rearrange {eq}`normineq1` for a well-chosen particular value of $\mathbf{x}.$)
``````

``````{exercise}
:label: problem-norms-diagnorm
✍ Prove using @definition-inducednorm that if $\mathbf{D}$ is a diagonal matrix, then $\|\mathbf{D}\|_2 = \max_{i} |D_{ii}|.$ You may assume the matrix is real and square, but that does not affect the result or the proof in any significant way. (Hint: Let $M=\max_{i} |D_{ii}|.$ Proceed in two stages, showing that $\|\mathbf{D}\|_2\ge M$ and separately that $\|\mathbf{D}\|_2\le M.$)
``````

``````{exercise}
:label: problem-norms-neumann

✍ Suppose that $\mathbf{A}$ is ${n\times n}$ and that $\| \mathbf{A} \|<1$ in some induced matrix norm.

**(a)** Show that $(\mathbf{I}-\mathbf{A})$ is nonsingular. (Hint: Use the definition of an induced matrix norm to show that if $(\mathbf{I}-\mathbf{A})\mathbf{x}=\boldsymbol{0}$ for all nonzero $\mathbf{x}$, then $\| \mathbf{A} \|\ge 1$.)

**(b)** Show that $\displaystyle \lim_{m\rightarrow \infty} \mathbf{A}^m = \boldsymbol{0}$. (For matrices as with vectors, we say $\mathbf{B}_m \rightarrow \mathbf{L}$ if $\| \mathbf{B}_m-\mathbf{L} \| \rightarrow 0$.)

**(c)** Use (a) and (b) to show that we may obtain the geometric series

```{math}
(\mathbf{I}-\mathbf{A})^{-1} = \sum_{k=0}^\infty \mathbf{A}^k.
```

(Hint: Start with $\left(\displaystyle\sum_{k=0}^m \mathbf{A}^k\right)(\mathbf{I}-\mathbf{A})$ and take the limit.)
``````
