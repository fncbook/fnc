---
numbering:
  enumerator: 3.2.%s
---
(section-leastsq-normaleqns)=
# The normal equations

```{index} linear least-squares problem
```

We now solve the general linear least-squares problem in {numref}`Definition {number} <definition-fitting-linearls>`. That is, given $\mathbf{A}\in\mathbb{R}^{m \times n}$ and $\mathbf{b}\in\mathbb{R}^m$, with $m>n$, find the $\mathbf{x}\in\mathbb{R}^n$ that minimizes $\| \mathbf{b} - \mathbf{A}\mathbf{x} \|_2$. 

There is a concise explicit solution. In the following proof we make use of the elementary algebraic fact that for two vectors $\mathbf{u}$ and $\mathbf{v}$,

```{math}
  (\mathbf{u}+\mathbf{v})^T(\mathbf{u}+\mathbf{v}) = \mathbf{u}^T\mathbf{u} + \mathbf{u}^T\mathbf{v} + \mathbf{v}^T\mathbf{u}
  + \mathbf{v}^T\mathbf{v} = \mathbf{u}^T\mathbf{u} + 2\mathbf{v}^T\mathbf{u} + \mathbf{v}^T\mathbf{v}.
```

```{index} ! normal equations
```

(theorem-normaleqns)=
````{prf:theorem}

If $\mathbf{x}$ satisfies $\mathbf{A}^T(\mathbf{A}\mathbf{x}-\mathbf{b})=\boldsymbol{0}$, then $\mathbf{x}$ solves the linear least-squares problem, i.e., $\mathbf{x}$ minimizes $\| \mathbf{b}-\mathbf{A}\mathbf{x} \|_2$.
````

````{prf:proof}

Let $\mathbf{y}\in \mathbb{R}^n$ be any vector. Then
  $\mathbf{A}(\mathbf{x}+\mathbf{y})-\mathbf{b}=\mathbf{A}\mathbf{x}-\mathbf{b}+\mathbf{A}\mathbf{y}$, and
  
```{math}
\begin{split}
    \| \mathbf{A}(\mathbf{x}+\mathbf{y})-\mathbf{b} \|_2^2 &=
    [(\mathbf{A}\mathbf{x}-\mathbf{b})+(\mathbf{A}\mathbf{y})]^T[(\mathbf{A}\mathbf{x}-\mathbf{b})+(\mathbf{A}\mathbf{y})]\\
    &= (\mathbf{A}\mathbf{x}-\mathbf{b})^T(\mathbf{A}\mathbf{x}-\mathbf{b}) + 2(\mathbf{A}\mathbf{y})^T(\mathbf{A}\mathbf{x}-\mathbf{b}) + (\mathbf{A}\mathbf{y})^T(\mathbf{A}\mathbf{y})\\
    &= \| \mathbf{A}\mathbf{x}-\mathbf{b} \|_2^2 + 2\mathbf{y}^T \mathbf{A}^T(\mathbf{A}\mathbf{x}-\mathbf{b}) + \| \mathbf{A}\mathbf{y} \|_2^2 \\
    &= \| \mathbf{A}\mathbf{x}-\mathbf{b} \|_2^2 + \| \mathbf{A}\mathbf{y} \|_2^2 \\
    & \ge \| \mathbf{A}\mathbf{x}-\mathbf{b} \|_2^2.
  \end{split}
```
````

::::{prf:definition} Normal equations
Given $\mathbf{A}\in \real^{m\times n}$ and $\mathbf{b}\in \real^{m}$, the **normal equations** for the linear least-squares problem $\operatorname{argmin}\| \mathbf{b}- \mathbf{A} \mathbf{x}\|$ are $\mathbf{A}^T(\mathbf{A}\mathbf{x}-\mathbf{b})=\boldsymbol{0}$, or equivalently,

```{math}
:label: normaleqns
\mathbf{A}^T\mathbf{A}\mathbf{x}=\mathbf{A}^T\mathbf{b}.
```
::::

The normal equations have a geometric interpretation, as shown in {numref}`fig-normaleqns2d`. The vector in the range (column space) of $\mathbf{A}$ that lies closest to $\mathbf{b}$ makes the vector difference $\mathbf{A}\mathbf{x}-\mathbf{b}$ perpendicular to the range. Thus for any $\mathbf{z}$, we must have $(\mathbf{A} \mathbf{z})^T(\mathbf{A}\mathbf{x}-\mathbf{b})=0$, which is satisfied if $\mathbf{A}^T(\mathbf{A}\mathbf{x}-\mathbf{b})=\boldsymbol{0}$.

```{figure} figures/normaleqns2d.svg
:name: fig-normaleqns2d
Geometry of the normal equations. The smallest residual is orthogonal to the range of the matrix $\mathbf{A}$. 
```

## Pseudoinverse and definiteness

If we associate the left-hand side of the normal equations as $(\mathbf{A}^T\mathbf{A})\,\mathbf{x}$, we recognize {eq}`normaleqns` as a *square* $n\times n$ linear system to solve for $\mathbf{x}$. 

```{index} ! pseudoinverse
```

::::{prf:definition} Pseudoinverse
If $\mathbf{A}\in\real^{m\times n}$  with $m>n$, its {term}`pseudoinverse` is the $n\times m$ matrix

```{math}
:label: pinv
\mathbf{A}^+ = (\mathbf{A}^T\mathbf{A})^{-1}\mathbf{A}^T.
```
::::

Mathematically, the overdetermined least-squares problem $\mathbf{A}\mathbf{x}\approx \mathbf{b}$ has the solution $\mathbf{x}=\mathbf{A}^+\mathbf{b}$. 

Computationally we can generalize the observation about Julia from Chapter 2: backslash is equivalent mathematically to left-multiplication by the inverse (in the square case) or pseudoinverse (in the rectangular case) of a matrix. One can also compute the pseudoinverse directly using `pinv`, but as with matrix inverses, this is rarely necessary in practice.

The matrix $\mathbf{A}^T\mathbf{A}$ appearing in the pseudoinverse has some important properties.

```{index} normal equations
```

(theorem-ATA)=
::::{prf:theorem} 

For any real $m\times n$ matrix $\mathbf{A}$ with $m\ge n$, the following are true:

1. $\mathbf{A}^T\mathbf{A}$ is symmetric. 
2. $\mathbf{A}^T \mathbf{A}$ is singular if and only if the columns of $\mathbf{A}$ are linearly dependent. (Equivalently, if and only if the rank of $\mathbf{A}$ is less than $n$.)
3. If $\mathbf{A}^T\mathbf{A}$ is nonsingular, then it is positive definite. 
::::

````{prf:proof}
The first part is left as [Exercise 3](#problem-normaleqns-symmetry). For the second part, suppose that $\mathbf{A}^T\mathbf{A}\mathbf{z}=\boldsymbol{0}$. Note that $\mathbf{A}^T\mathbf{A}$ is singular if and only if $\mathbf{z}$ may be nonzero. Left-multiplying by $\mathbf{z}^T$, we find that
  
```{math}
0 = \mathbf{z}^T\mathbf{A}^T\mathbf{A}\mathbf{z}=(\mathbf{A}\mathbf{z})^T(\mathbf{A}\mathbf{z}) = \| \mathbf{A}\mathbf{z} \|_2^2,
```

which is equivalent to $\mathbf{A}\mathbf{z}=\boldsymbol{0}$. Then $\mathbf{z}$ may be nonzero if and only if the columns of $\mathbf{A}$ are linearly dependent.

Finally, we can repeat the manipulations above to show that for any nonzero $n$-vector $\mathbf{v}$, $\mathbf{v}^T(\mathbf{A}^T\mathbf{A})\mathbf{v}=\| \mathbf{A}\mathbf{v} \|_2^2\ge 0$, and equality is not possible thanks to the second part of the theorem.
````

## Implementation

The definition of the pseudoinverse involves taking the inverse of a matrix, so it is not advisable to use the pseudoinverse computationally. Instead, we use the definition of the normal equations to set up a linear system, which we already know how to solve. In summary, the steps for solving the linear least squares problem $\mathbf{A}\mathbf{x}\approx\mathbf{b}$ are:

(algorithm-normaleqns-solve)=
::::{prf:algorithm} Solution of linear least squares by the normal equations
1. Compute $\mathbf{N}=\mathbf{A}^T\mathbf{A}$.
2. Compute $\mathbf{z} = \mathbf{A}^T\mathbf{b}$.
3. Solve the $n\times n$ linear system $\mathbf{N}\mathbf{x} = \mathbf{z}$ for $\mathbf{x}$.
::::

Steps 1 and 3 of {numref}`Algorithm {number} <algorithm-normaleqns-solve>` dominate the flop count.

In the last step we can exploit the fact, proved in {numref}`Theorem %s <theorem-ATA>`, that $\mathbf{N}$ is symmetric and positive definite, and use Cholesky factorization as in {numref}`section-linsys-structure`. This detail is included in {numref}`Function {number} <function-lsnormal>`.

(function-lsnormal)=
``````{prf:algorithm} lsnormal
`````{tab-set} 
````{tab-item} Julia
:sync: julia
:::{embed} #function-lsnormal-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #function-lsnormal-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #function-lsnormal-python
:::
````
`````
``````

(theorem-normaleqns-flops)=
```{prf:theorem}
Solution of linear least squares by the normal equations takes $\sim (mn^2 + \frac{1}{3}n^3)$ flops.
```

## Conditioning and stability

We have already used `A\b` as the native way to solve the linear least-squares problem $\mathbf{A}\mathbf{x}\approx\mathbf{b}$ in Julia. The algorithm employed by the backslash does *not* proceed through the normal equations, because of instability.

The conditioning of the linear least-squares problem relates changes in the solution $\mathbf{x}$ to those in the data, $\mathbf{A}$ and $\mathbf{b}$. A full accounting of the condition number is too messy to present here, but we can get the main idea. We start by generalizing our previous definition of the matrix condition number.

```{index} condition number; of linear least squares
```

```{index} condition number; of a matrix
```

````{prf:definition} Matrix condition number (rectangular case)
If $\mathbf{A}$ is $m\times n$ with $m > n$, then its condition number is defined to be

```{math}
:label: rectcond
\kappa(\mathbf{A}) = \|\mathbf{A}\|_2 \cdot \|\mathbf{A}^{+}\|_2.
```

If the rank of $ \mathbf{A}$ is less than $n$ (i.e., if it has linearly dependent columns), then $\kappa(\mathbf{A})=\infty$.
````

Provided that the minimum residual norm $\|\mathbf{b}-\mathbf{A}\mathbf{x}\|$ is relatively small, the conditioning of the linear least-squares problem is close to $\kappa(\mathbf{A})$. 

As an algorithm, the normal equations begin by computing the data for the $n\times n$ system $(\mathbf{A}^T\mathbf{A})\mathbf{x} = \mathbf{A}^T \mathbf{b}$. When these equations are solved, perturbations to the data can be amplified by a factor $\kappa(\mathbf{A}^T\mathbf{A})$. 

```{index} condition number; of normal equations
```

The following can be proved using results in Chapter 7.

````{prf:theorem} Condition number in the normal equations
If $\mathbf{A}$ is $m\times n$ with $m > n$, then

```{math}
:label: condATA
\kappa(\mathbf{A}^T\mathbf{A}) = \kappa(\mathbf{A})^2.
```
````
This squaring of the condition number in the normal equations is the cause of instability. If $\kappa(\mathbf{A})$ is large, the squaring of it can destabilize the normal equations: while the solution of the least-squares problem is sensitive, finding it via the normal equations makes it doubly so.


(demo-normaleqns-instab)=
::::{prf:example} Instability in the normal equations
`````{tab-set} 
````{tab-item} Julia
:sync: julia
:::{embed} #demo-normaleqns-instab-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-normaleqns-instab-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-normaleqns-instab-python
:::
```` 
`````
::::

## Exercises

1. ✍ Work out the least-squares solution when
  
    ```{math}
    \mathbf{A} = \begin{bmatrix}
      2 & -1 \\
      0 & 1 \\
      -2 & 2
    \end{bmatrix}, \qquad \mathbf{b} =
    \begin{bmatrix}
      1\\-5\\6
    \end{bmatrix}.
    ```

(problem-normaleqns-pinvcompute)=
2. ✍ Use {eq}`pinv` to find the pseudoinverse $\mathbf{A}^+$ of the matrix $\mathbf{A}=\begin{bmatrix}1&-2&3\end{bmatrix}^T$.

(problem-normaleqns-symmetry)=
3. ✍ Prove the first statement of {numref}`Theorem %s <theorem-ATA>`: $\mathbf{A}^T\mathbf{A}$ is symmetric for any $m\times n$ matrix $\mathbf{A}$ with $m \ge n$.

(problem-normaleqns-pinveqinv)=
4. ✍ Prove that if $\mathbf{A}$ is an invertible square matrix, then $\mathbf{A}^+=\mathbf{A}^{-1}$.

5. **(a)** ✍ Show that for any $m\times n$ $\mathbf{A}$ with $m>n$ for which $\mathbf{A}^T\mathbf{A}$ is nonsingular, $\mathbf{A}^+\mathbf{A}$ is the $n\times n$ identity.
  
    **(b)** ⌨ Show using an example in Julia that $\mathbf{A}\mathbf{A}^+$ is not an identity matrix. (This matrix has rank no greater than $n$, so it can't be an $m\times m$ identity.)
  
6. ✍ Prove that the vector $\mathbf{A}\mathbf{A}^+\mathbf{b}$ is the vector in the column space (i.e., range) of $\mathbf{A}$ that is closest to $\mathbf{b}$ in the sense of the 2-norm.

7. ✍ Show that the flop count for {numref}`Function {number} <function-lsnormal>` is asymptotically $\sim 2m n^2 + \tfrac{1}{3}n^3$. (In finding the asymptotic count you can ignore terms like $m n$ whose total degree is less than 3.)

8. ⌨ Let $t_1,\ldots,t_m$ be $m$ equally spaced points in $[0,2\pi]$. In this exercise, use $m=500$.
  
    **(a)** Let $\mathbf{A}_\beta$ be the matrix in {eq}`vandersystemrect` that corresponds to fitting data with the function $c_1 + c_2 \sin(t) + c_3 \cos(\beta t)$. Using the identity {eq}`condATA`, make a table of the condition numbers of $\mathbf{A}_\beta$ for $\beta = 2,1.1,1.01,\ldots,1+10^{-8}$.

    **(b)** Repeat part (a) using the fitting function $c_1 + c_2 \sin^2(t) + c_3 \cos^2(\beta t).$

    **(c)** Why does it make sense that $\kappa\bigl(\mathbf{A}_\beta\bigr)\to \infty$ as $\beta\to 1$ in part (b) but not in part (a)?
  
9. ✍ ⌨  When $\mathbf{A}$ is $m\times n$ with rank less than $n$, the pseudoinverse is still defined and can be computed using `pinv` from `LinearAlgebra`. However, the behavior in this case is not always intuitive. Let
  
    ```{math}
    \mathbf{A}_s =
    \begin{bmatrix}
      1 & 1 \\ 0 & 0 \\ 0 & s
    \end{bmatrix}.
    ```

    Then $\mathbf{A}_0$ has rank equal to 1. Demonstrate experimentally that $\mathbf{A}_0^+\neq \lim_{s\to 0} \mathbf{A}_s^+$.


