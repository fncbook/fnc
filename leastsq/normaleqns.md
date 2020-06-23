# The normal equations

```{index} least squares
```

We seek to solve the general linear least squares problem: Given $\mathbf{A}\in\mathbb{R}^{m \times n}$ and $\mathbf{b}\in\mathbb{R}^m$, with $m>n$, find $\mathbf{x}\in\mathbb{R}^n$ such that $\| \mathbf{b} - \mathbf{A}\mathbf{x} \|_2$ is minimized. There is a concise explicit solution to the problem.

In the following proof we make use of the elementary algebraic fact that for two vectors $\mathbf{u}$ and $\mathbf{v}$,

```{math}
  (\mathbf{u}+\mathbf{v})^T(\mathbf{u}+\mathbf{v}) = \mathbf{u}^T\mathbf{u} + \mathbf{u}^T\mathbf{v} + \mathbf{v}^T\mathbf{u}
  + \mathbf{v}^T\mathbf{v} = \mathbf{u}^T\mathbf{u} + 2\mathbf{v}^T\mathbf{u} + \mathbf{v}^T\mathbf{v}.
```

```{index} normal equations
```

(theorem-normaleqns)=

````{proof:theorem}
  If $\mathbf{x}$ satisfies $\mathbf{A}^T(\mathbf{A}\mathbf{x}-\mathbf{b})=\boldsymbol{0}$, then $\mathbf{x}$ solves the linear least squares problem, i.e., $\mathbf{x}$ minimizes $\| \mathbf{b}-\mathbf{A}\mathbf{x} \|_2$.
````

````{proof:proof}
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

The condition $\mathbf{A}^T(\mathbf{A}\mathbf{x}-\mathbf{b})=\boldsymbol{0}$ is often written as

```{math}
:label: normaleqns
\mathbf{A}^T\mathbf{A}\mathbf{x}=\mathbf{A}^T\mathbf{b},
```

called the {term}`normal equations`. They have a straightforward geometric interpretation, as shown in {numref}`fig-normaleqns2d`. The vector in the range (column space) of $\mathbf{A}$ that lies closest to $\mathbf{b}$ makes the vector difference $\mathbf{A}\mathbf{x}-\mathbf{b}$ perpendicular to the range. Thus for any $\mathbf{z}$, we must have $(\mathbf{A} \mathbf{z})^T(\mathbf{A}\mathbf{x}-\mathbf{b})=0$, which is satisfied if $\mathbf{A}^T(\mathbf{A}\mathbf{x}-\mathbf{b})=\boldsymbol{0}$.

```{figure} figures/normaleqns2d.svg
:name: fig-normaleqns2d
Geometry of the normal equations.
```

```{margin}
The normal equations express an $m\times n$ linear least squares problem as an $n\times n$ linear system of equations.
```

If we group the left-hand side of the normal equations as $(\mathbf{A}^T\mathbf{A})\,\mathbf{x}$, we recognize {eq}`normaleqns` as a *square* $n\times n$ linear system to solve for $\mathbf{x}$.

## Pseudoinverse and definiteness

The $n\times m$ matrix

```{math}
:label: pinv
\mathbf{A}^+ = (\mathbf{A}^T\mathbf{A})^{-1}\mathbf{A}^T
```

```{margin}
In Julia, backslash is equivalent mathematically to left-multiplication by the inverse (square case) or pseudoinverse (rectangular case) of a matrix.
```

```{index} pseudoinverse
```

is called the {term}`pseudoinverse` of $\mathbf{A}$. Mathematically, the overdetermined least squares problem $\mathbf{A}\mathbf{x}\approx \mathbf{b}$ has the solution $\mathbf{x}=\mathbf{A}^+\mathbf{b}$. Hence we can generalize our earlier observation: backslash is equivalent mathematically to left-multiplication by the inverse (square case) or pseudoinverse (rectangular case) of a matrix. One may also compute the pseudoinverse directly using {term}`pinv`, but as with matrix inverses, this is rarely necessary in practice.

The matrix $\mathbf{A}^T\mathbf{A}$ appearing in the pseudoinverse has some important properties.

(theorem-ATA)=

````{proof:theorem} (AtA)
  For any real $m\times n$ matrix $\mathbf{A}$ with $m\ge n$, the following are true:

1. $\mathbf{A}^T\mathbf{A}$ is symmetric. 
2. $\mathbf{A}^T \mathbf{A}$ is singular if and only if the columns of $\mathbf{A}$ are linearly dependent. (Equivalently, the rank of $\mathbf{A}$ is less than $n$.)
3. If $\mathbf{A}^T\mathbf{A}$ is nonsingular, then it is positive definite. 
````

````{proof:proof}
The first part is left as an [exercise](problem-normalsymmetry). For the second part, suppose that $\mathbf{A}^T\mathbf{A}\mathbf{z}=\boldsymbol{0}$. Note that $\mathbf{A}^T\mathbf{A}$ is singular if and only if $\mathbf{z}$ may be nonzero. Left-multiplying by $\mathbf{z}^T$, we find that
  
```{math}
0 = \mathbf{z}^T\mathbf{A}^T\mathbf{A}\mathbf{z}=(\mathbf{A}\mathbf{z})^T(\mathbf{A}\mathbf{z}) = \| \mathbf{A}\mathbf{z} \|_2^2,
```

which is equivalent to $\mathbf{A}\mathbf{z}=\boldsymbol{0}$. Then $\mathbf{z}$ may be nonzero if and only if the columns of $\mathbf{A}$ are linearly dependent.

Finally, we can repeat the manipulations above to show that for any nonzero $n$-vector $\mathbf{v}$, $\mathbf{v}^T(\mathbf{A}^T\mathbf{A})\mathbf{v}=\| \mathbf{A}\mathbf{v} \|_2^2\ge 0$, and equality is not possible thanks to the second part of the theorem.
````

The definition of the pseudoinverse involves taking the inverse of a matrix and is therefore not advisable to use computationally. Instead, we simply use the definition of the normal equations to set up a linear system, which we already know how to solve. In summary, the steps for solving the linear least squares problem $\mathbf{A}\mathbf{x}\approx\mathbf{b}$ are:

1. Compute $\mathbf{N}=\mathbf{A}^T\mathbf{A}$.
1. Compute $\mathbf{z} = \mathbf{A}^T\mathbf{b}$.
1. Solve the $n\times n$ linear system $\mathbf{N}\mathbf{x} = \mathbf{z}$ for $\mathbf{x}$.

In the last step we can exploit the fact, proved in [the AtA theorem](theorem-ATA), that $\mathbf{N}$ is symmetric and positive definite, and use Cholesky factorization as in {ref}`sec-SPD`. (The backslash command does this automatically.)

```{index} normal equations
```

(function-lsnormal)=

```{proof:function} lsnormal

```{code-block} julia
:lineno-start: 1
"""
lsnormal(A,b)

Solve a linear least squares problem by the normal equations.
Returns the minimizer of ||b-Ax||.
"""
function lsnormal(A,b)

N = A'*A;  z = A'*b;
R = cholesky(N).U
w = forwardsub(R',z)                   # solve R'z=c
x = backsub(R,w)                       # solve Rx=z

return x
end
```

## Conditioning and stability

We have already used `A\b` as the native way to solve the linear least squares problem $\mathbf{A}\mathbf{x}\approx\mathbf{b}$ in Julia. The algorithm employed by the backslash does *not* proceed through the normal equations, because of instability.

The conditioning of the linear least-squares problem relates changes in the solution $\mathbf{x}$ to those in the data, $\mathbf{A}$ and $\mathbf{b}$. A full accounting of the condition number is too messy to present here, but we can generalize from the linear system problem $\mathbf{A} \mathbf{x} =\mathbf{b}$, where $m=n$. Recall that the condition number of solving $\mathbf{A} \mathbf{x}=\mathbf{b}$ is $\kappa(\mathbf{A})=\|\mathbf{A}\| \cdot \|\mathbf{A}^{-1}\|$.

Provided that the residual norm $\|\mathbf{b}-\mathbf{A}\mathbf{x}\|$ at the least-squares solution is relatively small, the conditioning of the linear least squares problem is similar. The condition number of the problem generalizes via the pseudoinverse:

```{index} condition number; of linear least squares
```

```{math}
:label: rectcond
\kappa(\mathbf{A}) = \|\mathbf{A}\| \cdot \|\mathbf{A}^{+}\|.
```

These are rectangular matrices, but the induced matrix norm is defined by {eq}`matrixnorm` just as in the square case. If $\mathbf{A}$ has rank less than $n$, then $\kappa(\mathbf{A})=\infty$. The Julia function `cond` computes condition numbers of rectangular matrices in the 2-norm.

```{margin}
When the normal equations are solved, perturbations to the data can be amplified by a factor $\kappa(\mathbf{A}^T\mathbf{A})$.
```

As an algorithm, the normal equations begin by computing the $n\times n$ system
$(\mathbf{A}^T\mathbf{A})\mathbf{x} = \mathbf{A}^T \mathbf{b}$. When these equations are solved, perturbations to the data can be amplified by a factor $\kappa(\mathbf{A}^T\mathbf{A})$. It turns out that

```{index} condition number; of a matrix
```

```{index} condition number; of normal equations
```

```{math}
:label: condATA
\kappa(\mathbf{A}^T\mathbf{A}) = \kappa(\mathbf{A})^2
```

````{sidebar} Demo
:class: demo
{doc}`demos/normaleqns-instab`
````

in the 2-norm. If $\kappa(\mathbf{A})$ is large, the squaring of it can destabilize the normal equations: while the solution of the least squares problem is sensitive, finding it via the normal equations makes it doubly so.


## Exercises

1. ✍ Work out the least squares solution when
  
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

    (problem-pinvcompute)=
2. ✍ Find the pseudoinverse $\mathbf{A}^+$ of the matrix $\mathbf{A}=\begin{bmatrix}1&-2&3\end{bmatrix}^T$.

      <!-- A = [1 -2 3]';
      ATA = A'*A
      ATAinv = inv(ATA)
      ATAinvAT = ATA\A' -->

    (problem-normalsymmetry)=
3. ✍ Prove the first statement of [the AtA theorem](theorem-ATA): $\mathbf{A}^T\mathbf{A}$ is symmetric for any $m\times n$ matrix $\mathbf{A}$ with $m \ge n$.

    (problem-pinveqinv)=
4. ✍ Prove that if $\mathbf{A}$ is a nonsingular square matrix, then $\mathbf{A}^+=\mathbf{A}^{-1}$.

5. **(a)** ✍ Show that for any $m\times n$ $\mathbf{A}$ ($m>n$) for which $\mathbf{A}^T\mathbf{A}$ is nonsingular, $\mathbf{A}^+\mathbf{A}$ is the $n\times n$ identity.
  
    **(b)** ⌨ Show using an example in Julia that $\mathbf{A}\mathbf{A}^+$ is not an identity matrix. (This matrix has rank no greater than $n$, so it can't be an $m\times m$ identity.)
  
6. ✍ Prove that the vector $\mathbf{A}\mathbf{A}^+\mathbf{b}$ is the vector in the column space (i.e., range) of $\mathbf{A}$ that is closest to $\mathbf{b}$, in the sense of the 2-norm.

7. ✍ Show that the flop count for {ref}`function-lsnormal` is asymptotically $\sim 2m n^2 + \tfrac{1}{3}n^3$. (In finding the asymptotic count you can ignore terms like $mn$ whose total degree is less than three.)

8. ⌨ Let $t_1,\ldots,t_m$ be $m+1$ equally spaced points in $[0,2\pi]$.
  
    **(a)** Let $\mathbf{A}_\beta$ be the matrix in {eq}`vandersystemrect` that corresponds to fitting data with the function $c_1 + c_2 \sin(t) + c_3 \cos(\beta t)$. Using the identity {eq}`condATA`, make a table of the condition numbers of $\mathbf{A}_\beta$ for $\beta = 2,1.1,1.01,\ldots,1+10^{-8}$.

    **(b)** Repeat part (a) using the fitting function $c_1 + c_2 \sin^2(t) + c_3 \cos^2(\beta t).$

    **(c)** Why does it make sense that $\kappa\bigl(\mathbf{A}_\beta\bigr)\to \infty$ as $\beta\to 1$ in part (b) but not in part (a)?
  
    <!-- t = linspace(0,2*pi,401)';
    A = [ t.^0 sin(t).^2 cos((1+1e-16)*t).^2 ];
    kappa = sqrt(cond(A'*A)) -->

9. ✍ ⌨  When $\mathbf{A}$ is $m\times n$ with rank less than $n$, the pseudoinverse is still defined and can be computed using `pinv` from `LinearAlgebra`. However, the behavior in this case is not always intuitive. Let
  
    ```{math}
    \mathbf{A}_s =
    \begin{bmatrix}
      1 & 1 \\ 0 & 0 \\ 0 & s
    \end{bmatrix}.
    ```

    Then $\mathbf{A}_0$ has rank equal to one. Demonstrate experimentally that $\mathbf{A}_0^+\neq \lim_{s\to 0} \mathbf{A}_s^+$.
