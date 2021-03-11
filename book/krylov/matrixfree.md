# Matrix-free iterations

A primary reason for our interest in matrices is their relationship to linear transformations. If we define $\mathbf{f}(\mathbf{x})=\mathbf{A}\mathbf{x}$, then for all $\mathbf{x}$, $\mathbf{y}$, and $\alpha$,

:::{math}
:label: lintrans
\begin{split}
\mathbf{f}(\mathbf{x} + \mathbf{y} ) &= \mathbf{f}(\mathbf{x}) + \mathbf{f}(\mathbf{y} ), \\
\mathbf{f}(\alpha \mathbf{x} ) & = \alpha \mathbf{f}(\mathbf{x}).
\end{split}
:::

These properties define a linear transformation. Moreover, *every* linear transformation between finite-dimensional vector spaces can be represented the same way, as a matrix-vector multiplication.

In scientific computing one often has a procedure for computing a transformation at any given input. That is, for some multidimensional function $\mathbf{f}$ you have the ability to compute $\mathbf{f}(\mathbf{x})$ at any given value of $\mathbf{x}$. In the [rootfinding chapter](../nonlineqn) we solved the rootfinding problem $\mathbf{f}(\mathbf{x})=\boldsymbol{0}$ with secant and quasi-Newton methods that needed exactly this capability. By repeatedly evaluating $\mathbf{f}$ at cleverly chosen points, these algorithms were able to return a value for $\mathbf{f}^{-1}(\boldsymbol{0})$.

A close examination reveals that Krylov subspace methods have a similar feel when the transformation is a linear map $\mathbf{f}(\mathbf{x})=\mathbf{A}\mathbf{x}$. While we presumed that we had access to $\mathbf{A}$ directly, in fact the only use of it was to compute $\mathbf{A}\mathbf{x}$ for a sequence of iteratively chosen $\mathbf{x}$ values. This detail goes back to where we started with the power iteration: what can we compute about $\mathbf{A}$ by using only matrix-vector multiplications? In the cases of GMRES, MINRES, and CG, the answer is that we can find $\mathbf{A}^{-1}\mathbf{b}$.

Bringing all of these observations together leads us to a cornerstone of modern scientific computation: *matrix-free iterations*. Krylov subspace methods can be used to invert a linear transformation if one provides code for the transformation, even if its associated matrix is not known explicitly. That is, it's possible to solve $\mathbf{A}\mathbf{x}=\mathbf{b}$ without even knowing $\mathbf{A}$! We look at an extended example of this capability in the rest of this section.
## Blurring images

```{index} matrix; as image
```
In an [earlier chapter](../matrixanaly/insight.md) we saw that a grayscale image can be represented as an $m\times n$ matrix $\mathbf{X}$ of pixel intensity values. Now consider a simple model for blurring the image. Define $\mathbf{B}$ as the $m\times m$ tridiagonal matrix

:::{math}
:label: blurmatrix
B_{ij} =
\begin{cases}
\tfrac{1}{2}, & \text{if $i=j$},\\
\tfrac{1}{4}, & \text{if $|i-j|=1$},\\
0, & \text{otherwise.}
\end{cases}
:::

Now $\mathbf{B}\mathbf{X}$ applies $\mathbf{B}$ to each column of $\mathbf{X}$. Within that column it does a weighted average of the values of each pixel and its two neighbors. That has the effect of blurring the image vertically. We can increase the amount of blur by applying $\mathbf{B}$ repeatedly.

In order to blur horizontally, we can transpose the image and apply blurring in the same way. We need a blurring matrix defined as in {eq}`blurmatrix` but with size $n\times n$. We call this matrix $\mathbf{C}$. Altogether the horizontal blurring is done by transposing, applying $\mathbf{C}$, and transposing back to the original orientation. That is,

:::{math}
\bigl(\mathbf{C} \mathbf{X}^T\bigr)^T = \mathbf{X}\mathbf{C}^T = \mathbf{X}\mathbf{C},
:::

using the symmetry of $\mathbf{C}$. So we can describe blur in both directions as the function

:::{math}
:label: blurfunction
\operatorname{blur}(\mathbf{X}) = \mathbf{B}^k \mathbf{X} \mathbf{C}^k,
:::

for a positive integer $k$.

::::{prf:example} Julia demo
:class: demo
:label: demos-matrixfree-blur
{doc}`demos/matrixfree-blur`
::::

## Deblurring

A more interesting operation is *deblurring*: given an image blurred by poor focus, can we reconstruct the true image? Conceptually, we want to invert the function $\operatorname{blur}(\mathbf{X})$.

It's easy to see from {eq}`blurfunction` that the blur operation is a linear transformation on image matrices. But an $m\times n$ image matrix is equivalent to a length-$mn$ vector—it's just a matter of interpreting the shape of the same data. Let $\operatorname{vec}(\mathbf{X})=\mathbf{x}$ and $\operatorname{unvec}(\mathbf{x})=\mathbf{X}$ be the mathematical statements of such reshaping operations. Now say $\mathbf{X}$ is the original image and $\mathbf{Z}=\operatorname{blur}(\mathbf{X})$ is the blurred one. Then by linearity there is some matrix $\mathbf{A}$ such that

:::{math}
\mathbf{A} \operatorname{vec}(\mathbf{X}) = \operatorname{vec}(\mathbf{Z}),
:::

or $\mathbf{A}\mathbf{x}=\mathbf{z}$.

The matrix $\mathbf{A}$ is $mn\times mn$; for a 12-megapixel image, it would have $1.4\times 10^{14}$ entries! It is true that the matrix is extremely sparse. But more to our point, it's unnecessary. Instead, given any vector $\mathbf{u}$ we can compute $\mathbf{v}=\mathbf{A}\mathbf{u}$ through the steps

\begin{align*}
  \mathbf{U} &= \operatorname{unvec}(\mathbf{u}),\\
  \mathbf{V} &= \operatorname{blur}(\mathbf{U}),\\
  \mathbf{v} &= \operatorname{vec}(\mathbf{V}).
\end{align*}

The following example shows how to put these ideas into practice with MINRES.

::::{prf:example} Julia demo
:class: demo
:label: demos-matrixfree-deblur
{doc}`demos/matrixfree-deblur`
::::

## Exercises

1. ✍ Show using {eq}`lintrans` and {eq}`blurfunction` that the blur operation is a linear transformation. 

2. ✍ In each case, state with reasons whether the given transformation on $n$-vectors is linear. 

    **(a)** $\mathbf{f}(\mathbf{x}) = \begin{bmatrix}
      x_2\\x_3 \\\vdots\\ x_n \\ x_1 
    \end{bmatrix},\qquad$
    **(b)** $\mathbf{f}(\mathbf{x}) = \begin{bmatrix}
      x_1\\x_1+x_2\\x_1+x_2+x_3\\\vdots\\x_1+\cdots+x_n
    \end{bmatrix}, \qquad$
    **(c)** $\mathbf{f}(\mathbf{x}) = \begin{bmatrix}
      x_1 + 1 \\x_2 + 2 \\ x_3 + 3 \\\vdots \\ x_n+n
    \end{bmatrix}$.

3. ✍ Suppose that code for the linear transformation $\mathbf{f}(\mathbf{x})=\mathbf{A}\mathbf{x}$ is given for an unknown matrix $\mathbf{A}$. Explain carefully how one could construct $\mathbf{A}$.

4. ⌨ The matrix of the blur operation happens to be symmetric and positive definite. Repeat {prf:ref}`demos-matrixfree-deblur`, adding lines to do MINRES and CG (setting the maximum number of iterations to get convergence to tolerance $10^{-5}$). Report the running time required by each of the three Krylov methods.

    ::::{only} solutions
    load mandrill
    [m,n] = size(X)
    v = [1/4 1/2 1/4];
    B = spdiags( repmat(v,m,1), -1:1, m,m);
    C = spdiags( repmat(v,n,1), -1:1, n,n);
    blur = @(X) B^12 * X * C^12;
    Z = blur(X);
    vec = @(X) reshape(X,m*n,1);
    unvec = @(x) reshape(x,m,n);
    T = @(x) vec( blur(unvec(x)) );
    tic, y = gmres(T,vec(Z),50,1e-5); time_gmres = toc
    tic, y = minres(T,vec(Z),1e-5,300); time_minres = toc 
    tic, y = pcg(T,vec(Z),1e-5,300); time_pcg = toc 
    ::::

5. The condition number of the unknown matrix associated with blurring is plausibly related to the condition numbers of the single-dimension matrices $\mathbf{B}^k$ and $\mathbf{C}^k$ in {eq}`blurfunction`.

    **(a)** ⌨  Let $m=50$. Show that $\mathbf{B}$ has a Cholesky factorization and thus is SPD. Find $\kappa(\mathbf{B})$. 

    **(b)** ✍ Explain why part (a) implies $\kappa( \mathbf{B}^k ) = \kappa(\mathbf{B})^k$.

    **(c)** ✍ Explain two important effects of the severity of blur (that is, as $k\to \infty$) on deblurring by Krylov methods. 
