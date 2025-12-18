---
numbering:
  enumerator: 8.7.%s
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.12
---
```{code-cell}
:tags: [remove-cell]
import Pkg; Pkg.activate("/Users/driscoll/Documents/GitHub/fnc")

using FNCFunctions

using Plots
default(
    titlefont=(11,"Helvetica"),
    guidefont=(11,"Helvetica"),
    linewidth = 2,
    markersize = 3,
    msa = 0,
    size=(500,320),
    label="",
    html_output_format = "svg"
)

using PrettyTables, LaTeXStrings, Printf
using LinearAlgebra

@ptconf backend = Val(:html) tf = tf_html_simple
```

(section-krylov-matrixfree)=

# Matrix-free iterations

In Chapter 4, we solved the nonlinear rootfinding problem $\mathbf{f}(\mathbf{x})=\boldsymbol{0}$ with methods that needed only the ability to evaluate $\mathbf{f}$ at any known value of $\mathbf{x}$. By repeatedly evaluating $\mathbf{f}$ at cleverly chosen points, these algorithms were able to return an estimate for $\mathbf{f}^{-1}(\boldsymbol{0})$.

We can explore the same idea in the context of linear algebra by shifting our viewpoint from matrices to *linear transformations*. If we define $\mathbf{f}(\mathbf{x})=\mathbf{A}\mathbf{x}$, then for all vectors $\mathbf{x}$, $\mathbf{y}$, and scalars $\alpha$,

```{math}
:label: lintrans
\begin{split}
\mathbf{f}(\mathbf{x} + \mathbf{y} ) &= \mathbf{f}(\mathbf{x}) + \mathbf{f}(\mathbf{y} ), \\
\mathbf{f}(\alpha \mathbf{x} ) & = \alpha\, \mathbf{f}(\mathbf{x}).
\end{split}
```

These properties define a linear transformation. Moreover, *every* linear transformation between finite-dimensional vector spaces can be represented as a matrix-vector multiplication.

A close examination reveals that in the power iteration and Krylov subspace methods, the only appearance of the matrix $\mathbf{A}$ is to apply it a known vector, i.e., to evaluate the linear transformation $\mathbf{f}(\mathbf{x})=\mathbf{A}\mathbf{x}$. If we have access to $\mathbf{f}$, we don't need the matrix at all! That is, Krylov subspace methods can be used to invert a linear transformation if one provides code for the transformation, even if its associated matrix is not known explicitly. That may sound like a strange situation, but it is not uncommon.

## Blurring images

```{index} image (as a matrix)
```

In {numref}`section-matrixanaly-insight` we saw that a grayscale image can be represented as an $m\times n$ matrix $\mathbf{X}$ of pixel intensity values. Now consider a simple model for blurring the image. Define $\mathbf{B}$ as the $m\times m$ tridiagonal matrix

```{math}
:label: blurmatrix
B_{ij} =
\begin{cases}
\tfrac{1}{2} & \text{if $i=j$},\\
\tfrac{1}{4} & \text{if $|i-j|=1$},\\
0 & \text{otherwise.}
\end{cases}
```

The product $\mathbf{B}\mathbf{X}$ applies $\mathbf{B}$ to each column of $\mathbf{X}$. Within that column it does a weighted average of the values of each pixel and its two neighbors. That has the effect of blurring the image vertically. We can increase the amount of blur by applying $\mathbf{B}$ repeatedly.

In order to blur horizontally, we can transpose the image and apply blurring in the same way. We need a blurring matrix defined as in {eq}`blurmatrix` but with size $n\times n$. We call this matrix $\mathbf{C}$. Altogether the horizontal blurring is done by transposing, applying $\mathbf{C}$, and transposing back to the original orientation. That is,

```{math}
\bigl(\mathbf{C} \mathbf{X}^T\bigr)^T = \mathbf{X}\mathbf{C}^T = \mathbf{X}\mathbf{C},
```

using the symmetry of $\mathbf{C}$. So we can describe blur in both directions as the function

```{math}
:label: blurfunction
\operatorname{blur}(\mathbf{X}) = \mathbf{B}^k \mathbf{X} \mathbf{C}^k
```

for a positive integer $k$.

::::{prf:example} Blurring an image
:label: demo-matrixfree-blur

We use a readily available test image.

```{code-cell}
using Images, TestImages
img = testimage("mandrill")
m, n = size(img)
X = @. Float64(Gray(img))
plot(Gray.(X), title="Original image", aspect_ratio=1)
```

We define the one-dimensional tridiagonal blurring matrices.

```{code-cell}
using SparseArrays
function blurmatrix(d)
    v1 = fill(0.25, d - 1)
    return spdiagm(0 => fill(0.5, d), 1 => v1, -1 => v1)
end
B, C = blurmatrix(m), blurmatrix(n);
```

Finally, we show the results of using $k=12$ repetitions of the blur in each direction.

```{code-cell}
using Plots
blur = X -> B^12 * X * C^12;
Z = blur(X)
plot(Gray.(Z), title="Blurred image")
```

::::

## Deblurring

A more interesting operation is *deblurring*: given an image blurred by poor focus, can we reconstruct the true image? Conceptually, we want to invert the function $\operatorname{blur}(\mathbf{X})$.

It's easy to see from {eq}`blurfunction` that the blur operation is a linear transformation on image matrices. But an $m\times n$ image matrix is equivalent to a length-$mn$ vector—it's just a matter of interpreting the shape of the same data. Let $\operatorname{vec}(\mathbf{X})=\mathbf{x}$ and $\operatorname{unvec}(\mathbf{x})=\mathbf{X}$ be the mathematical statements of such reshaping operations. Now say $\mathbf{X}$ is the original image and $\mathbf{Z}=\operatorname{blur}(\mathbf{X})$ is the blurred one. Then by linearity there is some matrix $\mathbf{A}$ such that

```{math}
\mathbf{A} \operatorname{vec}(\mathbf{X}) = \operatorname{vec}(\mathbf{Z}),
```

or $\mathbf{A}\mathbf{x}=\mathbf{z}$.

The matrix $\mathbf{A}$ is $mn\times mn$; for a 12-megapixel image, it would have $1.4\times 10^{14}$ entries! Admittedly, it is extremely sparse, but the point is that we don't need it at all.

Instead, given any vector $\mathbf{u}$ we can compute $\mathbf{v}=\mathbf{A}\mathbf{u}$ through the steps

\begin{align*}
  \mathbf{U} &= \operatorname{unvec}(\mathbf{u}), \\
  \mathbf{V} &= \operatorname{blur}(\mathbf{U}), \\
  \mathbf{v} &= \operatorname{vec}(\mathbf{V}).
\end{align*}

The following example shows how to put these ideas into practice with MINRES.

::::{prf:example} Deblurring an image
:label: demo-matrixfree-deblur

We repeat the earlier process to blur an original image $\mathbf{X}$ to get $\mathbf{Z}$.

```{code-cell}
:tags: hide-input

img = testimage("lighthouse")
m, n = size(img)
X = @. Float64(Gray(img))

B = spdiagm(0 => fill(0.5, m),
    1 => fill(0.25, m - 1), -1 => fill(0.25, m - 1))
C = spdiagm(0 => fill(0.5, n),
    1 => fill(0.25, n - 1), -1 => fill(0.25, n - 1))
blur = X -> B^12 * X * C^12
Z = blur(X)
plot(Gray.(Z), title="Blurred image")
```

Now we imagine that $\mathbf{X}$ is unknown and that we want to recover it from $\mathbf{Z}$. We first need functions that translate between vector and matrix representations.

```{code-cell}
# vec (built-in) converts matrix to vector
unvec = z -> reshape(z, m, n);  # convert vector to matrix
```

```{index} ! Julia; LinearMap
```

Now we declare the three-step blur transformation as a `LinearMap`, supplying also the size of the vector form of an image.

```{code-cell}
using LinearMaps
T = LinearMap(x -> vec(blur(unvec(x))), m * n);
```

The blurring operators are symmetric, so we apply `minres` to the composite blurring transformation `T`.
```{tip}
:class: dropdown
The function `clamp01` in `Images` restricts values to be in the interval $[0,1]$.
```

```{code-cell}
using IterativeSolvers
y = minres(T, vec(Z), maxiter=50, reltol=1e-5);
Y = unvec(clamp01.(y))

plot(Gray.(X), layout=2, title="Original")
plot!(Gray.(Y), subplot=2, title="Deblurred")
```

::::

## Exercises

``````{exercise}
:label: problem-matrixfree-blur
✍ Show using {eq}`lintrans` and {eq}`blurfunction` that the blur operation is a linear transformation. 
``````

``````{exercise}
:label: problem-matrixfree-linearity
✍ In each case, state with reasons whether the given transformation on $n$-vectors is linear. 

**(a)** $\,\mathbf{f}(\mathbf{x}) = \begin{bmatrix} x_2\\x_3 \\\vdots\\ x_n \\ x_1 \end{bmatrix}\qquad$
**(b)** $\,\mathbf{f}(\mathbf{x}) = \begin{bmatrix} x_1\\x_1+x_2\\x_1+x_2+x_3\\\vdots\\x_1+\cdots+x_n \end{bmatrix} \qquad$
**(c)** $\,\mathbf{f}(\mathbf{x}) = \begin{bmatrix} x_1 + 1 \\x_2 + 2 \\ x_3 + 3 \\\vdots \\ x_n+n \end{bmatrix} \qquad$
**(d)** $\,\mathbf{f}(\mathbf{x}) = \|\mathbf{x}\|_\infty\, \mathbf{e}_1$
``````

``````{exercise}
:label: problem-matrixfree-reconstruct
✍ Suppose that code for the linear transformation $\mathbf{f}(\mathbf{x})=\mathbf{A}\mathbf{x}$ is given for an unknown matrix $\mathbf{A}$. Explain carefully how one could construct $\mathbf{A}$.
``````

``````{exercise}
:label: problem-matrixfree-cg
⌨ The matrix of the blur transformation happens to be symmetric and positive definite. Repeat @demo-matrixfree-deblur using CG for the deblurring.
``````

``````{exercise}
:label: problem-matrixfree-cond
The condition number of the matrix of the blur transformation is related to the condition numbers of the single-dimension matrices $\mathbf{B}^k$ and $\mathbf{C}^k$ in {eq}`blurfunction`.

**(a)** ⌨  Let $m=50$. Show that $\mathbf{B}$ has a Cholesky factorization and thus is SPD. Find $\kappa(\mathbf{B})$. (Note: `cond` requires a regular dense matrix, not a sparse matrix.)

**(b)** ✍ Explain why part (a) implies $\kappa( \mathbf{B}^k ) = \kappa(\mathbf{B})^k$.

**(c)** ✍ Explain two important effects of the limit $k\to \infty$ on deblurring by Krylov methods. 
``````

``````{exercise}
:label: problem-matrixfree-cumsum
The cumulative summation function `cumsum` is defined as

$$
\mathbf{f}(\mathbf{x}) = \begin{bmatrix} x_1 \\ x_1+x_2 \\ \vdots \\ x_1 + x_2 + \cdots + x_n \end{bmatrix}.
$$

**(a)** ✍ Show that $\mathbf{f}$ is a linear transformation.

**(b)** ⌨ Define vector $\mathbf{b}$ by $b_i = (i/100)^2$ for $i=1,\ldots,100$. Then use `gmres` to find $\mathbf{x}=\mathbf{f}^{-1}(\mathbf{b})$. 

**(c)** ⌨ Plot $\mathbf{x}$, and explain why the result looks as it does.
``````
