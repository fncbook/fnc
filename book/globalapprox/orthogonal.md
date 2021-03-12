# Orthogonal polynomials

Interpolation is not the only way to use polynomials for global approximation of functions. We have also seen how to find least-squares polynomial fits to data, by solving linear least-squares matrix problems. This idea can be extended from fitting data directly to fitting functions.

::::{prf:example} Julia demo
:class: demo
:label: demos-orthogonal-approx
{doc}`demos/orthogonal-approx`
::::

```{index} inner product; of functions
```
We can extend least-squares fitting from data to functions by extending several familiar finite-dimensional definitions to functions. Let $S$ be the set of continuous real-valued functions on the interval $[-1,1]$. By analogy with the inner product of two vectors, $\mathbf{u}^T\mathbf{v} = \sum u_iv_i$, we can define the {term}`inner product` of any functions $f$ and $g$ in $S$:

:::{math}
:label: leginner
\langle f,g \rangle  = \int_{-1}^1 f(x)g(x)\,dx.
:::

The inner product assigns a real number to any pair of functions in $S$. In abstract usage we refer to $S$ as an **inner product space**. This new definition of inner product allows us to extend to $S$ familiar vector concepts such as the 2-norm,

:::{math}
:label: fun2norm
\|f\|_2^2 = \langle f,f \rangle,
:::

and orthogonality,

$$
\langle f,g \rangle = 0.
$$


## Quasimatrices

If we are extending our notion of vectors to include continuous functions, what should serve as an extension of a matrix? One of our most important interpretations of a matrix is as a collection of its columns. Consider, for instance, the Vandermonde system $\mathbf{V} \mathbf{c} = \mathbf{y}$ from {eq}`vandersystem`. We can connect it to the polynomial interpolation problem by writing the product as a linear combination of columns:

$$
  \mathbf{y} = \mathbf{V} \mathbf{c} = c_0 \mathbf{v}_0 + c_1 \mathbf{v}_1 + \cdots c_n \mathbf{v}_n,
$$

in which the $i$th row of $\mathbf{v}_j$ is $t_i^j$. We can interpret the columns $\mathbf{v}_0,\ldots,\mathbf{v}_n$ as discrete versions of the monomial functions $x^j$ for $j=0,\ldots,n$. The identical linear combination expansion holds in the case of a linear least-squares fit, except that there are more rows (evaluation nodes) than columns (monomial terms), and the equality is replaced by approximation.

Now consider the linear combination

:::{math}
  :label: quasimatvec
  c_0\,\underline{1} + c_1\, \underline{x} + \cdots + c_n \,\underline{x^n},
:::

```{index} quasimatrix
```
where the "column vectors" are actually the monomial functions in the function space $S$. We underline the monomials here just to emphasize the connection to abstract vector spaces. We are motivated to define the {term}`quasimatrix`

:::{math}
  :label: quasimat
  \mathbf{V} =
  \begin{bmatrix}
     \underline{1} & \underline{x} & \cdots &  \underline{x^n}
  \end{bmatrix}.
:::

The quasimatrix-vector product $\mathbf{V}\mathbf{c}$ is to be interpreted as {eq}`quasimatvec`. Moreover, given a function $g(x)$, we define

:::{math}
  :label: quasimatfun
  \mathbf{V}^T g =
  \begin{bmatrix}
    \langle 1,g \rangle \\ \vdots  \\  \langle x^n,g \rangle
  \end{bmatrix}.
:::

Finally, we define the operation $\mathbf{W}=\mathbf{V}^T\mathbf{V}$  as the $(n+1)\times (n+1)$ matrix having $W_{ij} = \langle {x^i},{x^j}\rangle$.

We are not limited to monomial functions as the columns of a quasimatrix. If $f_1(x),\ldots,f_k(x)$ are any $k$ functions in the space $S$, then we can define

$$
\mathbf{F}(x) = \begin{bmatrix}
     \underline{f_1} & \underline{f_2} & \cdots &  \underline{f_k}
  \end{bmatrix}.
$$

Then $\mathbf{F}\mathbf{c}$ is a linear combination of the functions, $\mathbf{F}^T{g}$ is a vector of inner products, and $\mathbf{F}^T\mathbf{F}$ is a matrix of inner products.  We consider any other expressions involving $\mathbf{F}$ to be undefined. It might help to think of $\mathbf{F}$ as an "$\infty\times k$" matrix, which is consistent with the definitions that $\mathbf{F}\mathbf{c}$ is a function ($\infty\times 1$), $\mathbf{F}^T g$ is a vector ($k\times 1$), and $\mathbf{F}^T\mathbf{F}$ is a matrix ($k \times k$).

## Normal equations

Let us return to the general discrete linear least squares problem (see {eq}`linls`) of minimizing  $\| \mathbf{f} - \mathbf{V} \mathbf{c} \|_2$ over all possible $\mathbf{c}$, given matrix $\mathbf{V}$ and vector $\mathbf{f}$. Its solution is given by the normal equations {eq}`normaleqns`,

:::{math}
:label: contLSsol
\mathbf{c} = \bigl[\mathbf{V}^T\mathbf{V}\bigr]^{-1} \mathbf{V}^T \mathbf{f}.
:::

Now suppose $\mathbf{V}$ is the quasimatrix {eq}`quasimat` and $f\in S$. Thanks to our extended definitions, we can interpret $\mathbf{V}\mathbf{c}$ as a polynomial, and $\| f-\mathbf{V}\mathbf{c} \|_2$ has a precise meaning. Minimizing this norm gives the best polynomial approximation to $f$ in the least squares sense. Remarkably, equation {eq}`contLSsol` still provides the solution in this context! The underlying reason is that with our new definitions, $S$ captures the only ingredients necessary to prove the optimality of the normal equations: linear combination and inner product.

::::{prf:example}
:label: example-lsfitexpfun
We revisit approximation of $e^x$ as suggested in {prf:ref}`demos-orthogonal-approx`. With $\mathbf{V}= \begin{bmatrix} \underline{1} & \underline{x} \end{bmatrix}$, we get
  
  $$
    \mathbf{V}^Te^x =
    \begin{bmatrix}
      \langle 1,e^x \rangle \\[1mm] \langle x,e^x \rangle
    \end{bmatrix} =
    \begin{bmatrix}
      \int_{-1}^1 e^x\, dx \\[1mm] \int_{-1}^1 x e^x\, dx
    \end{bmatrix} =
    \begin{bmatrix}
      e-e^{-1} \\ 2 e^{-1}
    \end{bmatrix},
  $$

  and

  $$
    \mathbf{V}^T \mathbf{V} =
    \begin{bmatrix}
      \langle 1,1 \rangle & \langle 1,x \rangle \\[1mm] \langle x,1 \rangle & \langle x,x \rangle
    \end{bmatrix} =
    \begin{bmatrix}
      2 & 0 \\ 0 & 2/3
    \end{bmatrix}.
  $$

  The normal equations {eq}`contLSsol` therefore have solution
  
  $$
    \mathbf{c} = \begin{bmatrix}
      2 & 0 \\ 0 & 2/3
    \end{bmatrix}^{-1}
    \begin{bmatrix}
      e-e^{-1} \\ 2 e^{-1}
    \end{bmatrix}
    =
    \begin{bmatrix}
      \sinh(1) \\ 3e^{-1}
    \end{bmatrix} \approx
    \begin{bmatrix}
      1.175201\\ 1.103638
    \end{bmatrix},
  $$

  which is well in line with the values found in {prf:ref}`demos-orthogonal-approx`.
::::

## Legendre polynomials

```{index} orthogonal; polynomials
```
Equation {eq}`contLSsol` becomes much simpler if $\mathbf{V}^T\mathbf{V}$ is diagonal. By our definitions, this would imply that the columns of $\mathbf{V}$ are mutually orthogonal in the sense of the function inner product. Such is not true of the monomial functions $x^j$. But there are {term}`orthogonal polynomials` which do satisfy this property.

Let $\mathcal{P}_n \subset S$ be the set of polynomials of degree $n$ or less. Define a sequence of polynomials by

:::{math}
 :label: legendre
 \begin{split}
     P_0(x) &= 1, \\
     P_1(x) &= x, \\
     P_{k}(x) &= \frac{2k-1}{k}xP_{k-1}(x) - \frac{k-1}{k}P_{k-2}(x), \qquad k = 2,3,\ldots.
 \end{split}
:::

```{index} Legendre polynomials
```
These are known as the **Legendre polynomials**. One can show that $P_0,\ldots,P_n$ form a basis for $\mathcal{P}_n$. Most significantly, the Legendre polynomials are orthogonal, because

```{margin}
The Legendre polynomials are orthogonal in the continuous inner product over $[-1,1]$.
```

:::{math}
 :label: legendreortho
 \langle P_i,P_j \rangle =
 \begin{cases}
   0, & i \neq j \\
   \alpha_i^2 = \bigl(i+\tfrac{1}{2}\bigr)^{-1}, & i=j.
 \end{cases}
:::

If we define the quasimatrix

:::{math}
  :label: legquasi
  \mathbf{L}_n(x) =
  \begin{bmatrix}
    \alpha_0^{-1} \underline{P_0} & \alpha_1^{-1} \underline{P_1} & \cdots
    & \alpha_{n}^{-1} \underline{P_{n}}
  \end{bmatrix},
:::

then $\mathbf{L}_n^T\mathbf{L}_n=\mathbf{I}$. The normal equations {eq}`contLSsol` thus simplify accordingly. Unraveling the definitions, we find the least-squares solution

:::{math}
  :label: funlssoln
  \mathbf{L}_n \bigl( \mathbf{L}_n^T f \bigr) = \sum_{k=0}^n c_k P_k(x), \quad \text{where }
  c_k = \frac{1}{\alpha_k^2} \langle P_k,f \rangle.
:::

## Roots of orthogonal polynomials

Interesting properties can be deduced from the orthogonality conditions. The following result will be relevant in an [upcoming section](integration.md).

::::{prf:theorem}
:label: theorem-legroot
All $n$ roots of the Legendre polynomial $P_n(x)$ are simple and real, and they lie in the open interval $(-1,1)$.
::::


::::{prf:proof}
  Let $x_1,\ldots,x_m$ be all of the distinct roots of $P_n(x)$ between $-1$ and $1$ at which $P_n(x)$ changes sign (in other words, all roots of odd multiplicity). Define
  
$$
    r(x) = \prod_{i=1}^m (x-x_i).
  $$

  By definition, $r(x)P_n(x)$ does not change sign over $(-1,1)$. Therefore
  
:::{math}
    :label: legrootint
    \int_{-1}^1 r(x)P_n(x) \, dx \neq 0.
  :::

  Because $r$ is a degree-$m$ polynomial, we can express it as a combination of $P_0,\ldots,P_m$. If $m<n$, the integral {eq}`legrootint` would be zero, by the orthogonality property of Legendre polynomials. So $m\ge n$. Since $P_n(x)$ has at most $n$ real roots, $m=n$. All of the roots must therefore be simple, and this completes the proof.
::::

## Chebyshev polynomials

Equation {eq}`leginner` is not the only reasonable way to define an inner product on a function space. It can be generalized to

:::{math}
  :label: weightedinner
  \langle f,g \rangle = \int_{-1}^1f(x)g(x)w(x)\,dx,
:::

for a positive function $w(x)$ called the **weight function** of the inner product. An important special case is

:::{math}
  :label: chebinner
  \langle f,g \rangle  = \int_{-1}^1 \frac{f(x)g(x)}{\sqrt{1-x^2}}\,dx.
:::

```{index} Chebyshev polynomials
```
The polynomials that are orthogonal with respect to this weight function are the 
 **Chebyshev polynomials**, which have their own recursive definition:

 ```{margin}
 The Chebyshev polynomials are orthogonal with respect to a particular weight function over $[-1,1]%.
 ```

:::{math}
 :label: chebyshev
 \begin{split}
     T_0(x) &= 1, \\
     T_1(x) &= x, \\
     T_{k}(x) &= 2xT_{k-1}(x) - T_{k-2}(x) \qquad k = 2,3,\ldots.
 \end{split}
:::

As in the Legendre case, each $T_k$ is a polynomial of degree exactly $k$. Chebyshev polynomials also have a surprising  alternate form,

:::{math}
:label: chebaltform
T_k(x) = \cos\left( k \arccos(x) \right).
:::


In the weighted inner product, $\langle T_i,T_j \rangle=0$ for $i\neq j$, and $\langle T_i,T_i \rangle=\gamma_i^2$, where $\gamma_0^2=\pi$ and $\gamma_i^2=\pi/2$ for $i>0$. The rest of the definitions relating to least-squares problems remain the same under the new inner product, if we replace Legendre polynomials with Chebyshev polynomials. Note that the least-squares solution is not the same in the Legendre and Chebyshev cases; both find the closest approximation to a given $f(x)$, but the norm used to measure distances is not the same.

{prf:ref}`theorem-legroot` also applies to Chebyshev polynomials. In fact, thanks to {eq}`chebaltform`, the roots of $T_n$ are known explicitly:

:::{math}
  :label: chebfirstpoints
  t_k = \cos\left(\frac{2k-1}{2n}\pi\right), \qquad k=1,\ldots,n.
:::

```{index} Chebyshev points!first kind
```
These are known as the **Chebyshev points of the first kind**. The chief difference between first-kind and second-kind points is that the latter type include the endpoints $\pm 1$. Both work well for polynomial interpolation and give spectral convergence.

## Least squares versus interpolation

Both interpolation and the solution of a linear least-squares problem produce a **projection** of the original function $f$ onto the space of polynomials $\mathcal{P}_n$. A projection is a linear operator or matrix satisfying the identity $\m{P}^2=\m{P}$. In words, applying a projection a second time has no effect.

```{index} orthogonal projection
```
The least-squares case is special in that it produces an *orthogonal projection*. Minimization of the residual over the range of a matrix makes it orthogonal to the range, as shown in {prf:ref}`theorem-normaleqns`. The close connection with inner products and orthogonality makes the 2-norm (the norm established by {eq}`fun2norm`, perhaps with a weight function) a natural setting for analysis. Because the unit weight function is the simplest choice for the inner product, Legendre polynomials are commonly used for least squares.

Interpolation has no easy connection to inner products or the 2-norm. With interpolants it's more fruitful to perform a different kind of approximation analysis, often involving the complex plane, in which the infinity-norm (max-norm) is the natural choice. For reasons beyond the scope of this text, Chebyshev polynomials are typically the most convenient to work with in this context.

There are many other families of orthogonal polynomials defined by the inner product {eq}`weightedinner` for different weight functions that may be favorable in particular problems. While Legendre and Chebyshev polynomials are the most prominent cases, any orthogonal family can be used with either least squares or interpolation.

## Exercises

1. ✍ Let $\mathbf{F}$ be the quasimatrix $\Bigl[ \underline{1} \quad \underline{\cos(\pi x)}\quad \underline{\sin(\pi x)}\Bigr]$ (all on $x\in[-1,1]$). Find $\mathbf{F}^T x$ and $\mathbf{F}^T \mathbf{F}$.

2. ✍ Find the best linear approximation (in the least squares sense) to the function $\sin(x)$ on $[-1,1]$. 

3. 
    **(a)** ✍ ⌨ Use {eq}`legendre` to write out $P_2(x)$ and $P_3(x)$. Plot $P_0,P_1,P_2,P_3$ on one graph for $-1\le x \le 1$. 
    
    **(b)** ✍ ⌨  Use {eq}`chebyshev` to write out $T_2(x)$ and $T_3(x)$. Plot $T_0,T_1,T_2,T_3$ on one graph for $-1\le x \le 1$.

    ::::{only} solutions
    %% (a)
    % $$P_2(x) = \frac{3}{2}x^2 - \frac{1}{2},\quad P_3(x) = \frac{5}{2}x^3 - \frac{3}{2}x$$
    fplot(@(x) 1,[-1 1])
    hold on
    fplot(@(x) x,[-1 1])
    fplot(@(x) (3*x^2-1)/2,[-1 1])
    fplot(@(x) (5*x^3-3*x)/2,[-1 1])
    %% (b)
    % $$T_2(x) = 2x^2 - 1,\quad T_3(x) = 4x^3 - 3x$$
    clf
    fplot(@(x) 1,[-1 1])
    hold on
    fplot(@(x) x,[-1 1])
    fplot(@(x) 2*x^2-1,[-1 1])
    fplot(@(x) 4*x^3-3*x,[-1 1])
    ::::

4. ✍ Use {eq}`legendre` to show that $P_n(x)$ is an odd function if $n$ is odd and an even function if $n$ is even.

5. ⌨ Using {eq}`legendre`, write a function `legpoly(x,n)` that returns a matrix whose columns are the Legendre polynomials $P_0,P_1,\ldots,P_n$ evaluated at all the points in the vector $x$. Then use your function to plot $P_0,P_1,P_2,P_3$ on one graph. 

6. ⌨ (Continuation of previous problem.) Choose 1600 evenly spaced points in $[-1,1]$. For $n=1,2,\ldots,16$, use this vector of points and the function `legpoly` to construct a $1600\times (n+1)$ matrix that discretizes the quasimatrix
   
    $$
     \mathbf{A}_n = 
    \begin{bmatrix}
     \underline{P_0} &  \underline{P_1} & \cdots 
    &  \underline{P_{n}}
    \end{bmatrix}.
    $$

    Make a table of the matrix condition number $\kappa(\mathbf{A}_n)$ as a function of $n$. (These will not be much larger than one, thanks to the connection to the quasimatrix $\mathbf{L}_n$ in {eq}`legquasi`.) 

7. ⌨ Using {eq}`chebaltform`, write a function `chebpoly` that returns a matrix whose columns are the Chebyshev polynomials $T_0,T_1,\ldots,T_n$ evaluated at all the points in the vector $x$. Then use your function to plot $T_0,T_1,T_2,T_3$ on one graph. 

    ::::{only} solutions
    %%
    n = 16;
    x = linspace(-1,1,1600)';
    L = [ x.^0, x ];
    for j = 2:n
    L(:,j+1) = ( (2*j-1)*x.*L(:,j) - (j-1)*L(:,j-1) ) / (j);
    end
    ::::

8. 
    **(a)** ✍ Use {eq}`chebaltform` to show that the first-kind points {eq}`chebfirstpoints` are roots of $T_n$. 
    
    **(b)** ✍ Use {eq}`chebaltform` to show that the second-kind points {eq}`chebextreme` are local extreme points of $T_n$.

    ::::{only} solutions
    %%
    %% (a)
    % First,
    %% 
    % $$\arccos(t_k) = \frac{2k-1}{2n}\pi.$$ 
    %%
    % Thus for any integer $k$,
    %%
    % $$T_n(x_k) = \cos\left( n \frac{2k-1}{2n}\pi \right) = \cos\left[ \left(k-\frac{1}{2}\right)\pi \right]=0$$
    %% (b)
    % Now 
    %%
    % $$\arccos(t_k) = \frac{k}{n}\pi.$$ 
    %%
    % Thus for any integer $k$,
    %%
    % $$T_n(x_k) = \cos\left( n \frac{k}{n}\pi \right) = \cos(k\pi) = (-1)^k,$$
    %%
    % which is always an extreme point.
    ::::


9. ✍ Show that the definition {eq}`chebaltform` satisfies the recursion relation in {eq}`chebyshev`.

10. ✍ Use {eq}`chebaltform` to show that $\langle T_0,T_0 \rangle=\pi$ and $\langle T_k,T_k \rangle=\pi/2$ for $k>0$ in the Chebyshev-weighted inner product. 

