# Polynomial interpolation

```{index} interpolation; by polynomials
```

In {ref}`../linsys/polyinterp` and {ref}`../localapprox/interpolation` we encountered polynomial interpolation for the $n+1$ data points $(t_0,y_0),\ldots, (t_n,y_n)$.[^nodevar] Theoretically at least, we can always construct an interpolating polynomial, and the result is unique among polynomials whose degree is less than $n+1$.

[^nodevar]: As before, we use $t_i$ to denote interpolation or data nodes, and $x$ to denote the independent variable.

::::{prf:theorem}
:label: theorem-polyinterp
If the nodes $t_0,\dots,t_n$ are all distinct, there exists a unique polynomial $p$ of degree at most $n$ that satisfies $p(t_k)=y_k$ for all $k=0,\dots,n$.
::::

::::{prf:proof}
We defer the existence part to equation {eq}`lagrangeinterp`.  As for uniqueness, if $p$ and $q$ are two interpolating polynomials, then $p-q$ is a  polynomial of degree at most $n$ that is zero at the $n+1$ points  $t_0,\dots,t_n$. By the Fundamental Theorem of Algebra, which states that a $k$th degree polynomial has no more than $k$ roots, we conclude that $p-q\equiv 0$, so $p=q$.
::::

## Lagrange formula

```{index} cardinal functions
```

In our earlier encounters with polynomial interpolation, we found the interpolant by solving a linear system of equations with a Vandermonde matrix. The first step was to express the polynomial in the natural monomial basis $1,x,x^2,\ldots$. However, as we saw in {ref}`../localapprox/pwlin`, no basis is more convenient than a cardinal basis, in which each member is one at a single node and zero at all of the other nodes. 

It is surprisingly easy to construct a cardinal basis for global polynomial interpolation. By definition, each member $\ell_{k}$ of the basis, for $k=0,\ldots,n$, is an $n$th degree polynomial satisfying the cardinality conditions

:::{math}
:label: lagrangecond
  \ell_{k}(t_j) = \begin{cases}
  1, &\text{if $j=k$,}\\
  0, & \text{otherwise.}
  \end{cases}
:::

Recall that any polynomial of degree $n$ can be expressed as

$$
c(x-r_1)(x-r_2)\dots(x-r_n) = c\prod_{k=1}^n(x-r_k),
$$

where $r_1,\dots,r_n$ are the roots of the polynomial and $c$ is a constant. The conditions {eq}`lagrangecond` give all $n$ roots of $\ell_{k}$, and the normalization $\ell_{k}(t_k)=1$ tells us how to find $c$. The result is

:::{math}
:label: lagrange
\ell_{k}(x) = \frac{(x-t_0)\dots(x-t_{k-1})(x-t_{k+1})\dots(x-t_n)}{(t_k-t_0)\dots(t_k-t_{k-1})(t_k-t_{k+1})\dots(t_k-t_n)}
= \prod_{\substack{i=0\\i\ne k}}^n \frac{(x-t_i)}{(t_k-t_i)},
:::

which is called a {term}`Lagrange polynomial`.

::::{prf:example} Julia demo
:class: demo
:label: demos-polynomial-lagrange
{doc}`demos/polynomial-lagrange`
::::

Because they are a cardinal basis, the Lagrange polynomials lead to a simple expression for the polynomial interpolating the $(t_k,y_k)$ points:

:::{margin}
Lagrange polynomials lead to a simple expression for the interpolating polynomial.
:::

:::{math}
:label: lagrangeinterp
p(x) = \sum_{k=0}^n y_k \ell_k(x).
:::

```{index} Lagrange interpolation formula
```
This is called the {term}`Lagrange formula` for the interpolating polynomial. At this point we can say that we have completed the proof of {prf:ref}`theorem-polyinterp`.

::::{prf:example}
:label: example-ClassicalLagrange
We construct the Lagrange interpolating polynomials of degrees $n=1$ and 2 to interpolate samples of $f(x) = \tan (x)$.  For $n=1$, we use $t_0= 0$ and $t_1 = \pi/3$. The Lagrange formula then gives

\begin{align*}
  P_1(x) & = y_0 \ell_0(x)  +  y_1 \ell_1(x) \\
      & = y_0 \frac{x-t_1}{t_0-t_1} + y_1 \frac{x-t_0}{t_1-t_0} \\
      & = 0 \cdot \frac{x-\frac{\pi}{3}}{0-\frac{\pi}{3}} + \sqrt{3} \cdot \frac{x-0}{\frac{\pi}{3}-0} \\
      & = \frac{3 \sqrt{3}}{\pi} x.
\end{align*}

This is the unique linear function passing through $(0,0)$ and $(\pi/3,\sqrt{3})$.

For $n=2$, we use $t_0= 0$, $_1 = \pi/6$ and $t_2 = \pi/3$. We now have

\begin{align*}
  P_2(x) & = y_0 \ell_0(x) +  y_1 \ell_1(x) +  y_2 \ell_2(x) \\
    & = y_0 \frac{(x-t_1)(x-t_2)}{(t_0-t_1)(t_0-t_2)}  +
        y_1 \frac{(x-t_0)(x-t_2)}{(t_1-t_0)(t_1-t_2)}  +
        y_2  \frac{(x-t_0)(x-t_1)}{(t_2-t_0)(t_2-t_1)} \\
    & = 0
      + \frac{1}{\sqrt{3}} \frac{\left(x-0\right)\left(x-\frac{\pi}{3}\right)}
      {\left(\frac{\pi}{6}-0\right)\left(\frac{\pi}{6}-\frac{\pi}{3}\right)}
      + \sqrt{3} \frac{\left(x-0\right)\left(x-\frac{\pi}{6}\right)}
      {\left(\frac{\pi}{3}-0\right)\left(\frac{\pi}{3}-\frac{\pi}{6}\right)}
           = \frac{6\sqrt{3}}{\pi^2}x^2 + \frac{\sqrt{3}}{\pi} x
\end{align*}

%\inputexample{globalfuncapprox}{laginterp}
%The two interpolating polynomials are shown in Figure~\ref{fig:twointerpolants}, along
%with the function $\tan x$.  Observe that $ P_1(\pi/4) = 1.2990$ and
%$P_2(\pi/4) = 1.0825$ are approximations to $\tan (\pi/4)=1$. The
%approximation improves as $n$ increases, which is typical. But when
%the interpolating polynomials are evaluated outside the interval
%$[t_0,t_n]$, they have little approximation power.
%\begin{figure}[htb]
%  \centering
%  \psfrag{n=1}[Br][B][2]{$n=1$}
%  \psfrag{tan x}[Bl][l]{$\tan x$}
%  \psfrag{x}[b][c]{$x$}
%  \psfrag{n=2}[tr][Bl][2]{$n=2$}
%  \includegraphics[width=4in]{twointerpolants}
%  \caption{The degree 1 and 2 interpolants using data from $f(x)=\tan x$ indicated by the circles.
%  }
%  \label{fig:twointerpolants}
%\end{figure}
::::


## Error formula

In addition to existence, uniqueness,  and the constructive Lagrange formula, we have a useful formula for the error in a polynomial interpolant when the data are samples of a smooth function.

::::{prf:theorem} Polynomial interpolation error
:label: theorem-interperror
Let $t_0,\dots,t_n$ be distinct points in $[a,b]$, and suppose $f$ has at least $n+1$ continuous derivatives in that interval. Let $p(x)$ be the unique polynomial of degree at most $n$ interpolating $f$ at $t_0,\dots,t_n$. Then for each $x\in[a,b]$, there exists a number $\xi(x)\in(a,b)$ such that
  
:::{math}
:label: interperror
f(x) - p(x) = \frac{f^{(n+1)}(\xi)}{(n+1)!} \prod_{i=0}^n (x-t_i).
:::
::::

::::{prf:proof}
Define
  
:::{math}
:label: lagrange-phi
\Phi(x) = \prod_{i=0}^n (x-t_i).
:::

If $x=t_i$ for some $i$, the statement of the theorem is trivially true. Otherwise, we define a new function $g(s)$ by
  
$$
g_x(s) = \Phi(s)[f(x)-p(x)] - \Phi(x)[f(s)-p(s)].
$$

Note that $x$ is now arbitrary but fixed. Clearly $g_x(t_i)=0$ for each $i=0,\dots,n$, because both $\Phi$ and the error $f-p$ have that property. Also, $g_x(x)=0$. So $g_x$ has at least $n+2$ zeros in $[a,b]$. This is possible only if $g_x$ has at least $n+1$ local minima in $(a,b)$; i.e., $g_x'$ has at least $n+1$ zeros. But that implies that $g_x''$ must have at least $n$ zeros, etc. Eventually we conclude that $g_x^{(n+1)}$ has at least one zero in $(a,b)$.[^Rolle] Let $\xi(x)$ be such a zero.

[^Rolle]: This deduction on $g_x^{(n+1)}$ is known as Rolle's Theorem in calculus.

Observe that $\Phi$ is a monic polynomial (i.e., its leading coefficient is 1) of degree $n+1$. Hence $\Phi^{(n+1)}(t)=(n+1)!$. Since $p$ has degree at most $n$, $p^{(n+1)}=0$. Finally, we write
  
\begin{align*}
0 = g_x^{(n+1)}(\xi) &= \Phi^{(n+1)}(\xi)[f(x)-p(x)] - \Phi(x)[f^{(n+1)}(\xi)-p^{(n+1)}(\xi)]\\
&= (n+1)!\,[f(x)-p(x)] - \Phi(x)f^{(n+1)}(\xi),
\end{align*}

which is a restatement of {eq}`interperror`.
::::

Usually $f^{(n+1)}$ and the function $\xi(x)$ are unknown. The importance of the formula {eq}`interperror` is how it helps to express the error as a function of $x$, and its dependence on the nodes $t_0,\dots,t_n$. We will exploit this knowledge later.

::::{prf:example} Julia demo
:class: demo
:label: demos-polynomial-error
{doc}`demos/polynomial-error`
::::

For equispaced nodes, {prf:ref}`theorem-polyinterp` has an immediate consequence. 

::::{prf:corollary}
:label: theorem-polyequi
Suppose $t_i=i h$ for constant step size $h$ and all $i=0,1,\ldots,n$, and that $f$ has $n+1$ continuous derivatives in $(t_0,t_n)$. If $x\in[t_0,t_n]$, then there exists $\xi(x)\in(t_0,t_n)$ and $C$ independent of $x$ such that
  
:::{math}
:label: equisperror
|f(x) - p(x)| \le C f^{(n+1)}(\xi) h^{n+1}.
:::

In particular, $|f(x)-p(x)|=O(h^{n+1})$ as $h\to 0$.
::::

::::{prf:proof}
If $x\in[t_0,t_n]$, then $|x-t_i|<nh$ for all $i$, and {eq}`interperror` implies {eq}`equisperror`. As $h\to 0$, $\xi\to x$, and the continuity of $f^{(n+1)}$ allows us to make the asymptotic conclusion.
::::

```{index} order of accuracy
```

Thus, linear (and piecewise linear) interpolation on an interval of width $O(h)$ is $O(h^2)$, and a finite difference method based on $n+1$ nodes is $O(h^n)$ (because of division by $h$ in the finite difference formula).

## Instability

As presented in {eq}`lagrangeinterp`, the Lagrange formula is not a good choice for numerical computation, because it is unstable (see [this exercise](problem-lagrange-instability).) In the next section we derive an algebraically equivalent formula that is numerically stable and faster to apply as well.

## Exercises

(problem-lagrange)=
1. ✍ Write out the Lagrange form of the interpolating polynomial of degree $n$ for the given functions and nodes. Using a calculator, evaluate the polynomial at $x=\pi/4$ and compute the error there.

    **(a)** $f(x) = \sin(x), \ n=1, \ t_0=0, t_1 = \pi/2$

    **(b)** $f(x) = \sin(x), \ n=2, \ t_0=0, t_1 = \pi/6, t_2 = \pi/2$

    **(c)** $f(x) = \cos(x), \ n=2, \ t_0=0, t_1 = \pi/3, t_2 = \pi/2$

    **(d)** $f(x) = \tanh(x), \ n=2, \ t_0=0, t_1 = \pi/3, t_2 = \pi/2$

    ::::{only} solutions
    %%
    %% (a)
    %% (b)
    % The points are $(0,0)$, $(\pi/6,1/2)$, and $(\pi/2,1)$.
    %
    % $$ 0 + \frac{x(x-\pi/2)}{(\pi/6)(-\pi/3)}\cdot \frac{1}{2} + \frac{x(x-\pi/6)}{(\pi/2)(2\pi/3)} $$
    %%
    p = polyinterp([0 pi/6 pi/2]',[0 1/2 1]');
    p(pi/4)
    ::::

2. ⌨ For each case, plot the requested Lagrange cardinal polynomial for the given set of nodes over the interval $[t_0,t_n]$. Superimpose dots or circles for the points represented by the cardinal conditions {eq}`lagrangecond`. 

    **(a)** $n=2,\quad t_0=-1, \, t_1=-0.2,\, t_2=0, \quad \ell_2(x)$
    
    **(b)** $n=4,\quad t_0=0, \, t_1=1,\, t_2=1.5,\, t_3=2.5,\, t_4=3, \quad \ell_3(x)$

    **(c)** $n=20, \quad t_i=i/n \text{ for } i=0,\ldots,n, \quad \ell_0(x)$

    **(d)** $n=20, \quad t_i=i/n \text{ for } i=0,\ldots,n, \quad \ell_{10}(x)$ 

    **(e)** $n=40, \quad t_i=i/n \text{ for } i=0,\ldots,n, \quad \ell_{20}(x)$ 

    ::::{only} solutions
    ::::

3. ✍ Suppose $p$ is the quadratic polynomial interpolating the points $(-2,12)$, $(1,3a)$, and $(2,0)$. Use {eq}`lagrangeinterp` to compute $p'(0)$. 

4. ✍ Explain carefully why using {eq}`lagrangeinterp` to compute $p(x)$ at a single value of $x$ takes $O(n^2)$ floating point operations.

    ::::{only} solutions
    ::::

5. ✍  Explain why for any distribution of nodes and all $x$,
    
    $$
    1 = \sum_{k=0}^n \ell_k(x).
    $$
   
    (Hint: This problem does not require any computation or formula manipulation.) 

6. ✍ Show that
    
    $$
    \ell_k(x) = \frac{\Phi(x)}{(x-t_k)\Phi'(t_k)},
    $$

    where $\Phi$ is the function defined in {eq}`lagrange-phi`. 

    (problem-lagrange-instability)=
7. ✍ Consider the nodes $t_0=0$, $t_1=1$, $t_2=\beta$, where $\beta>1$.

    **(a)** Write out the Lagrange cardinal polynomials $\ell_0$, $\ell_1$, and $\ell_2$.      
    
    **(b)** Suppose the data are $y_0=y_1=y_2=1$. What is the unique interpolating polynomial of degree no greater than 2? 

    **(c)** By letting $x=1/2$ and $\beta\to 1$, explain how parts (a) and (b) demonstrate a numerical instability in the Lagrange formula. 
