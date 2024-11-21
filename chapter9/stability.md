---
jupytext:
  cell_metadata_filter: -all
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.10.3
kernelspec:
  display_name: Julia 1.7.1
  language: julia
  name: julia-fast
---
```{code-cell} 
:tags: [remove-cell]
using FundamentalsNumericalComputation
FNC.init_format()
```

(section-globalapprox-stability)=
# Stability of polynomial interpolation

```{index} stability; of polynomial interpolation
```

With  barycentric interpolation available in the form of {numref}`Function {number} <function-polyinterp>`, we can explore polynomial interpolation using a numerically stable algorithm. Any remaining sensitivity to error is due to the conditioning of the interpolation process itself.

(demo-stability-equispaced)=
```{proof:demo}
```

```{raw} html
<div class='demo'>
```

```{raw} latex
%%start demo%%
```

We choose a function over the interval $[0,1]$. 

```{code-cell} 
f = x -> sin(exp(2*x));
```

Here is a graph of $f$ and its polynomial interpolant using seven equally spaced nodes.

```{code-cell} 
:tags: [hide-input]
plot(f,0,1,label="function",legend=:bottomleft)
t = (0:6)/6
y = f.(t)
scatter!(t,y,label="nodes")

p = FNC.polyinterp(t,y)
plot!(p,0,1,label="interpolant",title="Equispaced interpolant, n=6")
```

This looks pretty good. We want to track the behavior of the error as $n$ increases. We will estimate the error in the continuous interpolant by sampling it at a large number of points and taking the max-norm.

```{code-cell} 
:tags: [hide-input]
n = 5:5:60;  err = zeros(size(n))
x = range(0,1,length=2001)      # for measuring error
for (i,n) in enumerate(n) 
  t = (0:n)/n                   # equally spaced nodes
  y = f.(t)                     # interpolation data
  p = FNC.polyinterp(t,y)
  err[i] = norm( (@. f(x)-p(x)), Inf )
end
plot(n,err,m=:o,title="Interpolation error for equispaced nodes",
    xaxis=(L"n"),yaxis=(:log10,"max error"),)
```

The error initially decreases as one would expect but then begins to grow. Both phases occur at rates that are exponential in $n$, i.e., $O(K^n$) for a constant $K$, appearing linear on a semi-log plot.

```{raw} html
</div>
```

```{raw} latex
%%end demo%%
```

## Runge phenomenon

The disappointing loss of convergence in {numref}`Demo {number} <demo-stability-equispaced>` is a sign of ill conditioning due to the use of equally spaced nodes. We will examine this effect using the error formula {eq}`interperror` as a guide:

$$
f(x) - p(x) = \frac{f^{(n+1)}(\xi)}{(n+1)!} \Phi(x), \qquad \Phi(x) = \prod_{i=0}^n (x-t_i).
$$

While the dependence on $f$ is messy here, the error indicator $\Phi(x)$ can be studied as a function of the nodes only.

(demo-stability-errfun)=
```{proof:demo}
```

```{raw} html
<div class='demo'>
```

```{raw} latex
%%start demo%%
```

We plot $|\Phi(x)|$ over the interval $[-1,1]$ with equispaced nodes for different values of $n$. 

```{code-cell} 
:tags: [hide-input]
plot(xaxis=(L"x"),yaxis=(:log10,L"|\Phi(x)|",[1e-25,1]),legend=:bottomleft)

x = range(-1,1,length=2001)
for n in 10:10:50
    t = range(-1,1,length=n+1)
    Φ = [ prod(xₖ.-t) for xₖ in x ]
    scatter!(x,abs.(Φ),m=(1,stroke(0)),label="n=$n")
end
title!("Error indicator for equispaced nodes")
```

Each time $\Phi$ passes through zero at an interpolation node, the value on the log scale should go to $-\infty$, which explains the numerous cusps on the curves.

```{raw} html
</div>
```

```{raw} latex
%%end demo%%
```

Two observations from the result of {numref}`Demo {number} <demo-stability-errfun>` are important. First, $|\Phi|$ decreases exponentially at each fixed location in the interval (note that the spacing between curves is constant for constant increments of $n$). Second, $|\Phi|$ is larger at the ends of the interval than in the middle, by an exponentially growing factor. This gap is what can ruin the convergence of polynomial interpolation.

(demo-stability-runge)=

```{proof:demo}
```

```{raw} html
<div class='demo'>
```

```{raw} latex
%%start demo%%
```

This function has infinitely many continuous derivatives on the entire real line and looks easy to approximate over $[-1,1]$.

```{code-cell} 
f = x -> 1/(x^2 + 16)
plot(f,-1,1,title="Test function",legend=:none)
```

We start by doing equispaced polynomial interpolation for some small values of $n$.

```{code-cell} 
:tags: [hide-input]
plot(xaxis=(L"x"),yaxis=(:log10,L"|f(x)-p(x)|",[1e-20,1]))

x = range(-1,1,length=2501)
n = 4:4:12
for (k,n) in enumerate(n)
    t = range(-1,1,length=n+1)      # equally spaced nodes
    y = f.(t)                       # interpolation data
    p = FNC.polyinterp(t,y)
    err = @. abs(f(x)-p(x))
    plot!(x,err,m=(1,:o,stroke(0)),label="degree $n")
end
title!("Error for low degrees")
```

The convergence so far appears rather good, though not uniformly so. However, notice what happens as we continue to increase the degree.

```{code-cell} 
:tags: [hide-input]
n = @. 12 + 15*(1:3)
plot(xaxis=(L"x"),yaxis=(:log10,L"|f(x)-p(x)|",[1e-20,1]))

for (k,n) in enumerate(n)
    t = range(-1,1,length=n+1)      # equally spaced nodes
    y = f.(t)                       # interpolation data
    p = FNC.polyinterp(t,y)
    err = @. abs(f(x)-p(x))
    plot!(x,err,m=(1,:o,stroke(0)),label="degree $n")
end
title!("Error for higher degrees")
```

The convergence in the middle can't get any better than machine precision relative to the function values. So maintaining the growing gap between the center and the ends pushes the error curves upward exponentially fast at the ends, wrecking the convergence.
```{raw} html
</div>
```

```{raw} latex
%%end demo%%
```


```{index} ! Runge phenomenon
```

The observation of instability in {numref}`Demo {number} <demo-stability-runge>` is known as the **Runge phenomenon**. The Runge phenomenon is an instability manifested when the nodes of the interpolant are equally spaced and the degree of the polynomial increases. We reiterate that the phenomenon is rooted in the interpolation convergence theory and not a consequence of the algorithm chosen to implement polynomial interpolation.

Significantly, the convergence observed in {numref}`Demo {number} <demo-stability-runge>` is stable within a middle portion of the interval. By redistributing the interpolation nodes, we will next sacrifice a little of the convergence in the middle portion in order to improve it near the ends and rescue the process globally.

## Chebyshev nodes

```{index} ! Chebyshev points; second kind
```

The observations above hint that we might find success by having more nodes near the ends of the interval than in the middle. Though we will not give the details, it turns out that there is a precise asymptotic sense in which this must be done to make polynomial interpolation work over the entire interval. One especially important node family that gives stable convergence for polynomial interpolation is the **Chebyshev points of the second kind** (or *Chebyshev extreme points*) defined by

:::{math}
:label: chebextreme
 t_k = - \cos\left(\frac{k \pi}{n}\right), \qquad k=0,\ldots,n.
:::

These are the projections onto the $x$-axis of $n$ equally spaced points on a unit circle. They are densely clustered near the ends of $[-1,1]$, and this feature turns out to overcome the Runge phenomenon.

(demo-stability-errcheb)=
```{proof:demo}
```

```{raw} html
<div class='demo'>
```

```{raw} latex
%%start demo%%
```

Now we look at the error indicator function $\Phi$ for Chebyshev node sets.

```{code-cell} 
:tags: [hide-input]

plot(xaxis=(L"x"),yaxis=(:log10,L"|\Phi(x)|",[1e-18,1e-2]))
x = range(-1,1,length=2001)
for n in 10:10:50
    t = [ -cos(π*k/n) for k in 0:n ]                  
    Φ = [ prod(xₖ.-t) for xₖ in x ]
    plot!(x,abs.(Φ),m=(1,:o,stroke(0)),label="n=$n")
end
title!("Error indicator for Chebyshev nodes")
```

In contrast to the equispaced case, $|\Phi|$ decreases exponentially with $n$ almost uniformly across the interval.
```{raw} html
</div>
```

```{raw} latex
%%end demo%%
```

(demo-stability-rungefix)=
```{proof:demo}
```

```{raw} html
<div class='demo'>
```

```{raw} latex
%%start demo%%
```

Here again is the function from {numref}`Demo {number} <demo-stability-runge>` that provoked the Runge phenomenon when using equispaced nodes.

```{code-cell} 
f = x -> 1/(x^2 + 16);
```

```{code-cell} 
:tags: [hide-input]

plot(label="",xaxis=(L"x"),yaxis=(:log10,L"|f(x)-p(x)|",[1e-20,1]))
x = range(-1,1,length=2001)
for (k,n) in enumerate([4,10,16,40])
    t = [ -cos(pi*k/n) for k in 0:n ]
    y = f.(t)                           # interpolation data
    p = FNC.polyinterp(t,y)
    err = @.abs(f(x)-p(x))
    plot!(x,err,m=(1,:o,stroke(0)),label="degree $n")
end
title!("Error for Chebyshev interpolants")
```

By degree 16 the error is uniformly within machine epsilon, and, importantly, it stays there as $n$ increases. Note that as predicted by the error indicator function, the error is uniform over the interval at each value of $n$.
```{raw} html
</div>
```

```{raw} latex
%%end demo%%
```


As a bonus, for Chebyshev nodes the barycentric weights are simple:

:::{math}
:label: weightcheb
  w_k = (-1)^k d_k, \qquad d_k =
  \begin{cases}
    1/2 & \text{if $k=0$ or $k=n$},\\
    1 & \text{otherwise}.
  \end{cases}
:::

## Spectral convergence

If we take $n\rightarrow \infty$ and use polynomial interpolation on Chebyshev nodes, the convergence rate is exponential in $n$. The following is typical of the results that can be proved.

(theorem-spectral)=
::::{proof:theorem}
Suppose $f(x)$ is analytic in an open real interval containing $[-1,1]$. Then there exist constants $C>0$ and $K>1$ such that
  
:::{math}
:label: spectral
\max_{x\in[-1,1]} | f(x) - p(x) | \le C K^{-n},
:::

where $p$ is the unique polynomial of degree $n$ or less defined by interpolation on $n+1$ Chebyshev second-kind points.
::::

The condition "$f$ is analytic" means that the Taylor series of $f$ converges to $f(x)$ in an open interval containing $[-1,1]$.[^analytic] A necessary condition of analyticity is that $f$ is infinitely differentiable.

[^analytic]: Alternatively, analyticity means that the function is extensible to one that is differentiable in the complex plane.

```{index} ! convergence rate; spectral
```
In other contexts we refer to {eq}`spectral` as linear convergence, but here it is usual to say that the rate is exponential or that one has **spectral convergence**. It achieves constant reduction factors in the error by constant increments of $n$. By contrast, algebraic convergence in the form $O(n^{-p})$ for some $p>0$ requires *multiplying* $n$ by a constant factor in order to reduce error by a constant factor. Graphically, spectral error is a straight line on a log-linear scale, while algebraic convergence is a straight line on a log-log scale.

(demo-stability-spectral)=
:::{proof:demo}
:::

```{raw} html
<div class='demo'>
```

```{raw} latex
%%start demo%%
```

On the left, we use a log-log scale, which makes second-order algebraic convergence $O(n^{-4})$ a straight line. On the right, we use a log-linear scale, which makes spectral convergence $O(K^{-n})$ linear.

```{code-cell} 
:tags: [remove-cell]
disable_logging(Logging.Warn)
```

```{code-cell} 
:tags: [hide-input]
n = 20:20:400
algebraic = @. 100/n^4
spectral = @. 10*0.85^n
plot(n,[algebraic spectral],layout=(1,2),subplot=1,
    xaxis=(L"n",:log10),yaxis=(:log10,(1e-15,1)),
    label=["algebraic" "spectral"],title="Log-log")
plot!(n,[algebraic spectral],subplot=2,  
    xaxis=L"n",yaxis=(:log10,(1e-15,1)),
    label=["algebraic" "spectral"],title="log-linear")
```

```{raw} html
</div>
```

```{raw} latex
%%end demo%%
```

## Exercises

1. ⌨ Revisit {numref}`Demo %s <demo-stability-equispaced>` and determine an approximate value for the convergent phase of the constant $K$ mentioned in the comments there. 

2. ⌨ For each case, compute the polynomial interpolant using $n$ second-kind Chebyshev nodes in $[-1,1]$ for $n=4,8,12,\ldots,60$. At each value of $n$, compute the infinity-norm error (that is, $\max |p(x)-f(x)|$ evaluated for at least 4000 values of $x$). Using a log-linear scale, plot the error as a function of $n$, then determine a good approximation to the constant $K$ in {eq}`spectral`. 

    **(a)** $f(x) = 1/(25x^2+1)\qquad$ 
    **(b)** $f(x) = \tanh(5 x+2)$

    **(c)** $f(x) = \cosh(\sin x)\qquad$
    **(d)** $f(x) = \sin(\cosh x)$

    (problem-chebinterp)=
3. ⌨ Write a function `chebinterp(f,n)` that returns a function representing the polynomial interpolant of the input function `f` using $n+1$ Chebyshev second kind nodes over $[-1,1]$. You should use {eq}`weightcheb` to compute the barycentric weights directly, rather than using the method in {numref}`Function {number} <function-polyinterp>`. Test your function by revisiting {numref}`Demo %s <demo-stability-runge>` to use Chebyshev rather than equally spaced nodes. 

4. {numref}`Theorem %s <theorem-spectral>` assumes that the function being approximated has infinitely many derivatives over $[-1,1]$. But now consider the family of functions $f_m(x)=|x|^m$. 

    **(a)** ✍ How many continuous derivatives over $[-1,1]$ does $f_m$ possess?

    **(b)** ⌨ Compute the polynomial interpolant using $n$ second-kind Chebyshev nodes in $[-1,1]$ for $n=10,20,30,\ldots,100$. At each value of $n$, compute the max-norm error (that is, $\max |p(x)-f_m(x)|$ evaluated for at least 41000 values of $x$). Using a single log-log graph, plot the error as a function of $n$ for all six values $m=1,3,5,7,9,11$.
    
    **(c)** ✍  Based on the results of parts (a) and (b), form a hypothesis about the asymptotic behavior of the error for fixed $m$ as $n\rightarrow \infty$. 

    (problem-stability-changeinterval)=
5. The Chebyshev points can be used when the interval of interpolation is $[a,b]$ rather than $[-1,1]$ by means of the change of variable

    :::{math}
    :label: changeinterval
      z = \psi(x) = a + (b-a)\frac{(x+1)}{2}.
    :::

    **(a)** ✍  Show that $\psi(-1) = a$, $\psi(1) = b$, and $\psi$ is strictly increasing on $[-1,1]$.

    **(b)** ✍ Invert the relation {eq}`changeinterval` to solve for $x$ in terms of $\psi^{-1}(z)$. 

    **(c)** ✍ Let $t_0,\ldots,t_n$ be standard second-kind Chebyshev points. Then a polynomial in $x$ can be used to interpolate the function values $f\bigl(\psi(t_i)\bigr)$. This in turn implies an interpolating function $\tilde{P}(z) =P\bigl(\psi^{-1}(z)\bigr)$. Show that $\tilde{P}$ is a polynomial in $z$. 
	
    **(d)** ⌨ Implement the idea of part (c) to plot a polynomial interpolant of $f(z) =\cosh(\sin z)$ over $[0,2\pi]$ using $n+1$ Chebyshev nodes with $n=40$. 

6. The Chebyshev points can be used for interpolation of functions defined on the entire real line by using the change of variable

    :::{math}
    :label: changeintervalinf
    z = \phi(x) = \frac{2x}{1-x^2},
    :::

    which maps the interval $(-1,1)$ in one-to-one fashion to the entire real line.

    **(a)** ✍ Find $\displaystyle \lim_{x\to 1^-} \phi(x)$ and $\displaystyle \lim_{x\to -1^+} \phi(x)$. 

    **(b)** ✍ Invert {eq}`changeintervalinf` to express $x=\phi^{-1}(z)$. (Be sure to enforce $-1\le x \le 1$.)

    **(c)** ⌨ Let $t_0,\ldots,t_n$ be standard second-kind Chebyshev points. These map to the $z$ variable as $\zeta_i=\phi(t_i)$ for all $i$. Suppose that $f(z)$ is a given function whose domain is the entire real line. Then the function values $y_i=f(\zeta_i)$ can be associated with the Chebyshev nodes $t_i$, leading to a polynomial interpolant $p(x)$. This in turn implies an interpolating function on the real line, defined as

    $$
    q(z)=p\bigl(\phi^{-1}(z)\bigr) = p(x).
    $$
    
    Implement this idea to plot an interpolant of $f(z)=(z^2-2z+2)^{-1}$ using $n=30$. Your plot should show $q(z)$ evaluated at 1000 evenly spaced points in $[-6,6]$, with markers at the nodal values (those lying within the $[-6,6]$ window).
