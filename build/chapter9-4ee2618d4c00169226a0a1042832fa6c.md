---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
# Chapter 9

## Functions

(function-polyinterp-julia)=
``````{dropdown} Barycentric polynomial interpolation
```{literalinclude} FNCFunctions/src/chapter09.jl
:filename: polyinterp.jl
:start-after: # begin polyinterp
:end-before: # end polyinterp
:language: julia
:linenos: true
```

````{admonition} About the code
:class: dropdown
As noted in {numref}`Example %s <example-writeoutbary2>`, a common scaling factor in the weights does not affect the barycentric formula {eq}`bary2`. In lines 9--10 this fact is used to rescale the nodes in order to avoid eventual tiny or enormous numbers that could go outside the bounds of double precision.

The return value is a function that evaluates the polynomial interpolant. Within this function, `isinf` is used to detect either `Inf` or `-Inf`, which occurs when $x$ exactly equals one of the nodes. In this event, the corresponding data value is returned.
````
``````

(function-triginterp-julia)=
``````{dropdown} Trigonometric interpolation
```{literalinclude} FNCFunctions/src/chapter09.jl
:filename: triginterp.jl
:start-after: # begin triginterp
:end-before: # end triginterp
:language: julia
:linenos: true
```

````{admonition} About the code
:class: dropdown
The construct on line 13 is known as a *ternary operator*. It is a shorthand for an `if`–`else` statement, giving two alternative results for the true/false cases. Line 19 uses `eachindex(y)`, which generalizes `1:length(y)` to cases where a vector might have a more exotic form of indexing.
````
``````

(function-ccint-julia)=
``````{dropdown} Clenshaw–Curtis integration
```{literalinclude} FNCFunctions/src/chapter09.jl
:filename: ccint.jl
:start-after: # begin ccint
:end-before: # end ccint
:language: julia
:linenos: true
```
``````

(function-glint-julia)=
``````{dropdown} Gauss–Legendre integration
```{literalinclude} FNCFunctions/src/chapter09.jl
:filename: glint.jl
:start-after: # begin glint
:end-before: # end glint
:language: julia
:linenos: true
```
``````

(function-intinf-julia)=
``````{dropdown} Integration over $(-\infty,\infty)$
```{literalinclude} FNCFunctions/src/chapter09.jl
:filename: intinf.jl
:start-after: # begin intinf
:end-before: # end intinf
:language: julia
:linenos: true
```

::::{admonition} About the code
:class: dropdown
The test `isinf(x(M))` in line 17 checks whether $x(M)$ is larger than the maximum double-precision value, causing it to *overflow* to `Inf`.
::::
``````

(function-intsing-julia)=
``````{dropdown} Integration with endpoint singularities
```{literalinclude} FNCFunctions/src/chapter09.jl
:filename: intsing.jl
:start-after: # begin intsing
:end-before: # end intsing
:language: julia
:linenos: true
```

::::{admonition} About the code
:class: dropdown
The test `iszero(x(M))` in line 17 checks whether $x(M)$ is less than the smallest positive double-precision value, causing it to *underflow* to zero.
::::
``````

## Examples

```{code-cell}
:tags: remove-output
include("FNC_init.jl")
```

### 9.1 @section-globalapprox-polynomial

(demo-polynomial-lagrange-julia)=
``````{dropdown} @demo-polynomial-lagrange
Here is a vector of nodes.

```{code-cell}
t = [ 1, 1.5, 2, 2.25, 2.75, 3 ]
n = length(t) - 1;
```

::::{grid} 1 1 2 2
Let's apply the definition of the cardinal Lagrange polynomial for $k=2$. First we define a polynomial $q$ that is zero at all the nodes except $i=k$. Then $\ell_2$ is found by normalizing $q$ by $q(t_k)$.
:::{card}
Character ℓ is typed as `\ell`<kbd>Tab</kbd>.
:::
::::

```{code-cell}
k = 2
q(x) = prod(x - t[i] for i in [0:k-1; k+1:n] .+ 1)
ℓₖ(x) = q(x) / q(t[k+1]);
```

A plot confirms the cardinal property of the result.

```{code-cell}
using Plots
plot(ℓₖ, 1, 3)
y = zeros(n+1);  y[k+1] = 1
scatter!(t, y, color=:black,
    xaxis=(L"x"),  yaxis=(L"\ell_2(x)"),
    title="Lagrange cardinal function")
```

Observe that $\ell_k$ is _not_ between zero and one everywhere, unlike a hat function.
``````

(demo-polynomial-error-julia)=
``````{dropdown} @demo-polynomial-error
Consider the problem of interpolating $\log(x)$ at these nodes:

```{code-cell}
t =  [ 1, 1.6, 1.9, 2.7, 3 ]
n = length(t) - 1;
```

::::{grid} 1 1 2 2
Here $n=4$ and $f^{(5)}(\xi) = 4!/\xi^5$. For $\xi\in[1,3]$ we can say that $|f^{(5)}(\xi)| \le 4!$. Hence 

$$ |f(x)-p(x)| \le \frac{1}{5} \left| \Phi(x) \right|.$$
:::{card}
Character Φ is typed as `\Phi`<kbd>Tab</kbd>. (Note the capitalization.)
:::
::::

```{code-cell}
using Polynomials
Φ(x) = prod(x - tᵢ for tᵢ in t)
plot(x -> 0.2 * abs(Φ(x)), 1, 3, label=L"\frac{1}{5}|\Phi(t)|")
p = Polynomials.fit(t, log.(t))
plot!(t -> abs(log(t) - p(t)), 1, 3, label=L"|f(x)-p(x)|")
scatter!(t, zeros(size(t)), color=:black,
    xaxis=(L"x"), title="Interpolation error and upper bound")
```

The error is zero at the nodes, by the definition of interpolation. The error bound, as well as the error itself, has one local maximum between each consecutive pair of nodes.
``````

### 9.2 @section-globalapprox-barycentric

(demo-barycentric-example-julia)=
``````{dropdown} @demo-barycentric-example
```{code-cell}
using Plots
f(x) = sin(exp(2x))
plot(f, 0, 1, label="function", legend=:bottomleft)
```

```{code-cell}
t = (0:3) / 3
y = f.(t)
scatter!(t, y, color=:black, label="nodes")
```

```{code-cell}
p = FNC.polyinterp(t, y)
plot!(p, 0, 1, label="interpolant", title="Interpolation on 4 nodes")
```

The curves must intersect at the interpolation nodes. For $n=6$ the interpolant is noticeably better.

```{code-cell}
plot(f, 0, 1, label="function", legend=:bottomleft)
t = (0:6) / 6
y = f.(t)
p = FNC.polyinterp(t, y)
scatter!(t, y, color=:black, label="nodes")
plot!(p, 0, 1, label="interpolant", title="Interpolation on 7 nodes")
```
``````

### 9.3 @section-globalapprox-stability

(demo-stability-equispaced-julia)=
``````{dropdown} @demo-stability-equispaced
We choose a function over the interval $[0,1]$. 

```{code-cell} 
f(x) = sin(exp(2x));
```

Here is a graph of $f$ and its polynomial interpolant using seven equally spaced nodes.

```{code-cell} 
:tags: hide-input
using Plots
plot(f, 0, 1, label="function", legend=:bottomleft)
t = range(0, 1, 7)    # 7 equally spaced nodes
y = f.(t)
scatter!(t, y, label="nodes")

p = FNC.polyinterp(t, y)
plot!(p, 0, 1, label="interpolant", title="Equispaced interpolant, n=6")
```

This looks pretty good. We want to track the behavior of the error as $n$ increases. We will estimate the error in the continuous interpolant by sampling it at a large number of points and taking the max-norm.

```{code-cell} 
:tags: hide-input
n = 5:5:60
err = zeros(size(n))
x = range(0, 1, 2001)             # for measuring error
for (k, n) in enumerate(n)
    t = range(0, 1, n+1)          # equally spaced nodes
    y = f.(t)                     # interpolation data
    p = FNC.polyinterp(t, y)
    err[k] = norm((@. f(x) - p(x)), Inf)
end
plot(n, err, m=:o, 
    xaxis=(L"n"), yaxis=(:log10, "max error"),
    title="Interpolation error for equispaced nodes")
```

The error initially decreases as one would expect but then begins to grow. Both phases occur at rates that are exponential in $n$, i.e., $O(K^n$) for a constant $K$, appearing linear on a semi-log plot.
``````

(demo-stability-errfun-julia)=
``````{dropdown} @demo-stability-errfun
We plot $|\Phi(x)|$ over the interval $[-1,1]$ with equispaced nodes for different values of $n$. 

```{code-cell} 
:tags: hide-input
plot(xaxis=(L"x"), yaxis=(:log10, L"|\Phi(x)|", [1e-25, 1]), legend=:bottomleft)
x = range(-1, 1, 2001)
for n in 10:10:50
    t = range(-1, 1, n+1)
    Φ(x) = prod(x - t for t in t)
    scatter!(x, abs.(Φ.(x)), m=(1, stroke(0)), label="n=$n")
end
title!("Error indicator for equispaced nodes")
```

Each time $\Phi$ passes through zero at an interpolation node, the value on the log scale should go to $-\infty$, which explains the numerous cusps on the curves.
``````

(demo-stability-runge-julia)=
``````{dropdown} @demo-stability-runge
This function has infinitely many continuous derivatives on the entire real line and looks easy to approximate over $[-1,1]$.

```{code-cell} 
f(x) = 1 / (x^2 + 16)
plot(f, -1, 1, title="Test function", legend=:none)
```

We start by doing equispaced polynomial interpolation for some small values of $n$.

```{code-cell} 
:tags: hide-input
plot(xaxis=(L"x"), yaxis=(:log10, L"|f(x)-p(x)|", [1e-20, 1]))
x = range(-1, 1, 2501)
n = 4:4:12
for (k, n) in enumerate(n)
    t = range(-1, 1, n+1)           # equally spaced nodes
    y = f.(t)                       # interpolation data
    p = FNC.polyinterp(t, y)
    err = @. abs(f(x) - p(x))
    plot!(x, err, m=(1, :o, stroke(0)), label="degree $n")
end
title!("Error for low degrees")
```

The convergence so far appears rather good, though not uniformly so. However, notice what happens as we continue to increase the degree.

```{code-cell} 
:tags: hide-input
n = @. 12 + 15 * (1:3)
plot(xaxis=(L"x"), yaxis=(:log10, L"|f(x)-p(x)|", [1e-20, 1]))
for (k, n) in enumerate(n)
    t = range(-1, 1, n+1)           # equally spaced nodes
    y = f.(t)                       # interpolation data
    p = FNC.polyinterp(t, y)
    err = @. abs(f(x) - p(x))
    plot!(x, err, m=(1, :o, stroke(0)), label="degree $n")
end
title!("Error for higher degrees")
```

The convergence in the middle can't get any better than machine precision relative to the function values. So maintaining the growing gap between the center and the ends pushes the error curves upward exponentially fast at the ends, wrecking the convergence.
``````

(demo-stability-errcheb-julia)=
``````{dropdown} @demo-stability-errcheb
Now we look at the error indicator function $\Phi$ for Chebyshev node sets.

```{code-cell} 
:tags: hide-input
plot(xaxis=(L"x"), yaxis=(:log10, L"|\Phi(x)|", [1e-18, 1e-2]))
x = range(-1, 1, 2001)
for n in 10:10:50
    t = [-cospi(k / n) for k in 0:n]
    Φ(x) = prod(x - t for t in t)
    plot!(x, abs.(Φ.(x)), m=(1, :o, stroke(0)), label="n=$n")
end
title!("Error indicator for Chebyshev nodes")
```

In contrast to the equispaced case, $|\Phi|$ decreases exponentially with $n$ almost uniformly across the interval.
``````

(demo-stability-rungefix-julia)=
``````{dropdown} @demo-stability-rungefix
Here again is the function from {numref}`Demo {number} <demo-stability-runge>` that provoked the Runge phenomenon when using equispaced nodes.

```{code-cell} 
f(x) = 1 / (x^2 + 16);
```

```{code-cell} 
:tags: hide-input

plot(label="", xaxis=(L"x"), yaxis=(:log10, L"|f(x)-p(x)|", [1e-20, 1]))
x = range(-1, 1, 2001)
for (k, n) in enumerate([4, 10, 16, 40])
    t = [-cospi(k / n) for k in 0:n]
    y = f.(t)                           # interpolation data
    p = FNC.polyinterp(t, y)
    err = @. abs(f(x) - p(x))
    plot!(x, err, m=(1, :o, stroke(0)), label="degree $n")
end
title!("Error for Chebyshev interpolants")
```

By degree 16 the error is uniformly within machine epsilon, and, importantly, it stays there as $n$ increases. Note that as predicted by the error indicator function, the error is uniform over the interval at each value of $n$.
``````

(demo-stability-spectral-julia)=
``````{dropdown} @demo-stability-spectral
```{code-cell} 
:tags: remove-cell
using Logging
disable_logging(Logging.Warn);
```
On the left, we use a log-log scale, which makes second-order algebraic convergence $O(n^{-4})$ a straight line. On the right, we use a log-linear scale, which makes spectral convergence $O(K^{-n})$ linear.

```{code-cell} 
:tags: hide-input
n = 20:20:400
algebraic = @. 100 / n^4
spectral = @. 10 * 0.85^n
plot(n, [algebraic spectral], layout=(1, 2), subplot=1,
    xaxis=(L"n", :log10),  yaxis=(:log10, (1e-15, 1)),
    label=["algebraic" "spectral"],  title="Log-log")
plot!(n, [algebraic spectral], subplot=2,
    xaxis=L"n",  yaxis=(:log10, (1e-15, 1)),
    label=["algebraic" "spectral"],  title="log-linear")
```
``````

### 9.4 @section-globalapprox-orthogonal

(demo-orthogonal-approx-julia)=
``````{dropdown} @demo-orthogonal-approx
Let's approximate $e^x$ over the interval $[−1,1]$. We can sample it at, say, 15 points, and find the best-fitting straight line to that data.

```{code-cell}
using Plots
plot(exp, -1, 1, label="function")
t = range(-1, 1, 15)
y = exp.(t)
V = [ti^j for ti in t, j in 0:1]  # Vandermonde-ish
c = V \ y
plot!(t -> c[1] + c[2] * t, -1, 1;
    label="linear fit for 15 points", legend=:bottomright,
    xaxis=("x"),  yaxis=("value"),
    title="Least-squares fit of exp(x)")
```

There's nothing special about 20 points. Choosing more doesn't change the result much.

```{code-cell}
t = range(-1, 1, 150)
y = exp.(t)
V = [ ti^j for ti in t, j=0:1 ]
c = V \ y
plot!(t -> c[1] + c[2]*t, -1, 1,
    label="linear fit for 150 points",  legend=:bottomright,
    xaxis=("x"),  yaxis=("value"),
    title="Least-squares fit of exp(x)")
```

This situation is unlike interpolation, where the degree of the interpolant increases with the number of nodes. Here, the linear fit is apparently approaching a limit that we may think of as a continuous least-squares fit.

```{code-cell}
n = 40:60:400
slope = zeros(size(n))
intercept = zeros(size(n))

for (k, n) in enumerate(n)
    t = range(-1, 1, n)
    y = exp.(t)
    V = [ ti^j for ti in t, j in 0:1 ]
    c = V \ y
    intercept[k], slope[k] = c
end

labels = ["n", "intercept", "slope"]
@pt :header=labels, [n intercept slope]
```
``````

### 9.5 @section-globalapprox-trig

(demo-trig-interp-julia)=
``````{dropdown} @demo-trig-interp
::::{grid} 1 1 2 2
We will get a cardinal function without using an explicit formula, just by passing data that is 1 at one node and 0 at the others.
:::{card}
The operator `÷`, typed as `\div` then <kbd>Tab</kbd>, returns the quotient without remainder of two integers.
:::
::::

```{code-cell}
using Plots
N = 7
n = (N - 1) ÷ 2
t = 2 * (-n:n) / N
y = zeros(N)
y[n+1] = 1

p = FNC.triginterp(t, y);
plot(p, -1, 1)

scatter!(t, y, color=:black, 
    xaxis=(L"x"),  yaxis=(L"\tau(x)"),
    title="Trig cardinal function, N=$N")
```

Here is a 2-periodic function and one of its interpolants.

```{code-cell}
f(x) = exp( sinpi(x) - 2*cospi(x) )
y = f.(t)
p = FNC.triginterp(t, y)

plot(f, -1, 1, label="function",
    xaxis=(L"x"),  yaxis=(L"p(x)"),
    title="Trig interpolation, N=$N", legend=:top)
scatter!(t, y, m=:o, color=:black, label="nodes")
plot!(p, -1, 1, label="interpolant")
```

```{index} ! Julia; ÷
```

The convergence of the interpolant is spectral. We let $N$ go needlessly large here in order to demonstrate that unlike polynomials, trigonometric interpolation is stable on equally spaced nodes. Note that when $N$ is even, the value of $n$ is not an integer but works fine for defining the nodes.

```{code-cell}
:tags: hide-input
N = 2:2:60
err = zeros(size(N))
x = range(-1, 1, 2501)  # for measuring error
for (k,N) in enumerate(N)
    n = (N-1) / 2;   t = 2*(-n:n) / N;
    p = FNC.triginterp(t, f.(t))
    err[k] = norm(f.(x) - p.(x), Inf)
end

plot(N, err, m=:o,
    xaxis=(L"N"),  yaxis=(:log10, "max error"),
    title="Convergence of trig interpolation")
```
``````

(demo-trig-fft-julia)=
``````{dropdown} @demo-trig-fft
This function has frequency content at $2\pi$, $-2\pi$, and $\pi$. 

```{code-cell}
f(x) = 3 * cospi(2x) - cispi(x)    # cispi(x) := exp(1im * π * x)
```

To use `fft`, we set up nodes in the interval $[0,2)$. 

```{code-cell}
n = 4;  N = 2n+1;
t = [ 2j / N for j in 0:N-1 ]      # nodes in [0,2)
y = f.(t);
```

We perform Fourier analysis using `fft` and then examine the resulting coefficients.

```{code-cell}
using FFTW
c = fft(y) / N
freq = [0:n; -n:-1]
@pt :header=["k", "coefficient"] [freq round.(c, sigdigits=5)]
```

Note that $1.5 e^{2i\pi x}+1.5 e^{-2i\pi x} = 3 \cos(2\pi x)$, so this result is sensible.

Fourier's greatest contribution to mathematics was to point out that *every* periodic function is just a combination of frequencies—infinitely many of them in general, but truncated for computational use. Here we look at the magnitudes of the coefficients for $f(x) = \exp( \sin(\pi x) )$.

```{code-cell}
:tags: hide-input
f(x) = exp( sin(pi*x) )     # content at all frequencies
n = 9;  N = 2n+1;
t = [ 2j / N for j in 0:N-1 ]      # nodes in [0,2)
c = fft(f.(t)) / N

freq = [0:n; -n:-1]
scatter(freq, abs.(c);
    xaxis=(L"k", [-n, n]),  yaxis=(L"|c_k|", :log10), 
    title="Fourier coefficients",  legend=:none)
```

The Fourier coefficients of smooth functions decay exponentially in magnitude as a function of the frequency. This decay rate is determines the convergence of the interpolation error.
``````

### 9.6 @section-globalapprox-integration

(demo-integration-ellipse-julia)=
``````{dropdown} @demo-integration-ellipse
```{code-cell}
f(t) = π * sqrt( cospi(t)^2 + sinpi(t)^2 / 4 );
n = 4:4:48
perim = zeros(size(n))
for (k, n) in enumerate(n)
    h = 2 / n
    t = @. h * (0:n-1) - 1
    perim[k] = h * sum(f.(t))
end
err = @. abs(perim - perim[end])    # use last value as "exact"
@ptconf formatters=ft_printf(["%d", "%.15f", "%.2e"], 1:3)
@pt :header=["n", "perimeter", "error estimate"] [n perim err][1:end-1, :]
```
The approximations gain about one digit of accuracy for each constant increment of $n$, which is consistent with spectral convergence.
``````

(demo-integration-compare-julia)=
``````{dropdown} @demo-integration-compare
First consider the integral 

$$
\int_{-1}^1 \frac{1}{1+4x^2} \, dx = \arctan(2).
$$

```{code-cell}
f(x)= 1 / (1 + 4x^2);
exact = atan(2);
```

We compare the two spectral integration methods for a range of $n$ values.

```{code-cell}
:tags: hide-input
using Plots
n = 8:4:96
err = zeros(length(n), 2)
for (k, n) in enumerate(n)
  err[k, 1] = abs(exact - FNC.ccint(f, n)[1])
  err[k, 2] = abs(exact - FNC.glint(f, n)[1])
end

err[iszero.(err)] .= NaN    # remove from log-scale plot
plot(n, err, m=:o, label=["CC" "GL"],
    xaxis=("number of nodes"),  yaxis=(:log10, "error", [1e-16, 1]), 
    title="Spectral integration")
```

(The missing dots are where the error is exactly zero.) Gauss–Legendre does converge faster here, but at something less than twice the rate.

Now we try a more sharply peaked integrand:
 
 $$\int_{-1}^1 \frac{1}{1+16x^2} \, dx = \frac{1}{2}\arctan(4).$$ 

```{code-cell}
f(x) = 1 / (1 + 16x^2);
exact = atan(4) / 2;
```

```{code-cell}
:tags: hide-input
n = 8:4:96
err = zeros(length(n), 2)
for (k,n) in enumerate(n)
  err[k, 1] = abs(exact - FNC.ccint(f, n)[1])
  err[k, 2] = abs(exact - FNC.glint(f, n)[1])
end

err[iszero.(err)] .= NaN    # remove from log-scale plot
plot(n, err, m=:o, label=["CC" "GL"],
    xaxis=("number of nodes"),  yaxis=(:log10, "error", [1e-16, 1]), 
    title="Spectral integration")
```

The two are very close until about $n=40$, when the Clenshaw–Curtis method slows down.

Now let's compare the spectral performance to that of our earlier adaptive method in `intadapt`. We will specify varying error tolerances and record the error as well as the total number of evaluations of $f$.

```{code-cell}
:tags: hide-input
tol = 10 .^(-2.0:-2:-14)
n = zeros(size(tol))  
errAdapt = zeros(size(tol))
for (k, tol) in enumerate(tol)
  Q, t = FNC.intadapt(f, -1, 1, tol)
  errAdapt[k] = abs(exact - Q)
  n[k] = length(t)
end

errAdapt[iszero.(errAdapt)] .= NaN
plot!(n, errAdapt, m=:o, label="intadapt")
plot!(n, n.^(-4), l=:dash, label="4th order",
        xaxis=(:log10),  title="Spectral vs 4th order" )
```

At the core of `intadapt` is a fourth-order formula, and the results track that rate closely. For all but the most relaxed error tolerances, both spectral methods are far more efficient than the low-order counterpart. For other integrands, particularly those that vary nonuniformly across the interval, the adaptive method might be more competitive.

``````

### 9.7 @section-globalapprox-improper

(demo-improper-decay-julia)=
``````{dropdown} @demo-improper-decay
```{code-cell}
:tags: hide-input
using Plots
f(x) = 1 / (1 + x^2)
plot(f, -4, 4, layout=(2, 1),
    xlabel=L"x", 
    yaxis=(:log10, L"f(x)", (1e-16, 2)),
    title="Original integrand")

ξ(t) = sinh( π * sinh(t) / 2 )
dξ_dt(t) = π/2 * cosh(t) * cosh(π * sinh(t) / 2)
g(t) = f(ξ(t)) * dξ_dt(t)

plot!(g,-4, 4, subplot=2,
    xlabel=L"t",
    yaxis=(:log10, L"f(x(t))\cdot x'(t)", (1e-16, 2)),
    title="Transformed integrand")
```

This graph suggests that we capture all of the integrand values that are larger than machine epsilon by integrating in $t$ from $-4$ to $4$.
``````

(demo-improper-intinf-julia)=
``````{dropdown} @demo-improper-intinf
```{code-cell}
:tags: hide-input
f(x) = 1 / (1 + x^2)
tol = [1 / 10^d for d in 5:0.5:14]
err = zeros(length(tol), 2)
len = zeros(Int, length(tol), 2)
for (i, tol) in enumerate(tol)
    I1, x1 = FNC.intadapt(f, -2/tol, 2/tol, tol)
    I2, x2 = FNC.intinf(f, tol)
    @. err[i,:] = abs(π - [I1, I2])
    @. len[i,:] = length([x1, x2])
end
plot(len, err, m=:o, label=["direct" "double exponential"])
n = [100, 10000]
plot!(n, 1000n.^(-4), 
    color=:black,  l=:dash,
    label="fourth-order",  legend=:bottomleft,
    xaxis=(:log10, "number of nodes"), 
    yaxis=(:log10, "error"),
    title="Comparison of integration methods")
```

Both methods are roughly fourth-order due to Simpson's formula in the underlying adaptive integration method. At equal numbers of evaluation nodes, however, the double exponential method is consistently 2--3 orders of magnitude more accurate.
``````

(demo-improper-intsing-julia)=
``````{dropdown} @demo-improper-intsing

```{code-cell}
:tags: hide-input
f(x) = 1 / (10 * sqrt(x))
tol = [1 / 10^d for d in 5:0.5:14]
err = zeros(length(tol), 2)
len = zeros(Int, length(tol), 2)
for (i, tol) in enumerate(tol)
    I1, x1 = FNC.intadapt(f, (tol/20)^2, 1, tol)
    I2, x2 = FNC.intsing(f, tol)
    @. err[i, :] = abs(0.2 - [I1, I2])
    @. len[i, :] = length([x1, x2])
end
plot(len, err, m=:o, label=["direct" "double exponential"])
n = [30, 3000]
plot!(n, 30n.^(-4);
    color=:black,  l=:dash,
    label="fourth-order",  legend=:bottomleft,
    xaxis=(:log10, "number of nodes"),
    yaxis=(:log10, "error"),
    title="Comparison of integration methods")
```

As in {numref}`Demo {number} <demo-improper-intinf>`, the double exponential method is more accurate than direct integration by a few orders of magnitude. Equivalently, the same accuracy can be reached with many fewer nodes.
``````

