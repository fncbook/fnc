---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering: 
  headings: false
---
```{code-cell}
include("../../../julia/FNC_init.jl")
```
[**Demo %s**](#demo-int-extrap)

We estimate $\displaystyle\int_0^2 x^2 e^{-2x}\, dx$ using extrapolation. First we use `quadgk` to get an accurate value.

```{code-cell}
using QuadGK
f = x -> x^2 * exp(-2x);
a = 0;
b = 2;
Q, _ = quadgk(f, a, b, atol=1e-14, rtol=1e-14)
@show Q;
```

We start with the trapezoid formula on $n=N$ nodes.

```{code-cell}
N = 20;       # the coarsest formula
n = N;
h = (b - a) / n;
t = h * (0:n);
y = f.(t);
```

We can now apply weights to get the estimate $T_f(N)$.

```{code-cell}
T = [h * (sum(y[2:n]) + y[1] / 2 + y[n+1] / 2)]
```

Now we double to $n=2N$, but we only need to evaluate $f$ at every other interior node and apply {eq}`nc-doubling`.

```{code-cell}
n = 2n;
h = h / 2;
t = h * (0:n);
T = [T; T[end] / 2 + h * sum(f.(t[2:2:n]))]
```

We can repeat the same code to double $n$ again.

```{code-cell}
n = 2n;
h = h / 2;
t = h * (0:n);
T = [T; T[end] / 2 + h * sum(f.(t[2:2:n]))]
```

Let us now do the first level of extrapolation to get results from Simpson's formula. We combine the elements `T[i]` and `T[i+1]` the same way for $i=1$ and $i=2$.

```{code-cell}
S = [(4T[i+1] - T[i]) / 3 for i in 1:2]
```

With the two Simpson values $S_f(N)$ and $S_f(2N)$ in hand, we can do one more level of extrapolation to get a sixth-order accurate result.

```{code-cell}
R = (16S[2] - S[1]) / 15
```

We can make a triangular table of the errors:
```{tip}
:class: dropdown
The value `nothing` equals nothing except `nothing`.
```

```{code-cell}
err = [T .- Q [nothing; S .- Q] [nothing; nothing; R - Q]]
@pt :header=["order 2", "order 4", "order 6"] err
```

If we consider the computational time to be dominated by evaluations of $f$, then we have obtained a result with about twice as many accurate digits as the best trapezoid result, at virtually no extra cost.
