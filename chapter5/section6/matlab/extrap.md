---
kernelspec:
  display_name: MATLAB
  language: matlab
  name: jupyter_matlab_kernel
numbering:
  headings: false
---
```{code-cell}
:tags: [remove-cell]
cd  /Users/driscoll/Documents/GitHub/fnc/matlab
FNC_init
```
[**Demo %s**](#demo-int-extrap)

We estimate $\displaystyle\int_0^2 x^2 e^{-2x}\, dx$ using extrapolation. First we use `quadgk` to get an accurate value.

```{code-cell}
f = @(x) x.^2 .* exp(-2 * x);
a = 0;  b = 2;
format long
I = integral(f, a, b, abstol=1e-14, reltol=1e-14)
```

We start with the trapezoid formula on $n=N$ nodes.

```{code-cell}
N = 20;       % the coarsest formula
n = N;  h = (b - a) / n;
t = h * (0:n)';
y = f(t);
```

We can now apply weights to get the estimate $T_f(N)$.

```{code-cell}
T = h * ( sum(y(2:n)) + y(1) / 2 + y(n+1) / 2 )
```

Now we double to $n=2N$, but we only need to evaluate $f$ at every other interior node and apply {eq}`nc-doubling`.

```{code-cell}
n = 2*n;  h = h / 2;
t = h * (0:n)';
T(2) = T(1) / 2 + h * sum( f(t(2:2:n)) )
```

We can repeat the same code to double $n$ again.

```{code-cell}
n = 2*n;  h = h / 2;
t = h * (0:n)';
T(3) = T(2) / 2 + h * sum( f(t(2:2:n)) )
```

Let us now do the first level of extrapolation to get results from Simpson's formula. We combine the elements `T[i]` and `T[i+1]` the same way for $i=1$ and $i=2$.

```{code-cell}
S = (4 * T(2:3) - T(1:2)) / 3
```

With the two Simpson values $S_f(N)$ and $S_f(2N)$ in hand, we can do one more level of extrapolation to get a sixth-order accurate result.

```{code-cell}
R = (16*S(2) - S(1)) / 15
```

We can make a triangular table of the errors:

```{code-cell}
err2 = T(:) - I;
err4 = [NaN; S(:) - I];
err6 = [NaN; NaN; R - I];
format short e
disp(table(err2, err4, err6, variablenames=["order 2", "order 4", "order 6"]))
```

If we consider the computational time to be dominated by evaluations of $f$, then we have obtained a result with about twice as many accurate digits as the best trapezoid result, at virtually no extra cost.
