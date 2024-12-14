---
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

```{code-cell} ipython3
import numpy as np
import FNC
```

(demo-float-accuracy-python)=
``````{dropdown} Floating-point accuracy
Recall the grade-school approximation to the number $\pi$.

```{code-cell} ipython3
p = 22/7
print(p)
```
::::{grid} 1 1 2 2

Not all the digits displayed for `p` are the same as those of $\pi$. 

:::{card}
:columns: 7

The value of `pi` is predefined in the `numpy` package, which is typically abbreviated as `np`.
:::
::::

```{code-cell} ipython3
print(np.pi)
```

::::{grid} 1 1 2 2

The absolute and relative accuracies of the approximation are as follows:

:::{card}
:columns: 7

We often use [Python f-strings](https://docs.python.org/3/tutorial/inputoutput.html#tut-f-strings) to format numerical output. 

:::
::::


```{code-cell} ipython3
print(f"absolute accuracy: {abs(p - np.pi)}")
```

```{code-cell} ipython3
rel_acc = abs(p - np.pi) / np.pi
print("relative accuracy: {rel_acc:.4e}")
```

::::{grid} 1 1 2 2
Here we calculate the number of accurate digits in `p`:
:::{card}
The `np.log` function is for the natural log. For other common bases, use `log10` or `log2`.
:::
::::


```{code-cell} ipython3
print(f"accurate digits: {-np.log10(rel_acc):.1f}")
```
``````

(demo-float-python)=
``````{dropdown} Floating-point representation
:open: false

Python has native `int` and `float` types.

```{code-cell} ipython3
print(f"The type of {1} is {type(1)}")
print(f"The type of {float(1)} is {type(1.0)}")
```

The `numpy` package has its own `float` types:

```{code-cell} ipython3
one = np.float64(1)
print(f"The type of {one} is {type(one)}")
```

Both `float` and `np.float64` are double precision, using 64 binary bits per value. Although it is not normally necessary to do so, we can deconstruct a float into its significand and exponent:

```{code-cell} ipython3
significand, exponent = np.frexp(one)
print(f"significand: {significand}, exponent: {exponent}")
```

```{code-cell} ipython3
significand, exponent = np.frexp(one / 8)
print(f"significand: {significand}, exponent: {exponent}")
```

The spacing between floating-point values in $[2^n,2^{n+1})$ is $2^n \epsilon_\text{mach}$, where $\epsilon_\text{mach}$ is machine epsilon, given here for double precision:

```{code-cell} ipython3
mach_eps = np.finfo(float).eps
print(f"machine epsilon is {mach_eps:.4e}")
```

Because double precision allocates 52 bits to the significand, the default value of machine epsilon is $2^{-52}$.

```{code-cell} ipython3
print(f"machine epsilon is 2 to the power {np.log2(mach_eps)}")
```

A common mistake is to think that $\epsilon_\text{mach}$ is the smallest floating-point number. It's only the smallest *relative to 1*. The correct perspective is that the scaling of values is limited by the exponent, not the significand. The actual range of positive values in double precision is

```{code-cell}
finf = np.finfo(float)
print(f"range of positive values: [{finf.tiny}, {finf.max}]")
```

For the most part you can mix integers and floating-point values and get what you expect.

```{code-cell}
1/7
```

```{code-cell}
37.3 + 1
```

```{code-cell}
2**(-4)
```

You can convert a floating value to an integer by wrapping it in `int`.

```{code-cell}
int(3.14)
```
``````

(demo-float-arithmetic-python)=
``````{dropdown} Floating-point arithmetic oddity
:open: false

There is no double precision number between $1$ and $1+\varepsilon_\text{mach}$. Thus, the following difference is zero despite its appearance.

```{code-cell} ipython3
eps = np.finfo(float).eps
e = eps/2
print((1.0 + e) - 1.0)
```

However, $1-\varepsilon_\text{mach}/2$ is a double precision number, so it and its negative are represented exactly:

```{code-cell} ipython3
print(1.0 + (e - 1.0))
```

This is now the "correct" result. But we have found a rather shocking breakdown of the associative law of addition!

``````

<!-- SECTION 2 -->

(demo-condition-roots-python)=
``````{dropdown} Conditioning of polynomial roots
:open: false
::::{grid} 1 1 2 2

The polynomial $p(x) = \frac{1}{3}(x-1)(x-1-\epsilon)$ has roots $1$ and $1+\epsilon$. For small values of $\epsilon$, the roots are ill-conditioned. 

:::{card}
:columns: 5

The statement `x, y = 10, 20` makes individual assignments to both `x` and `y`.

:::
::::


```{code-cell}
ep = 1e-6   
a, b, c = 1/3, (-2-ep)/3, (1+ep)/3   # coefficients of p
```

Here are the roots as computed by the quadratic formula.

```{code-cell}
d = np.sqrt(b**2 - 4*a*c)
r1 = (-b - d) / 2*a
r2 = (-b + d) / 2*a
print(r1, r2)
```

The display of `r2` suggests that the last five digits or so are inaccurate. The relative error in the value is 

```{code-cell}
print(abs(r2 - (1+ep)) / (1+ep))
```

The condition number of the roots is inversely proportional to $2\epsilon$, the difference between them. Thus, roundoff error in the data can grow in the result to be roughly

```{code-cell}
print(np.finfo(float).eps / ep)
```

This matches the observation well.
``````

<!-- SECTION 3 -->

(function-horner-python)=
`````{dropdown} **Horner's algorithm for evaluating a polynomial**
```{code-block} python
:lineno-start: 1
def horner(c,x):
    """
    horner(c,x)

    Evaluate a polynomial whose coefficients are given in descending order in `c`, at the point `x`, using Horner's rule.
    """

    n = len(c)
    y = c[0]
    for k in range(1,n):
        y = x*y + c[k]   

    return y
```
`````

(demo-algorithms-horner-python)=
``````{dropdown} Using a function
:open: false

Here we show how to use `horner` to evaluate a polynomial. First, we have to ensure that the `FNC` package is imported.
  
```{code-cell} ipython3
import FNC
```

Here is the help string for the function:

```{code-cell} ipython3
help(FNC.horner)
```

We now define a vector of the coefficients of $p(x)=(x−1)^3=x^3−3x^2+3x−1$, in descending degree order. Note that the textbook's functions are all in a namespace called `FNC`, to help distinguish them from other Python commands and modules.

```{code-cell} ipython3
c = np.array([1, -3, 3, -1])
print(FNC.horner(c, 1.6))
```

The above is the value of $p(1.6)$, up to a rounding error.

``````

<!-- SECTION 4 -->

(demo-stability-quadbad-python)= 
``````{dropdown} Instability of the quadratic formula
:open: false

```{index} ! Python; scientific notation
```

::::{grid} 1 1 2 2
We apply the quadratic formula to find the roots of a quadratic via {eq}`quadunstable`. 

:::{card}
A number in scientific notation is entered as `1.23e4` rather than as `1.23*10^{4}`.
:::
::::

```{code-cell}
a = 1;  b = -(1e6 + 1e-6);  c = 1;
x1 = (-b + np.sqrt(b**2 - 4*a*c)) / 2*a
x2 = (-b - np.sqrt(b**2 - 4*a*c)) / 2*a
print(x1, x2)
```

The first value is correct to all stored digits, but the second has fewer than six accurate digits:

```{code-cell}
error = np.abs(1e-6 - x2) / 1e-6 
print(f"There are {-np.log10(error):.2f} accurate digits.")
```

 The instability is easily explained. Since $a=c=1$, we treat them as exact numbers. First, we compute the condition numbers with respect to $b$ for each elementary step in finding the "good" root:

| Calculation | Result | $\kappa$ |
|:------------|:-------|:---------|
|$u_1 = b^2$  | $1.000000000002000\times 10^{12}$ |  2 |
|$u_2 = u_1 - 4$ | $9.999999999980000\times 10^{11}$  | $\approx 1.00$ |
|$u_3 = \sqrt{u_2}$ | $999999.9999990000$ | 1/2 |
|$u_4 = u_3 - b$ | $2000000$ | $\approx 0.500$ |
|$u_5 = u_4/2$ | $1000000$  | 1 |

Using {eq}`condition-chain`, the chain rule for condition numbers, the conditioning of the entire chain is the product of the individual steps, so there is essentially no growth of relative error here. However, if we use the quadratic formula for the "bad" root, the next-to-last step becomes $u_4=(-u_3) - b$, and now  $\kappa=|u_3|/|u_4|\approx 5\times 10^{11}$. So we can expect to lose 11 digits of accuracy, which is what we observed. The key issue is the subtractive cancellation in this one step.
``````

(demo-stability-quadgood-python)=
````{dropdown} Stable alternative to the quadratic formula
We repeat the rootfinding experiment of {numref}`Demo %s <demo-stability-quadbad>` with an alternative algorithm.

```{code-cell}
a = 1;  b = -(1e6 + 1e-6);  c = 1;
```

First, we find the "good" root using the quadratic formula.

```{code-cell}
x1 = (-b + np.sqrt(b**2 - 4*a*c)) / 2*a
```

Then we use the identity $x_1x_2=\frac{c}{a}$ to compute the smaller root:

```{code-cell}
x2 = c / (a * x1)
print(x1, x2)
```

To be sure we have an accurate result, we compute its relative error.

```{code-cell}
print(abs(x2 - 1e-6) / 1e-6)
```
````

(demo-stability-roots-python)=
````{dropdown} Backward error
:open: false


Our first step is to construct a polynomial with six known roots.

```{code-cell} ipython3
r = [-2, -1, 1, 1, 3, 6]
p = np.poly(r)
print(p)
```

Now we use a standard numerical method for finding those roots, pretending that we don't know them already. This corresponds to $\tilde{y}$ in {numref}`Definition {number} <definition-stability-backward>`.

```{code-cell} ipython3
r_computed = np.sort(np.roots(p))
print(r_computed)
```

Here are the relative errors in each of the computed roots.

```{code-cell} ipython3
print(np.abs(r - r_computed) / r)
```

It seems that the forward error is acceptably close to machine epsilon for double precision in all cases except the double root at $x=1$. This is not a surprise, though, given the poor conditioning at such roots.

Let's consider the backward error. The data in the rootfinding problem is the polynomial coefficients. We can apply poly to find the coefficients of the polynomial (that is, the data) whose roots were actually computed by the numerical algorithm. This corresponds to $\tilde{x}$ in {numref}`Definition {number} <definition-stability-backward>`. 

```{code-cell} ipython3
p_computed = np.poly(r_computed)
print(p_computed)
```

We find that in a relative sense, these coefficients are very close to those of the original, exact polynomial:

```{code-cell} ipython3
print(np.abs(p - p_computed) / p)
```

In summary, even though there are some computed roots relatively far from their correct values, they are nevertheless the roots of a polynomial that is very close to the original.
