---
kernelspec:
  display_name: Python 3
  language: python
  name: python3
numbering:
  headings: false
---

# Chapter 1 

Python implementations

## Functions

(function-horner-python)=
`````{dropdown} **Horner's algorithm for evaluating a polynomial**
:open: true
````{literalinclude} fncbook/fncbook/chapter01.py
:lineno-start: 1
:start-at: def horner
:end-at: return
:filename: horner.py
:language: python
````
`````

## Examples

```{code-cell} ipython3
:tags: remove-cell
exec(open("FNC_init.py").read())
```

## 1.1 @section-intro-floating-point
(demo-float-accuracy-python)=
``````{dropdown} @demo-float-accuracy
```{tip} Getting started with Python
:class: dropdown
:open:
See @section-setup-python for guidance on how to set up Python for the demos in this book.
```

Recall the grade-school approximation to the number $\pi$.

```{code-cell} ipython3
p = 22/7
print(p)
```
Not all the digits displayed for `p` are the same as those of $\pi$. 
```{tip}
:class: dropdown
The value of `pi` is predefined in the `numpy` package.
```

```{code-cell} ipython3
print(pi)
```

The absolute and relative accuracies of the approximation are as follows:
```{tip}
:class: dropdown
We often use [Python f-strings](https://docs.python.org/3/tutorial/inputoutput.html#tut-f-strings) to format numerical output. 
```


```{code-cell} ipython3
print(f"absolute accuracy: {abs(p - pi)}")
```

```{code-cell} ipython3
rel_acc = abs(p - pi) / pi
print("relative accuracy: {rel_acc:.4e}")
```

Here we calculate the number of accurate digits in `p`:
```{tip}
:class: dropdown
The `log` function is for the natural log. For other common bases, use `log10` or `log2`.
```


```{code-cell} ipython3
print(f"accurate digits: {-log10(rel_acc):.1f}")
```
``````

(demo-float-double-python)=
``````{dropdown} @demo-float-double
:open: false

Python has native `int` and `float` types.

```{code-cell} ipython3
print(f"The type of {1} is {type(1)}")
print(f"The type of {float(1)} is {type(1.0)}")
```

The `numpy` package has its own `float` types:

```{code-cell} ipython3
one = float64(1)
print(f"The type of {one} is {type(one)}")
```

Both `float` and `float64` are double precision, using 64 binary bits per value. Although it is not normally necessary to do so, we can deconstruct a float into its significand and exponent:

```{code-cell} ipython3
x = 3.14
mantissa, exponent = frexp(x)
print(f"significand: {mantissa * 2}, exponent: {exponent - 1}")
```

```{code-cell} ipython3
mantissa, exponent = frexp(x / 8)
print(f"significand: {mantissa * 2}, exponent: {exponent - 1}")
```

The spacing between floating-point values in $[2^n,2^{n+1})$ is $2^n \epsilon_\text{mach}$, where $\epsilon_\text{mach}$ is machine epsilon, given here for double precision:

```{code-cell} ipython3
mach_eps = finfo(float).eps
print(f"machine epsilon is {mach_eps:.4e}")
```

Because double precision allocates 52 bits to the significand, the default value of machine epsilon is $2^{-52}$.

```{code-cell} ipython3
print(f"machine epsilon is 2 to the power {log2(mach_eps)}")
```

A common mistake is to think that $\epsilon_\text{mach}$ is the smallest floating-point number. It's only the smallest *relative to 1*. The correct perspective is that the scaling of values is limited by the exponent, not the significand. The actual range of positive values in double precision is

```{code-cell}
finf = finfo(float)
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
``````{dropdown} @demo-float-arithmetic
:open: false

There is no double precision number between $1$ and $1+\varepsilon_\text{mach}$. Thus, the following difference is zero despite its appearance.

```{code-cell} ipython3
eps = finfo(float).eps
e = eps/2
print((1.0 + e) - 1.0)
```

However, $1-\varepsilon_\text{mach}/2$ is a double precision number, so it and its negative are represented exactly:

```{code-cell} ipython3
print(1.0 + (e - 1.0))
```

This is now the "correct" result. But we have found a rather shocking breakdown of the associative law of addition!

``````

### 1.2 @section-intro-conditioning

(demo-condition-roots-python)=
``````{dropdown} @demo-condition-roots
:open: false

The polynomial $p(x) = \frac{1}{3}(x-1)(x-1-\epsilon)$ has roots $1$ and $1+\epsilon$. For small values of $\epsilon$, the roots are ill-conditioned. 

```{tip}
:class: dropdown
The statement `x, y = 10, 20` makes individual assignments to both `x` and `y`.
```


```{code-cell}
ep = 1e-6   
a, b, c = 1/3, (-2 - ep) / 3, (1 + ep) / 3   # coefficients of p
```

Here are the roots as computed by the quadratic formula.

```{code-cell}
d = sqrt(b**2 - 4*a*c)
r1 = (-b - d) / (2*a)
r2 = (-b + d) / (2*a)
print(r1, r2)
```

The display of `r2` suggests that the last five digits or so are inaccurate. The relative error in the value is 

```{code-cell}
print(abs(r1 - 1) / abs(1))
print(abs(r2 - (1 + ep)) / abs(1 + ep))
```

The condition number of each root is 
$$
\kappa(r_i) = \frac{|r_i|}{|r_1-r_2|} \approx \frac{1}{\epsilon}. 
$$
Thus, relative error in the data at the level of roundoff can grow in the result to be roughly
```{code-cell}
print(finfo(float).eps / ep)
```

This matches the observation pretty well.
``````

### 1.3 @section-intro-algorithms

(demo-algorithms-horner-python)=
``````{dropdown} @demo-algorithms-horner
:open: false

Here we show how to use `horner` to evaluate a polynomial. First, we have to ensure that the book's package is imported.
  
```{code-cell} ipython3
import fncbook as FNC
```

Here is the help string for the function:

```{code-cell} ipython3
help(FNC.horner)
```

We now define a vector of the coefficients of $p(x)=(x−1)^3=x^3−3x^2+3x−1$, in descending degree order. Note that the textbook's functions are all in a namespace called `FNC`, to help distinguish them from other Python commands and modules.

```{code-cell} ipython3
c = array([1, -3, 3, -1])
print(FNC.horner(c, 1.6))
```

The above is the value of $p(1.6)$, up to a rounding error.

``````

### 1.4 @section-intro-stability

(demo-stability-quadbad-python)= 
``````{dropdown} @demo-stability-quadbad
:open: false

```{index} ! Python; scientific notation
```

We apply the quadratic formula to find the roots of a quadratic via {eq}`quadunstable`. 
```{tip}
:class: dropdown
A number in scientific notation is entered as `1.23e4` rather than as `1.23*10^{4}`.
```

```{code-cell}
a = 1;  b = -(1e6 + 1e-6);  c = 1;
x1 = (-b + sqrt(b**2 - 4*a*c)) / 2*a
x2 = (-b - sqrt(b**2 - 4*a*c)) / 2*a
print(x1, x2)
```

The first value is correct to all stored digits, but the second has fewer than six accurate digits:

```{code-cell}
error = abs(1e-6 - x2) / 1e-6 
print(f"There are {-log10(error):.2f} accurate digits.")
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
````{dropdown} @demo-stability-quadgood
We repeat the rootfinding experiment of {numref}`Demo %s <demo-stability-quadbad>` with an alternative algorithm.

```{code-cell}
a = 1;  b = -(1e6 + 1e-6);  c = 1;
```

First, we find the "good" root using the quadratic formula.

```{code-cell}
x1 = (-b + sqrt(b**2 - 4*a*c)) / 2*a
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
````{dropdown} @demo-stability-roots
:open: false


Our first step is to construct a polynomial with six known roots.

```{code-cell} ipython3
r = [-2, -1, 1, 1, 3, 6]
p = poly(r)
print(p)
```

Now we use a standard numerical method for finding those roots, pretending that we don't know them already. This corresponds to $\tilde{y}$ in {numref}`Definition {number} <definition-stability-backward>`.

```{code-cell} ipython3
r_computed = sort(roots(p))
print(r_computed)
```

Here are the relative errors in each of the computed roots.

```{code-cell} ipython3
print(abs(r - r_computed) / r)
```

It seems that the forward error is acceptably close to machine epsilon for double precision in all cases except the double root at $x=1$. This is not a surprise, though, given the poor conditioning at such roots.

Let's consider the backward error. The data in the rootfinding problem is the polynomial coefficients. We can apply poly to find the coefficients of the polynomial (that is, the data) whose roots were actually computed by the numerical algorithm. This corresponds to $\tilde{x}$ in {numref}`Definition {number} <definition-stability-backward>`. 

```{code-cell} ipython3
p_computed = poly(r_computed)
print(p_computed)
```

We find that in a relative sense, these coefficients are very close to those of the original, exact polynomial:

```{code-cell} ipython3
print(abs(p - p_computed) / p)
```

In summary, even though there are some computed roots relatively far from their correct values, they are nevertheless the roots of a polynomial that is very close to the original.