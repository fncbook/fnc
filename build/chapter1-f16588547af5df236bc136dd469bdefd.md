---
kernelspec:
  display_name: MATLAB
  language: matlab
  name: jupyter_matlab_kernel
numbering:
  headings: false
---

# Chapter 1

## Functions 

(function-horner-matlab)=
`````{dropdown} Horner's algorithm for evaluating a polynomial
:open: true
```{literalinclude} FNC-matlab/horner.m
:lineno-start: 1
:language: matlab
```
`````

## Examples

```{code-cell}
:tags: [remove-cell]
cd  /Users/driscoll/Documents/GitHub/fnc/matlab
FNC_init;
pwd;
```

### 1.1 @section-intro-floating-point
(demo-float-accuracy-matlab)=
``````{dropdown} @demo-float-accuracy
:open:

:::{tip} Getting started in MATLAB
:class: dropdown
See @section-setup-matlab for instructions on how to install functions for MATLAB for this book.
:::

Recall the grade-school approximation to the number $\pi$.

```{index} MATLAB; format
```

```{tip}
:class: dropdown
The number of digits displayed is controlled by `format`, but the underlying values are not affected by it.
```

```{code-cell}
format long
p = 22/7
```
Not all the digits displayed for `p` are the same as those of $\pi$. 

```{tip}
:class: dropdown
The value of `pi` is predefined.
```

The absolute and relative accuracies of the approximation are as follows.

```{code-cell}
abs_accuracy = abs(p - pi)
rel_accuracy = abs(p - pi) / pi
```

Here we calculate the number of accurate digits in `p`.
```{tip}
:class: dropdown
The `log` function is for the natural log. For other common bases, use `log10` or `log2`.
```

```{code-cell}
format short
accurate_digits = -log10(rel_accuracy)
```

``````

(demo-float-double-matlab)=
``````{dropdown} @demo-float-double
:open:
In MATLAB, values are double-precision floats unless declared otherwise. 

```{code-cell}
fprintf('1 has type: %s', class(1))
fprintf('1.0 has type: %s', class(1.0))
```

The spacing between floating-point values in $[2^n,2^{n+1})$ is $2^n \epsilon_\text{mach}$, where $\epsilon_\text{mach}$ is machine epsilon. Its value is predefined as `eps`.

```{tip}
:class: dropdown
While you can assign a different value to `eps`, doing so does not change any arithmetic. It's generally a bad idea. 
```

```{code-cell} 
eps
```

Because double precision allocates 52 bits to the significand, the default value of machine epsilon is $2^{-52}$.

```{code-cell}
log2(eps)
```

The spacing between adjacent floating-point values is proportional to the magnitude of the value itself. This is how relative precision is kept roughly constant throughout the range of values. You can get the adjusted spacing by calling `eps` with a value.

```{code-cell}
eps(1.618)
```

```{code-cell}
eps(161.8)
```

```{code-cell}
x = 161.8 + 0.1*eps(161.8);
x - 161.8
```

A common mistake is to think that $\epsilon_\text{mach}$ is the smallest floating-point number. It's only the smallest *relative to 1*. The correct perspective is that the scaling of values is limited by the exponent, not the mantissa. The actual range of positive values in double precision is

```{code-cell}
format short e
[realmin, realmax]
```

``````

(demo-float-arithmetic-matlab)=
``````{dropdown} @demo-float-arithmetic
:open:

There is no double-precision number between $1$ and $1+\epsilon_\text{mach}$. Thus the following difference is zero despite its appearance.

```{code-cell}
( 1 + eps / 2 ) - 1
```

However, the spacing between floats in $[1/2,1)$ is $\macheps/2$, so both $1-\macheps/2$ and its negative are represented exactly:

```{code-cell}
1 - (1 - eps / 2)
```

This is now the expected result. But we have found a rather shocking breakdown of the associative law of addition!
``````

### 1.2 @section-intro-conditioning

(demo-condition-roots-matlab)=
``````{dropdown} @demo-condition-roots
:open: 

The polynomial $p(x) = \frac{1}{3}(x-1)(x-1-\epsilon)$ has roots $1$ and $1+\epsilon$. For small values of $\epsilon$, the roots are ill-conditioned. 


```{code-cell} matlab
ep = 1e-6;
a = 1/3;             % coefficients of p...
b = (-2 - ep) / 3;   % ...
c = (1 + ep) / 3;    % ...in ascending order
```

Here are the roots as computed by the quadratic formula.

```{code-cell}
d = sqrt(b^2 - 4*a*c);
format long   % show all digits
r1 = (-b - d) / (2*a)
r2 = (-b + d) / (2*a)
```

The display of `r2` suggests that the last five digits or so are inaccurate. The relative errors are
```{tip}
:class: dropdown
Putting values inside square brackets creates a vector.
```

```{code-cell}
format short e
err = abs(r1 - 1) ./ abs(1)
err = abs(r2 - (1 + ep)) ./ abs(1 + ep)
```

The condition number of each root is 
$$
\kappa(r_i) = \frac{|r_i|}{|r_1-r_2|} \approx \frac{1}{\epsilon}. 
$$
Thus, relative error in the data at the level of roundoff can grow in the result to be roughly

```{code-cell}
eps / ep
```

This matches the observation pretty well.
``````

### 1.3 @section-intro-algorithms

(demo-algorithms-horner-matlab)=
``````{dropdown} @demo-algorithms-horner
:open: 
Here we show how to use {numref}`Function {number} <function-horner>` to evaluate a polynomial. Let us define a vector of the coefficients of $p(x)=(x-1)^3=x^3-3x^2+3x-1$, in ascending degree order.

```{code-cell}
c = [1, -3, 3, 1]
```

Now we evaluate $p(1.6)$ using the function `horner`.

```{code-cell}
horner(c, 1.6)
```

The result above is the value of $p(1.6)$.

```{tip}
The comments at the start of {numref}`Function {number} <function-horner>` are documentation, which we can access using `help horner`.
```
``````

### 1.4 @section-intro-stability

(demo-stability-quadbad-matlab)= 
``````{dropdown} @demo-stability-quadbad
:open: 

```{index} ! MATLAB; scientific notation
```

We apply the quadratic formula to find the roots of a quadratic via {eq}`quadunstable`. 

```{tip}
:class: dropdown
A number in scientific notation is entered as `1.23e4` rather than as `1.23 * 10^4`.
```

```{code-cell}
a = 1;  b = -(1e6 + 1e-6);  c = 1;
x1 = (-b + sqrt(b^2 - 4*a*c)) / (2*a)
x2 = (-b - sqrt(b^2 - 4*a*c)) / (2*a)
```

The first value is correct to all stored digits, but the second has fewer than six accurate digits:

```{code-cell}
error = abs(x2 - 1e-6) / 1e-6 
accurate_digits = -log10(error)
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

(demo-stability-quadgood-matlab)=
``````{dropdown} @demo-stability-quadgood
:open:
We repeat the rootfinding experiment of {numref}`Demo %s <demo-stability-quadbad>` with an alternative algorithm.

```{code-cell}
a = 1;  b = -(1e6 + 1e-6);  c = 1;
```

First, we find the "good" root using the quadratic formula.

```{code-cell}
x1 = (-b + sqrt(b^2 - 4*a*c)) / (2*a);
```

Then we use the identity $x_1x_2=\frac{c}{a}$ to compute the smaller root:

```{code-cell}
x2 = c / (a*x1)
```

This matches the exact root to the displayed digits; to be sure we have an accurate result, we compute its relative error.

```{code-cell}
abs(x2 - 1e-6) / 1e-6
```
``````

(demo-stability-roots-matlab)=
``````{dropdown} @demo-stability-roots
:open: 

Our first step is to construct a polynomial with six known roots.
```{tip}
:class: dropdown
The `'` operator is used for transposition. Here, we want to make `r` a column vector.
```

```{code-cell}
r = [-2 ,-1, 1, 1, 3, 6]';
p = poly(r)
```

Now we use a standard numerical method for finding those roots, pretending that we don't know them already. This corresponds to $\tilde{y}$ in {numref}`Definition {number} <definition-stability-backward>`. 

```{code-cell}
rr = sort(roots(p))   
```

Here are the relative errors in each of the computed roots. 
```{tip}
:class: dropdown
The `./` operator is used for element-wise division.
```

```{code-cell}
disp("Root errors:") 
abs(r - rr) ./ r
```

It seems that the forward error is acceptably close to machine epsilon for double precision in all cases except the double root at $x=1$. This is not a surprise, though, given the poor conditioning at such roots.

Let's consider the backward error. The data in the rootfinding problem is the polynomial coefficients. We can apply `poly` to find the coefficients of the polynomial (that is, the data) whose roots were actually computed by the numerical algorithm. This corresponds to $\tilde{x}$ in {numref}`Definition {number} <definition-stability-backward>`. 

```{code-cell}
pp = poly(rr)
```

We find that in a relative sense, these coefficients are very close to those of the original, exact polynomial:

```{code-cell}
disp("Coefficient errors:") 
abs(p - pp) ./ abs(p)
```

In summary, even though there are some computed roots relatively far from their correct values, they are nevertheless the roots of a polynomial that is very close to the original.
``````