---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.12
numbering:
  headings: false
---
# Chapter 1

## Functions

(function-horner-julia)=
`````{dropdown} **Horner's algorithm for evaluating a polynomial**
:open:
```{literalinclude} FNCFunctions/src/chapter01.jl
:filename: horner.jl
:start-after: # begin horner
:end-before: # end horner
:linenos: true
:language: julia
```

```{index} ! Julia; length, ! Julia; for
```

:::{admonition} About the code
:class: dropdown
The quoted lines at the beginning are a documentation string. The function itself starts off with the keyword `function`, followed by a list of its input arguments. The first of these is presumed to be a vector, whose length can be obtained and whose individual components are accessed through square bracket notation. After the computation is finished, the `return` keyword indicates which value or values are to be returned to the caller.

The `length` function in line 8 returns the number of elements in vector `c`. Here, that value is one greater than the degree of the polynomial. The syntax `c[i]` accesses element `i` of a vector `c`. In Julia, the first index of a vector is 1 by default, so in line 9, the last element of `c` is accessed.

The `for` / `end` construct in lines 10–12 is a *loop*. The local variable `k` is assigned the value `n-1`, then the loop body is executed, then `k` is assigned `n-2`, the body is executed again, and so on until finally `k` is set to 1 and the body is executed for the last time.

The `return` statement in line 13 terminates the function and specifies one or more values to be returned to the caller.
:::

:::{important}
The `Polynomials` package for Julia provides its own fast methods for polynomial evaluation that supersede our simple `horner`. This will be the case for all the codes in this book because the problems we study are well-known and important. In a more practical setting, you would take implementations of basic methods for granted and build on top of them.
:::

`````

```{code-cell}
:tags: remove-cell
include("FNC_init.jl")
1+1;
```

## Examples

### 1.1 @section-intro-floating-point
(demo-float-accuracy-julia)=
``````{dropdown} @demo-float-accuracy
:open:
:::{tip} Getting started in Julia
:class: dropdown
See @section-setup-julia for instructions on how to install and use Julia for this book.
:::

Recall the grade-school approximation to the number $\pi$.

```{code-cell}
@show p = 22/7;
```
Not all the digits displayed for `p` are the same as those of $\pi$. 

```{tip}
:class: dropdown
The value of `pi` is predefined and equivalent to `π`, which is entered by typing `\pi` followed immediately by the <kbd>Tab</kbd> key.
```

```{code-cell}
@show float(π);
```

```{index} ! Julia; string interpolation
```


The absolute and relative accuracies of the approximation are as follows.

```{tip}
:class: dropdown
A dollar sign `$` in a string substitutes (or *interpolates*) the named variable or expression into the string.
```

```{code-cell}
acc = abs(p-π)
println("absolute accuracy = $acc")
println("relative accuracy = $(acc/π)")
```

Here we calculate the number of accurate digits in `p`.

```{tip}
:class: dropdown
The `log` function is for the natural log. For other common bases, use `log10` or `log2`.
```

```{code-cell}
println("Number of accurate digits = $(-log10(acc/π))")
```
This last value could be rounded down by using `floor`.

``````

(demo-float-double-julia)=
``````{dropdown} @demo-float-double
:open: 
In Julia, `1` and `1.0` are different values, because they have different types:

```{code-cell}
@show typeof(1);
@show typeof(1.0);
```

The standard choice for floating-point values is `Float64`, which is double precision using 64 binary bits. We can see all the bits by using `bitstring`.

```{code-cell}
bitstring(1.0)
```


The first bit determines the sign of the number:

```{tip}
:class: dropdown
Square brackets concatenate the contained values into vectors.
```

```{code-cell}
[bitstring(1.0), bitstring(-1.0)]
```

The next 11 bits determine the exponent (scaling) of the number, and so on.

```{code-cell}
[bitstring(1.0), bitstring(2.0)]
```

The sign bit, exponent, and significand in {eq}`floatpoint` are all directly accessible.

```{code-cell}
x = 3.14
@show sign(x), exponent(x), significand(x);
```

```{code-cell}
x = x / 8
@show sign(x), exponent(x), significand(x);
```

The spacing between floating-point values in $[2^n,2^{n+1})$ is $2^n \epsilon_\text{mach}$, where $\epsilon_\text{mach}$ is machine epsilon. You can get its value from the `eps` function in Julia. By default, it returns the value for double precision.

```{tip}
:class: dropdown
To call a function, including `eps`, you must use parentheses notation, even when there are no input arguments.
```

```{code-cell} julia
eps()
```

Because double precision allocates 52 bits to the significand, the default value of machine epsilon is $2^{-52}$.

```{code-cell}
log2(eps())
```

The spacing between adjacent floating-point values is proportional to the magnitude of the value itself. This is how relative precision is kept roughly constant throughout the range of values. You can get the adjusted spacing by calling `eps` with a value.

```{code-cell}
eps(1.618)
```

```{code-cell}
eps(161.8)
```

```{code-cell}
nextfloat(161.8)
```

```{code-cell}
ans - 161.8
```

A common mistake is to think that $\epsilon_\text{mach}$ is the smallest floating-point number. It's only the smallest *relative to 1*. The correct perspective is that the scaling of values is limited by the exponent, not the mantissa. The actual range of positive values in double precision is

```{code-cell}
@show floatmin(), floatmax();
```

For the most part you can mix integers and floating-point values and get what you expect.

```{code-cell}
1/7
```

```{code-cell}
37.3 + 1
```

```{code-cell}
2^(-4)
```

There are some exceptions. A floating-point value can't be used as an index into an array, for example, even if it is numerically equal to an integer. In such cases you use `Int` to convert it.

```{code-cell}
@show 5.0, Int(5.0);
```

If you try to convert a noninteger floating-point value into an integer you get an `InexactValue` error. This occurs whenever you try to force a type conversion that doesn't make clear sense.
``````

(demo-float-arithmetic-julia)=
``````{dropdown} @demo-float-arithmetic
:open: 

There is no double-precision number between $1$ and $1+\epsilon_\text{mach}$. Thus the following difference is zero despite its appearance.

```{code-cell}
e = eps()/2
(1.0 + e) - 1.0
```

However, the spacing between floats in $[1/2,1)$ is $\macheps/2$, so both $1-\macheps/2$ and its negative are represented exactly:

```{code-cell}
1.0 + (e - 1.0)
```

This is now the expected result. But we have found a rather shocking breakdown of the associative law of addition!
``````

### 1.2 @section-intro-conditioning
(demo-condition-roots-julia)=
``````{dropdown} @demo-condition-roots
:open: 
The polynomial $p(x) = \frac{1}{3}(x-1)(x-1-\epsilon)$ has roots $1$ and $1+\epsilon$. For small values of $\epsilon$, the roots are ill-conditioned. 

```{tip}
:class: dropdown
The statement `x,y = 10,20` makes individual assignments to both `x` and `y`.
```


```{code-cell}
ϵ = 1e-6   # type \epsilon and then press Tab
a,b,c = 1/3,(-2-ϵ)/3,(1+ϵ)/3   # coefficients of p
```

Here are the roots as computed by the quadratic formula.

```{code-cell}
d = sqrt(b^2 - 4a*c)
r₁ = (-b - d) / (2a)   # type r\_1 and then press Tab
r₂ = (-b + d) / (2a)
(r₁, r₂)
```

The relative errors in these values are 

```{code-cell}
@show abs(r₁ - 1) / abs(1);
@show abs(r₂ - (1+ϵ)) / abs(1+ϵ);
```

The condition number of each root is 
$$
\kappa(r_i) = \frac{|r_i|}{|r_1-r_2|} \approx \frac{1}{\epsilon}. 
$$
Thus, relative error in the data at the level of roundoff can grow in the result to be roughly


```{code-cell}
eps() / ϵ
```

This matches the observation pretty well.
``````

### 1.3 @section-intro-algorithms

(demo-algorithms-horner-julia)=
``````{dropdown} @demo-algorithms-horner
:open: 
Here we show how to use {numref}`Function {number} <function-horner>` to evaluate a polynomial. It's not a part of core Julia, so you need to download and install this text's package once, and load it for each new Julia session. The download is done by the following lines.

```{code-cell}
:tags: remove-output
#import Pkg
#Pkg.add("FNCBook");
```

```{index} ! Julia; using
```
Once installed, any package can be loaded with the `using` command, as follows.

```{tip}
:class: dropdown
Many Julia functions, including the ones in this text, are in packages that must be loaded via `using` or `import` in each session. Sometimes a `using` statement can take a few seconds or even minutes to execute, if packages have been installed or updated. 
```

```{code-cell}
:tags: remove-output
using FNCFunctions
```

```{tip}
:class: dropdown
If you are not sure where a particular function is defined, you can run `methods` on the function name to find all its definitions.
```

Returning to `horner`, let us define a vector of the coefficients of $p(x)=(x-1)^3=x^3-3x^2+3x-1$, in ascending degree order.

```{code-cell}
c = [-1, 3, -3, 1]
```

```{index} ! Julia; FNC, ! Julia; namespace
```

In order to avoid clashes between similarly named functions, Julia has boxed all the book functions into a **namespace** called `FNC`. We must use this namespace whenever we invoke one of the functions.

```{tip}
:class: dropdown
You must use the module name when a package is loaded by `import`, but when loaded via `using`, some functions may be available with no prefix.
```

```{code-cell}
FNC.horner(c, 1.6)
```

The above is the value of $p(1.6)$.

While the namespace does lead to a little extra typing, a nice side effect of using this paradigm is that if you type `FNC.` (including the period) and hit the <kbd>Tab</kbd> key, you will see a list of all the functions known in that namespace.

The multi-line string at the start of {numref}`Function {number} <function-horner>` is documentation, which we can access using `?FNC.horner`.
``````

### 1.4 @section-intro-stability
(demo-stability-quadbad-julia)= 
``````{dropdown} @demo-stability-quadbad
:open:

```{index} ! Julia; scientific notation
```

We apply the quadratic formula to find the roots of a quadratic via {eq}`quadunstable`. 

```{tip}
:class: dropdown
A number in scientific notation is entered as `1.23e4` rather than as `1.23*10^{4}`.
```

```{code-cell}
a = 1;  b = -(1e6 + 1e-6);  c = 1;
@show x₁ = (-b + sqrt(b^2 - 4a*c)) / 2a;
@show x₂ = (-b - sqrt(b^2 - 4a*c)) / 2a;
```

The first value is correct to all stored digits, but the second has fewer than six accurate digits:

```{code-cell}
error = abs(1e-6 - x₂) / 1e-6 
@show accurate_digits = -log10(error);
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

(demo-stability-quadgood-julia)=
````{dropdown} @demo-stability-quadgood
:open:
We repeat the rootfinding experiment of @demo-stability-quadbad with an alternative algorithm.

```{code-cell}
a = 1;  b = -(1e6 + 1e-6);  c = 1;
```

First, we find the "good" root using the quadratic formula.

```{code-cell}
@show x₁ = (-b + sqrt(b^2 - 4a*c)) / 2a;
```

Then we use the identity $x_1x_2=\frac{c}{a}$ to compute the smaller root:

```{code-cell}
@show x₂ = c / (a * x₁);
```

As you see in this output, Julia often suppresses trailing zeros in a decimal expansion. To be sure we have an accurate result, we compute its relative error.

```{code-cell}
abs(x₂ - 1e-6) / 1e-6
```
````

(demo-stability-roots-julia)=
````{dropdown} @demo-stability-roots
:open: 

For this example we will use the `Polynomials` package, which is installed by the `FNC` package.  

```{tip}
:class: dropdown
In the rest of the book, we do not show the `using` statement needed to load the book's package, but you will need to enter it if you want to run the codes yourself.
```

Our first step is to construct a polynomial with six known roots.

```{code-cell}
using Polynomials
r = [-2.0, -1, 1, 1, 3, 6]
p = fromroots(r)
```

Now we use a standard numerical method for finding those roots, pretending that we don't know them already. This corresponds to $\tilde{y}$ in {numref}`Definition {number} <definition-backward-error>`. 

```{code-cell}
r̃ = sort(roots(p))   # type r\tilde and then press Tab
```

```{index} ! Julia; @., Julia; broadcasting
```

Here are the relative errors in each of the computed roots. 

```{tip}
:class: dropdown
The `@.` notation at the start means to do the given operations on each element of the given vectors.
```

```{code-cell}
println("Root errors:") 
@. abs(r - r̃) / r
```

It seems that the forward error is acceptably close to machine epsilon for double precision in all cases except the double root at $x=1$. This is not a surprise, though, given the poor conditioning at such roots.

Let's consider the backward error. The data in the rootfinding problem is the polynomial coefficients. We can apply `fromroots` to find the coefficients of the polynomial (that is, the data) whose roots were actually computed by the numerical algorithm. This corresponds to $\tilde{x}$ in {numref}`Definition {number} <definition-backward-error>`. 

```{code-cell}
p̃ = fromroots(r̃)
```

We find that in a relative sense, these coefficients are very close to those of the original, exact polynomial:

```{code-cell}
c,c̃ = coeffs(p), coeffs(p̃)
println("Coefficient errors:") 
@. abs(c - c̃) / c
```

In summary, even though there are some computed roots relatively far from their correct values, they are nevertheless the roots of a polynomial that is very close to the original.
````