---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
---
```{code-cell}
:tags: [remove-cell]
using FundamentalsNumericalComputation
FNC.init_format()
```

(demo-float-accuracy-julia)=
:::::{prf:example}
Recall the grade-school approximation to the number $\pi$.

```{code-cell}
@show p = 22/7;
```
::::{grid} 1 1 2 2
:::{grid-item}
:columns: 5

Not all the digits displayed for `p` are the same as those of $\pi$. 

:::
:::{grid-item-card}
:columns: 7

The value of `pi` is predefined and equivalent to `π`, which is entered by typing `\pi` followed immediately by the <kbd>Tab</kbd> key.

:::
::::

```{code-cell}
@show float(π);
```

```{index} ! Julia; string interpolation
```

::::{grid} 1 1 2 2

:::{grid-item}
:columns: 5

The absolute and relative accuracies of the approximation are as follows.

:::
:::{grid-item-card}
:columns: 7

A dollar sign `$` in a string substitutes (or *interpolates*) the named variable or expression into the string.

:::
::::

```{code-cell}
acc = abs(p-π)
println("absolute accuracy = $acc")
println("relative accuracy = $(acc/π)")
```

::::{grid} 1 1 2 2

:::{grid-item}
:columns: 5

Here we calculate the number of accurate digits in `p`.

:::
:::{grid-item-card}
:columns: 7

The `log` function is for the natural log. For other common bases, use `log10` or `log2`.

:::
::::

```{code-cell}
println("Number of accurate digits = $(-log10(acc/π))")
```
This last value could be rounded down by using `floor`.

:::::

(demo-float-julia)=
:::::{prf:example}
In Julia, `1` and `1.0` are different values, because they have different types:

```{code-cell}
@show typeof(1);
@show typeof(1.0);
```

The standard choice for floating-point values is `Float64`, which is double precision using 64 binary bits. We can see all the bits by using `bitstring`.

```{code-cell}
bitstring(1.0)
```

::::{grid} 1 1 2 2

:::{grid-item}
:columns: 7

The first bit determines the sign of the number:

:::
:::{grid-item-card}
:columns: 5

Square brackets concatenate the contained values into vectors.

:::
::::

```{code-cell}
[bitstring(1.0), bitstring(-1.0)]
```

The next 11 bits determine the exponent (scaling) of the number, and so on.

```{code-cell}
[bitstring(1.0), bitstring(2.0)]
```

The sign bit, exponent, and significand in {eq}`floatpoint` are all directly accessible.

```{code-cell}
x = 1.0
@show sign(x),exponent(x),significand(x);
```

```{code-cell}
x = 0.125
@show sign(x),exponent(x),significand(x);
```

::::{grid} 1 1 2 2

:::{grid-item}
:columns: 6

The spacing between floating-point values in $[2^n,2^{n+1})$ is $2^n \epsilon_\text{mach}$, where $\epsilon_\text{mach}$ is machine epsilon. You can get its value from the `eps` function in Julia. By default, it returns the value for double precision.

:::
:::{grid-item-card}
:columns: 6

To call a function, including `eps`, you must use parentheses notation, even when there are no input arguments.

:::
::::

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
@show floatmin(),floatmax();
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
@show 5.0,Int(5.0);
```

If you try to convert a noninteger floating-point value into an integer you get an `InexactValue` error. This occurs whenever you try to force a type conversion that doesn't make clear sense.
:::::

(demo-float-arithmetic-julia)=
:::::{prf:example}

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
:::::

(demo-condition-roots-julia)=
:::::{prf:example}
::::{grid} 1 1 2 2

:::{grid-item}
:columns: 7

The polynomial $p(x) = \frac{1}{3}(x-1)(x-1-\epsilon)$ has roots $1$ and $1+\epsilon$. For small values of $\epsilon$, the roots are ill-conditioned. 

:::
:::{card}
:columns: 5

The statement `x,y = 10,20` makes individual assignments to both `x` and `y`.

:::
::::


```{code-cell}
ϵ = 1e-6   # type \epsilon and then press Tab
a,b,c = 1/3,(-2-ϵ)/3,(1+ϵ)/3   # coefficients of p
```

Here are the roots as computed by the quadratic formula.

```{code-cell}
d = sqrt(b^2-4a*c)
r₁ = (-b-d)/(2a)   # type r\_1 and then press Tab
r₂ = (-b+d)/(2a)
(r₁,r₂)
```

The display of `r₂` suggests that the last five digits or so are inaccurate. The relative error in the value is 

```{code-cell}
abs(r₂ - (1+ϵ)) / (1+ϵ)
```

The condition number of the roots is inversely proportional to $2\epsilon$, the difference between them. Thus roundoff error in the data can grow in the result to be roughly

```{code-cell}
eps()/ϵ
```

This matches the observation well.
:::::

(function-horner-julia)=
````{prf:function} horner

**Horner's algorithm for evaluating a polynomial**

```{code-block} julia
:lineno-start: 1
"""
    horner(c,x)

Evaluate a polynomial whose coefficients are given in ascending
order in `c`, at the point `x`, using Horner's rule.
"""
function horner(c,x)
    n = length(c)
    y = c[n]
    for k in n-1:-1:1
        y = x*y + c[k]
    end
    return y
end
```
````

(demo-algorithms-horner-julia)=
:::::{prf:example}
Here we show how to use {numref}`Function {number} <function-horner>` to evaluate a polynomial. It's not a part of core Julia, so you need to download and install this text's package once, and load it for each new Julia session. The download is done by the following lines.

```{code-cell}
:tags: [remove-output]
import Pkg
Pkg.add(Pkg.PackageSpec(url="https://github.com/fncbook/FundamentalsNumericalComputation"));
```

```{index} ! Julia; using
```

::::{grid} 1 1 2 2
:gutter: 2

:::{grid-item}
:columns: 5


Once installed, any package can be loaded with the `using` command, as follows.


:::
:::{grid-item-card}
:columns: 7


Many Julia functions, including the ones in this text, are in packages that must be loaded via `using` or `import` in each session. Sometimes a `using` statement can take a few seconds or even minutes to execute, if packages have been installed or updated. 

:::
::::

```{code-cell}
:tags: [remove-output]
using FundamentalsNumericalComputation
```

::::{grid} 1 1 2 2
:gutter: 2

:::{grid-item}
:columns: 7



For convenience, this package also imports many other packages used throughout the book and makes them available as though you had run a `using` command for each of them. 


:::
:::{grid-item-card}
:columns: 5



If you are not sure where a particular function is defined, you can run `methods` on the function name to find all its definitions.

:::
::::

Returning to `horner`, let us define a vector of the coefficients of $p(x)=(x-1)^3=x^3-3x^2+3x-1$, in ascending degree order.

```{code-cell}
c = [-1,3,-3,1]
```

```{index} ! Julia; FNC, ! Julia; namespace
```

::::{grid} 1 1 2 2
:gutter: 2

:::{grid-item}
:columns: 7



In order to avoid clashes between similarly named functions, Julia has boxed all the book functions into a **namespace** called `FNC`. We use this namespace whenever we invoke one of the functions.


:::
:::{grid-item-card}
:columns: 5
 


You must use the module name when a package is loaded by `import`, but when loaded via `using`, some functions may be available with no prefix.

:::
::::

```{code-cell}
FNC.horner(c,1.6)
```

The above is the value of $p(1.6)$.

While the namespace does lead to a little extra typing, a nice side effect of using this paradigm is that if you type `FNC.` (including the period) and hit the <kbd>Tab</kbd> key, you will see a list of all the functions known in that namespace.

The multi-line string at the start of {numref}`Function {number} <function-horner>` is documentation, which we can access using `?FNC.horner`.
:::::