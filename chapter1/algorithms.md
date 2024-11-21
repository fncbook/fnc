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

(section-intro-algorithms)=
# Algorithms

```{index} ! algorithm
```

An idealized mathematical problem $f(x)$ can usually only be approximated using a finite number of steps in finite precision. A complete set of instructions for transforming data into a result is called an **algorithm**. In most cases it is reasonable to represent an algorithm by another mathematical function, denoted here by $\tilde{f}(x)$. 

Even simple problems can be associated with multiple algorithms.

```{index} ! Horner's algorithm
```

````{prf:example}
Suppose we want to find an algorithm that maps a given $x$ to the value of the polynomial $f(x)= 5x^3 + 4x^2 + 3x + 2$. Representing $x^2$ as $(x)(x)$, we can find it with one multiplication. We can then find $x^3=(x)(x^2)$ with one more multiplication. We can then apply all the coefficients (three more multiplications) and add all the terms (three additions), for a total of 8 arithmetic operations.

There is a more efficient algorithm, however: organize the polynomial according to **Horner's algorithm**,

```{math}
f(x) = 2 + x \bigl( 3 + x( 4 + 5x) \bigr).
```

In this form you can see that evaluation takes only 3 additions and 3 multiplications. The savings represent 25% of the original computational effort, which could be significant if repeated billions of times.
````

## Algorithms as code

Descriptions of algorithms may be presented as a mixture of mathematics, words, and computer-style instructions called *pseudocode*, which varies in syntax and level of formality. In this book we use pseudocode to explain the outline of an algorithm, but the specifics are usually presented as working code. 

Of all the desirable traits of code, we emphasize clarity the most after correctness. We do not represent our programs as always being the shortest, fastest, or most elegant. Our primary goal is to illustrate and complement the mathematical underpinnings, while occasionally pointing out key implementation details.

As our first example, {numref}`Function {number} <function-horner>` implements an algorithm that applies Horner's algorithm to a general polynomial, using the identity

```{math}
:label: horner
p(x) &= c_1 + c_2 x + \cdots + c_n x^{n-1} \\
&= \Bigl( \bigl( (c_n x  + c_{n-1}) x + c_{n-2} \bigr) x + \cdots +c_{2} \Bigr)x + c_{1}.
```

(function-horner)=

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

```{index} Julia; length, Julia; indexing arrays
```

:::{admonition} About the code
The `length` function in line 1 returns the number of elements in vector `c`. The syntax `c[n]` accesses element `n` of a vector `c`. In Julia, the first index of a vector is 1 by default, so in line 2, the last element of `c` is accessed.

The `for` / `end` construct is a *loop*. The local variable `k` is assigned the value `n-1`, then the loop body is executed, then `k` is assigned `n-2`, the body is executed again, and so on until finally `k` is set to 1 and the body is executed for the last time.

The `return` statement in line 13 terminates the function and specifies one or more values to be returned to the caller. A function may have more than one `return` statement, in which case the first one encountered terminates the function; however, that coding style is mostly discouraged.
:::

```{index} ! Julia; for
```

(demo-algorithms-horner)=
```{prf:example}
```





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





The quoted lines at the beginning of {numref}`Function {number} <function-horner>` are a documentation string. The function itself starts off with the keyword `function`, followed by a list of its input arguments. The first of these is presumed to be a vector, whose length can be obtained and whose individual components are accessed through square bracket notation. After the computation is finished, the `return` keyword indicates which value or values are to be returned to the caller.

The `Polynomials` package for Julia provides its own fast methods for polynomial evaluation that supersede our simple {numref}`Function {number} <function-horner>` function. This will be the case for all the codes in this book because the problems we study are well-known and important. In a more practical setting, you would take implementations of basic methods for granted and build on top of them.

## Writing your own functions

```{index} ! Julia; functions
```

Functions are a primary way of working in Julia. Any collection of statements organized around solving a type of problem should probably be wrapped in a function. Functions can be defined in text files with the extension `.jl`, at the command line (called the *REPL prompt*), or in notebooks. 

As seen in {numref}`Function {number} <function-horner>`, one way to start a function definition is with the `function` keyword, followed by the function name and the input arguments in parentheses. For example, to represent the mathematical function $e^{\sin x}$, we could use

``` julia
function myfun(x)
    s = sin(x)
    return exp(s)
end
```

```{index} ! Julia; return
```

The `return` statement is used to end execution of the function and return one or more (comma-separated) values to the caller of the function. If an executing function reaches its `end` statement without encountering a `return` statement, then it returns the result of the most recent statement.

For a function with a short definition like the one above, there is a more compact syntax to do the same thing:

``` julia
myfun(x) = exp(sin(x))
```

You can also define **anonymous functions** or **lambda functions**, which are typically simple functions that are provided as inputs to other functions. This is done with an arrow notation. For example, to plot the function above (in the `Plots` package) without permanently creating it, you could enter

``` julia
plot( x->exp(sin(x)), 0, 6 )
```

As in most languages, input arguments and variables defined within a function have scope limited to the function itself. However, they can access values defined within an enclosing scope. For instance:

``` julia
mycfun(x) = exp(c*sin(x))
c = 1;  mycfun(3)   # returns exp(1*sin(3))
c = 2;  mycfun(3)   # returns exp(2*sin(3))
```

There's a lot more to be said about functions in Julia, but this is enough to get started.

## Exercises

1. ⌨ Write a Julia function

    ``` julia
    function poly1(p)
    ```

    that evaluates a polynomial $p(x) = c_1 + c_2 x + \cdots + c_n x^{n-1}$ at $x=-1$. You should do this directly, not by a call to or imitation of {numref}`Function {number} <function-horner>`. Test your function on $r(x)=3x^3-x+1$ and $s(x)=2x^2-x$.

    (problem-algorithms-samplevar)=
2. ⌨  In statistics, one defines the variance of sample values $x_1,\ldots,x_n$ by
  
    ```{math}
        :label: samplevar
        s^2 = \frac{1}{n-1} \sum_{i=1}^n (x_i - \overline{x})^2,
      \qquad \overline{x} = \frac{1}{n} \sum_{i=1}^n x_i.
    ```

    Write a Julia function

    ``` julia
    function samplevar(x)
    ```

    that takes as input a vector `x` of any length and returns $s^2$ as calculated by the formula. You should test your function on the vectors `ones(100)` and `rand(200)`. If you enter `using Statistics` in Julia, then you can compare to the results of the `var` function.

3. ⌨  Let `x` and `y` be vectors whose entries give the coordinates of the $n$ vertices of a polygon, given in counterclockwise order. Write a function

    ``` julia
    function polygonarea(x,y)
    ```

    that computes the area of the polygon, using this formula based on Green's theorem:
  
    ```{math}
    A = \frac{1}{2} \left| \sum_{k=1}^n x_k y_{k+1} - x_{k+1}y_k \right|.
    ```

    Here $n$ is the number of polygon vertices, and it's understood that $x_{n+1}=x_1$ and $y_{n+1}=y_1$. (Note: The function `abs` computes absolute value.) Test your functions on a square and an equilateral triangle.
