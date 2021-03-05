# Algorithms

When we idealize a problem as the mathematical function $f(x)$, we are explicitly stating that each input (data) has exactly one correct output. When it comes to implementing a computational version of $f$, though, we usually have some choices to make. A complete set of instructions for mapping data to a result is called an {term}`algorithm`. In most cases it is reasonable to represent each algorithm by another function, which for this section we denote $\tilde{f}(x)$.

Even simple problems can be associated with algorithms that have surprisingly different characteristics.

````{prf:example}
  Suppose we want to find an algorithm that maps a given $x$ to the
  value of the polynomial $f(x)= 5x^3 + 4x^2 + 3x + 2$.  Representing $x^2$ as $(x)(x)$, we can find it with one multiplication. We can then find $x^3=(x)(x^2)$ with one more multiplication. We can then apply all the coefficients (three more multiplications) and add all the terms (three
  additions), for a total of 8 arithmetic operations.

  There is a more efficient algorithm, however: organize the polynomial according to
  **Horner's rule**,
  
```{math}
    f(x) = 2 + x \bigl( 3 + x( 4 + 5x) \bigr).
  ```

  In this form you can see that evaluation takes only 3 additions
  and 3 multiplications. The savings represent 25\% of the original computational effort, which could be significant if repeated billions of times.
````

## Algorithms as code

Descriptions of algorithms vary widely. Sometimes they are presented as a mixture of mathematics, words, and computer-style instructions called pseudocode, which varies in syntax, precision, and formality. In this book we present nontrivial algorithms as code.  As a result, we sometimes sacrifice a little readability and brevity, and in some cases we have to address aspects of the algorithm that are peripheral to our mathematical interests. We think these drawbacks are more than outweighed by the lack of ambiguity and the ability to execute the algorithms yourself—not to mention the opportunity to address some issues that are important even if not very mathematical.

Of all the desirable traits of code, we emphasize clarity the most. We do not represent our programs as being the shortest, fastest, or most elegant. Our primary goal is to illustrate and complement the mathematical underpinnings. Hopefully the language and codes are clear enough that if you would rather use a different computing environment to implement them, you will experience few difficulties. We do recommend an *interactive* environment, however, because being able to tweak inputs and pick apart the results is crucial to understanding.

As our first example, [horner](function-horner) implements an algorithm that applies Horner's rule on a general polynomial, through the identity

```{math}
:label: horner
\begin{split}
  p(x) &= c_1 x^{n-1} + c_2 x^{n-2} + \cdots + c_n \\
  &= \Bigl( \bigl( (c_1 x  + c_2) x + c_3 \bigr) x + \cdots +
  c_{n-1} \Bigr)x + c_{n}.
\end{split}
```

(function-horner)=

````{proof:function} horner

**Horner's rule for evaluating a polynomial**

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

```{prf:example} Julia demo
:class: demo
{doc}`demos/algorithms-horner`
```

The quoted lines at the beginning of {ref}`function-horner` are an optional documentation string. The function itself starts off with the keyword `function`, followed by a list of its input arguments. The first of these is presumed to be a vector, whose length can be obtained and whose individual components are accessed through square bracket notation. After the computation is finished, the `return` keyword indicates which value or values are to be returned to the caller.

The `Polynomials` package for Julia provides its own fast methods for polynomial evaluation that supersede our simple [`horner`](function-horner) function. This will often be the case for codes in this book, because the problems we study are classic and important. In a more practical setting you would take implementations of well-known methods for granted and build on top of them.

## Writing your own functions

```{margin}
Julia code files are named with the extension `.jl`.
```

Functions are a primary way of working in Julia. Any collection of statements organized around solving a type of problem should probably be wrapped in a function. Functions can be defined in their own files or at the command line (i.e., REPL prompt). Typically multiple related functions are grouped into a single file with extension `.jl`.

As seen in [horner](function-horner), one way to start a function definition is with the `function` keyword, followed by the function name and the input arguments in parenthesis. For example, to represent the mathematical function $e^{\sin x}$, we could use

``` julia
function myfun(x)
    s = sin(x)
    return exp(s)
end
```

```{margin}
The `return` keyword stops execution of a function and declares value(s) to be returned.
```

The `return` statement is used to end execution of the function and return one or more (comma-separated) values to the caller of the function. A function may have more than one `return` statement. If an executing function reaches its `end` statement without encountering a `return` statement, then it returns the result of the most recent statement.

For a function with a short definition like the one above, there is a more compact syntax to do the same thing:

``` julia
myfun(x) = exp(sin(x))
```

```{margin}
Anonymous functions can be defined using arrow notation.
```

You can also define **anonymous functions**, which are usually functions so short-lived that they do not need a name. This is done with an arrow notation. For example, to plot the function above (in the `Plots` package) without permanently creating it, you could enter

``` julia
using Plots; plot( x->exp(sin(x)), 0, 6 )
```

As in most languages, input arguments and variables defined within a function have scope limited to the function itself. However, they can access values defined within an enclosing scope. For instance:

``` julia
mycfun(x) = exp(c*sin(x))
c = 1;  mycfun(3)   # returns exp(sin(3))
c = 2;  mycfun(3)   # returns exp(2*sin(3))
```

There's a lot more to be said about functions in Julia, but this is enough to get started.

## Exercises

1. ⌨ Write a Julia function

    ``` julia
    function polyadd(p,q)
    ```

    that returns the coefficient vector $r$ for the sum of two polynomials $p$ and $q$, specified as vectors of coefficients in increasing degree. You should *not* assume that `p` and `q` are vectors of the same length.

    (problem-samplevar)=
1. ⌨  In statistics, one defines the variance of sample values $x_1,\ldots,x_n$ by
  
    ```{math}
        :label: samplevar
        s^2 = \frac{1}{n-1} \sum_{i=1}^n (x_i - \overline{x})^2,
      \qquad \overline{x} = \frac{1}{n} \sum_{i=1}^n x_i.
    ```

    Write a Julia function

    ``` julia
    function samplevar(x)
    ```

    that takes as input a vector `x` of any length and returns $s^2$ as calculated by the formula. You should test your function on some artificial data. If you enter `using Statistics` in Julia, then you can compare to results of the `var` function.

    ```{only} solutions
    function Problem_samplevar

    %% Problem 1.3.3
    % Computing the sample variance of a vector x:
    function s2 = samplevar(x)
    % We can use loops, but the following is cleaner.
    xbar = mean(x);
    n = length(x);
    s2 = sum( (x-xbar).^2 ) / (n-1);
    end

    %%
    % Testing:
    x = rand(10,1);
    exact = var(x);
    value = samplevar(x);
    RelativeError = (exact-value)./exact

    end
    ```

    (problem-polyarea)=
1. ⌨  Let `x` and `y` be vectors whose entries give the coordinates of the vertices of a polygon, given in counterclockwise order. Write a function

    ``` julia
    function polyarea(x,y)
    ```

    that computes the area of the polygon, using this formula based on Green's theorem:
  
    ```{math}
    A = \frac{1}{2} \left| \sum_{k=1}^n x_k y_{k+1} - x_{k+1}y_k \right|.
    ```

    Here $n$ is the number of polygon vertices, and by definition, $x_{n+1}=x_1$ and $y_{n+1}=y_1$.  Test your functions on a square and an equilateral triangle.
