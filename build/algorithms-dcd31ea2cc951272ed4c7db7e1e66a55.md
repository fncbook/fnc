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
\begin{split}
p(x) &= c_1 + c_2 x + \cdots + c_n x^{n-1} \\
&= \Bigl( \bigl( (c_n x  + c_{n-1}) x + c_{n-2} \bigr) x + \cdots +c_{2} \Bigr)x + c_{1}.
\end{split}
```

(function-horner)=
``````{prf:algorithm} horner
`````{tab-set} 
````{tab-item} Julia
:sync: julia
:::{embed} #function-horner-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
matlab
```` 

````{tab-item} Python
:sync: python
:::{embed} #function-horner-python
:::
```` 
`````
``````

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
``````{prf:example}
`````{tab-set} 
````{tab-item} Julia
:sync: julia
:::{embed} #demo-algorithms-horner-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
matlab
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-algorithms-horner-python
:::
```` 
`````
``````


The quoted lines at the beginning of {numref}`Function {number} <function-horner>` are a documentation string. The function itself starts off with the keyword `function`, followed by a list of its input arguments. The first of these is presumed to be a vector, whose length can be obtained and whose individual components are accessed through square bracket notation. After the computation is finished, the `return` keyword indicates which value or values are to be returned to the caller.

The `Polynomials` package for Julia provides its own fast methods for polynomial evaluation that supersede our simple {numref}`Function {number} <function-horner>` function. This will be the case for all the codes in this book because the problems we study are well-known and important. In a more practical setting, you would take implementations of basic methods for granted and build on top of them.

## Writing your own functions


Any collection of statements organized around solving a type of problem should probably be wrapped in a function. One clue is that if you find yourself copying and pasting code, perhaps with small changes in each instance, you should probably be writing a function instead.

`````{tab-set} 
````{tab-item} Julia
:sync: julia

```{index} ! Julia; functions
```
Functions can be defined in text files with the extension `.jl`, at the command line (called the *REPL prompt*), or in notebooks. 

As seen in {numref}`Function {number} <function-horner>`, one way to start a function definition is with the `function` keyword, followed by the function name and the input arguments in parentheses. For example, to represent the mathematical function $e^{\sin x}$, we could use

``` julia
function myfun(x)
    s = sin(x)
    return exp(s)
end
```

```{index} ! Julia; return
```

The `return` statement is used to end execution of the function and return one or more (comma-separated) values to the caller of the function. If an executing function reaches its `end` statement without encountering a `return` statement, then it returns the result of the most recent statement, but this is considered poor style.

For a function with a short definition like the one above, there is a more compact syntax to do the same thing:

``` julia
myfun(x) = exp(sin(x))
```

```{index} ! Julia; anonymous functions
```
You can also define **anonymous functions** or **lambda functions**, which are typically simple functions that are provided as inputs to other functions. This is done with an arrow notation. For example, to plot the function above (in the `Plots` package) without permanently creating it, you could enter

``` julia
plot(x -> exp(sin(x)), 0, 6)
```

As in most languages, input arguments and variables defined within a function have scope limited to the function itself. However, they can access values defined within an enclosing scope. For instance:

``` julia
mycfun(x) = exp(c*sin(x))
c = 1;  mycfun(3)   # returns exp(1*sin(3))
c = 2;  mycfun(3)   # returns exp(2*sin(3))
```

There's a lot more to be said about functions in Julia, but this is enough to get started.

```` 

````{tab-item} MATLAB
:sync: matlab

```{index} ! MATLAB; functions
```
Functions can be defined in text files with the extension `.m`, at the command line, or in Live Scripts. 

As seen in {numref}`Function {number} <function-horner>`, one way to start a function definition is with the `function` keyword, followed one or more output arguments in brackets, the function name, and the input arguments in parentheses. For example, to represent the mathematical function $e^{\sin x}$, we could use the following in a file called `myfun.m`:

``` matlab
function [y] = myfun(x)
    s = sin(x);
    y = exp(s);
end
```

Whatever value is assigned to `y` when the function terminates will be returned as the output of the function. 

For a function with a short definition like the one above, there is a more compact syntax to do the same thing:

``` matlab
myfun = @(x) exp(sin(x));
```

```{index} ! MATLAB; anonymous functions
```

The syntax on the right of the `=` above defines an **anonymous function** (called a **lambda function** in computer science), which can be used in place without giving it a name as we did here. We'll have examples of doing this later on. 

As in most languages, input arguments and variables defined within a function have scope limited to the function itself. However, they can access values defined within an enclosing scope, with those values being locked in at the time of creation. For instance:

``` matlab
c = 1;
mycfun = @(x) exp(c*sin(x));
mycfun(3)   % returns exp(1*sin(3))
c = 2;  
mycfun(3)   % also returns exp(1*sin(3))
mycfun = @(x) exp(c*sin(x));   % redefines mycfun
mycfun(3)   % now returns exp(2*sin(3))

```

There's a lot more to be said about functions, but this is enough to get started.

```` 

````{tab-item} Python
:sync: python

Functions can be defined in text files with the extension `.py`, at the command line (called the *REPL prompt*), or in notebooks. 

As seen in {numref}`Function {number} <function-horner>`, one way to start a function definition is with the `def` keyword, followed by the function name and the input arguments in parentheses, ending with a colon. The statements for the body of the function must then all be indented. For example, to represent the mathematical function $e^{\sin x}$, we could use

``` python
def myfun(x):
    s = np.sin(x)
    return np.exp(s)
```

```{index} ! Python; return
```

The `return` statement is used to end execution of the function and return one or more (comma-separated) values to the caller of the function. 

```{tip}
If an executing function reaches its `end` statement without encountering a `return` statement, then the output is undefined, which is a common source of bugs.
```

For a function with a short definition like the one above, there is a more compact syntax to do the same thing:

``` python
myfun = lambda x : np.exp(np.sin(x))
```

The syntax on the right of the `=` above defines an **anonymous function** (called a **lambda function** in computer science), which can be used in place without giving it a name as we did here. We'll have examples of doing this later on. 

As in most languages, input arguments and variables defined within a function have scope limited to the function itself. However, they can access values defined within an enclosing scope. For instance:

``` python
mycfun = lambda x : np.exp(c * np.sin(x))
c = 1;  print(mycfun(3))   # exp(1*sin(3))
c = 2;  print(mycfun(3))   # exp(2*sin(3))
```

There's a lot more to be said about functions in Python, but this is enough to get started.

```` 
``````


## Exercises

1. ⌨ Write a function `poly1(p)` that returns the value of a polynomial $p(x) = c_1 + c_2 x + \cdots + c_n x^{n-1}$ at $x=-1$. You should do this directly, not by a call to or imitation of {numref}`Function {number} <function-horner>`. Test your function on $r(x)=3x^3-x+1$ and $s(x)=2x^2-x$.

    (problem-algorithms-samplevar)=
2. ⌨  In statistics, one defines the variance of sample values $x_1,\ldots,x_n$ by
  
    ```{math}
    :label: samplevar
    s^2 = \frac{1}{n-1} \sum_{i=1}^n (x_i - \overline{x})^2,
    \qquad \overline{x} = \frac{1}{n} \sum_{i=1}^n x_i.
    ```

    Write a function `samplevar(x)` that takes as input a vector `x` of any length and returns $s^2$ as calculated by the formula. You should test your function on the vectors `ones(100)` and `rand(200)`. If you enter `using Statistics` in Julia, then you can compare to the results of the `var` function.

3. ⌨  Let `x` and `y` be vectors whose entries give the coordinates of the $n$ vertices of a polygon, given in counterclockwise order. Write a function `polygonarea(x,y)` that computes and returns the area of the polygon using this formula based on Green's theorem:
  
    ```{math}
    A = \frac{1}{2} \left| \sum_{k=1}^n x_k y_{k+1} - x_{k+1}y_k \right|.
    ```

    Here $n$ is the number of polygon vertices, and it's understood that $x_{n+1}=x_1$ and $y_{n+1}=y_1$. (Note: The function `abs` computes absolute value.) Test your functions on a square and an equilateral triangle.
