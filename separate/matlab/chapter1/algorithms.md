---
numbering:
  enumerator: 1.3.%s
kernelspec:
  display_name: MATLAB
  language: matlab
  name: jupyter_matlab_kernel
---
```{code-cell}
:tags: [remove-cell]
clear all
format short
set(0, 'defaultaxesfontsize', 12)
set(0, 'defaultlinelinewidth', 1.5)
set(0, 'defaultFunctionLinelinewidth', 1.5)
set(0, 'defaultscattermarkerfacecolor', 'flat')
gcf;
set(gcf, 'Position', [0 0 600 350])
addpath FNC-matlab
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
\begin{split}
p(x) &= c_1 + c_2 x + \cdots + c_n x^{n-1} \\
&= \Bigl( \bigl( (c_n x  + c_{n-1}) x + c_{n-2} \bigr) x + \cdots +c_{2} \Bigr)x + c_{1}.
\end{split}
```

```{warning}
A polynomial is represented as a vector of coefficients in all three languages covered by this book. However, in Julia they are given in *ascending* degree order, which is most convenient programmatically, while in MATLAB and NumPy they are given in *descending* order, which is the way we usually write them. 
```

``````{prf:algorithm} horner
:label: function-horner

```{literalinclude} ../FNC_matlab/horner.m
:lineno-start: 1
:language: matlab
```
``````

``````{prf:example} Using a function
:label: demo-algorithms-horner

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

## Writing your own functions

Any collection of statements organized around solving a type of problem should probably be wrapped in a function. One clue is that if you find yourself copying and pasting code, perhaps with small changes in each instance, you should probably be writing a function instead.


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


## Exercises

``````{exercise}
:label: problem-algoruthms-poly1

⌨ Write a function `poly1(p)` that returns the value of a polynomial $p(x) = c_1 + c_2 x + \cdots + c_n x^{n-1}$ at $x=-1$. You should do this directly, not by a call to or imitation of {numref}`Function {number} <function-horner>`. Test your function on $r(x)=3x^3-x+1$ and $s(x)=2x^2-x$.
``````

``````{exercise}
:label: problem-algorithms-samplevar
⌨  In statistics, one defines the variance of sample values $x_1,\ldots,x_n$ by

```{math}
:label: samplevar
s^2 = \frac{1}{n-1} \sum_{i=1}^n (x_i - \overline{x})^2,
\qquad \overline{x} = \frac{1}{n} \sum_{i=1}^n x_i.
```

Write a function `samplevar(x)` that takes as input a vector `x` of any length and returns $s^2$ as calculated by the formula. You should test your function on a vector of 100 ones and a vector of 200 normally distributed random numbers. 
``````

``````{exercise}
:label: problem-algoruthms-area

⌨  Let `x` and `y` be vectors whose entries give the coordinates of the $n$ vertices of a polygon, given in counterclockwise order. Write a function `polygonarea(x, y)` that computes and returns the area of the polygon using this formula based on Green's theorem:

```{math}
:numbered: false
A = \frac{1}{2} \left| \sum_{k=1}^n x_k y_{k+1} - x_{k+1}y_k \right|.
```

Here $n$ is the number of polygon vertices, and it's understood that $x_{n+1}=x_1$ and $y_{n+1}=y_1$. (Note: The function `abs` computes absolute value.) Test your functions on a square and an equilateral triangle.
``````
