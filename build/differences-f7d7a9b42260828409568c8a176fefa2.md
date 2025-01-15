---
numbering: false
---

(section-key-differences)=
# Key language differences 

This book makes extensive use of computing to both illustrate and demonstrate important facts. The three languages used are Julia, MATLAB, and Python (specifically, NumPy and SciPy). These languages share a great deal of syntax and functionality because they address the same needs. As the youngest, Julia adopted some conventions from both of the others, while NumPy itself was strongly influenced by MATLAB. 

However, there are some major differences in the languages that you should be aware of. You probably aren't here for a history lesson, so I won't get into the whys, but please do understand that some different choices were made for good reasons in context.

We focus here on issues that come up often in this book. There is a broader quick reference [here](https://cheatsheets.quantecon.org/) and a much more comprehensive analysis [here](https://docs.julialang.org/en/v1/manual/noteworthy-differences/).

## Indexes

Let's get this out of the way.

```{caution}
In MATLAB and Julia, indexing starts at 1. In Python, it starts at 0.
```

Wars have been fought over less. In math, it's often—but not always—convenient to start at 1. In this book, we try to make an index in code equivalent to its mathematical definition, regardless of the language. 

Also, MATLAB uses parentheses `()` for indexing, while Julia and Python use brackets `[]`.[^brackets]

[^brackets]: Yes, we use the American convention of "brackets" to mean "square brackets." Take that, England!

## Linear algebra

The most consequential difference has to do with how they represent the key objects from linear algebra, matrices and vectors. In MATLAB, every numerical value can be considered a matrix. A 1×1 matrix is also often considered to be a scalar (number), and a 1×n or n×1 matrix is often interpretable as a vector. This flexibility can be convenient, though it can also lead to confusion. Row and column vectors are often interchangeable, but not always. For instance,

```matlab
% MATLAB
x = [1, 2, 3];       % row vector
y = [3; 4; 5; 6];    % column vector
x(2) == x(1, 2)      % true
x(3) == y(1, 1)      % true
x(3, 1) == y(1)      % error
```

In Julia, scalars, vectors, and matrices are different things. In particular, a vector is not identical to a matrix with one row or one column. However, there are contexts in which a vector is interpreted to have an implicit column shape. Moreover, vectors can hold any type of data, not just numbers. These rules have some subtle consequences. For example,

```julia
# Julia
[ [1, 2], 3 ]
```

is a vector of two elements, the first of which is itself a vector. However,

```julia
# Julia
[ [1, 2]; 3 ]
```

is a vector with three numerical elements, because it is the "vertical stacking" of a 2-vector and a scalar.

In Python, scalars, vectors, and matrices are also different things. Moreover, Python has lists that are more fundamental than vectors. Thus, while `[1,2,3]` is a vector in MATLAB and Julia, it is a list in Python. In NumPy, often lists can be used in place of vectors, but not always. When vector shape matters, the implicit NumPy interpretation is as a row.

### Vectors of vectors

One recurring situation that highlights the different approaches is when data takes the form of a sequence of identically-sized vectors. This could be an iteration in a vector space, or the solution of a system of ODEs at selected times, for example. 

MATLAB's choice is to arrange the data as an array. Say that `x1`, `x2`, and `x3` are 2×1 vectors. Then a compact representation is

```matlab
% MATLAB
X = [x1, x2, x3]    % 2×3 
X = [x1 x2 x3]      % same thing 
```

MATLAB would call `X` a matrix, though it might or might not require the mathematical properties of linear algebra. One downside is that you must decide on its shape in advance: are the members of the sequence rows or columns? It's an arbitrary choice, but a consequential one.

In Julia, the two syntaxes above are valid but different:

```julia
# Julia
X = [x1, x2, x3]    # 3-element vector of 2-element vectors
X = [x1 x2 x3]      # 2×3 matrix
```

Each representation has advantages and drawbacks, and you have to know which you are working with. 

In NumPy, all roads lead to the same place:

```python
# Python
X = np.array([x1, x2, x3])    # 3-element vector of 2-element vectors AND a 3×2 matrix
```

There is nothing to decide or remember here, because the two representations are interchangeable. 

## Functions

In MATLAB, functions are traditionally defined in files[^mfiles] with the same name as the function. The function is called from the command line or from other functions by using the root name of the file. So, if you have a function defined in a file `foo.m` on the MATLAB path, `foo` calls it. If you want to refer to the function *without* calling it, so that it can be passed as an argument to another function, you refer to it as `@foo`.

In Julia and Python, functions can be defined anywhere. The name of the function is a reference to it, not a call. To call the function, you use `()` after the name, even if it does not take any arguments.

All three languages also have the idea of an *anonymous function*, which is never given a name. These are equivalent:

Julia | MATLAB | Python
--- | --- | ---
x -> x^2 | @(x) x^2 | lambda x: x**2


[^mfiles]: MATLAB has blurred things in recent years by allowing function definitions within functions and script files, but not at the command line. Because of the way this book is constructed, we have to define functions in separate files.

## Vectorization

Often we want to apply a function to every element of a vector or matrix. This can be called *vectorization* or *broadcasting*, though both terms can have other related meanings. In MATLAB and NumPy, most of the functions provided, such as `cos` and `exp`, are defined to work elementwise on vectors and matrices. NumPy takes it further, defining operators like `*`, `/`' and `**` to work elementwise as well. MATLAB, on the other hand, requires the use of `.*`, `./`, and `.^` for these, reserving `*`, `/`, and `^` for operations defined in the sense of linear algebra.

When it comes to `*`, `/`' and `^`, Julia is like MATLAB. But it is more extreme in a different way. Mathematical functions such as `cos` and `exp` are *not* defined to work elementwise. Instead, you have to use a dot `.` after the function name. For example,

```julia
# Julia
x = [-1, 0, 1] * pi
cos(x)    # error
cos.(x)   # [-1, 1, -1]
```

At first, this seems like an enormous pain in the booty, especially when you want to combine multiple operations. But there is a nice shortcut for that:

```julia
# Julia
cos.(x.^3 + 2x.^2 .- 1)    # 🤮
@. cos(x^3 + 2x^2 - 1)     # 💯
```

And there is a huge reward:[^reward] you can use a vectorizing dot with *any* function, including ones that you define. So you only ever need to write a function as though it will operate on a scalar. 

[^reward]: Julia fanboys would claim that the error you get for `x - 1` rather than `x .- 1` when `x` is a vector is itself a reward, because it can help you catch subtle mistakes. We concede the point, but it's hard to feel appreciative at the moment you receive one of these errors.

By contrast, MATLAB and NumPy make you implement elementwise behavior yourself whenever you write a function that needs it. In MATLAB, we often end up with expressions such as

``` matlab
% MATLAB
@(x) cos(x.^3 + 2*x.^2 - 1)
```

so that the function is available elementwise. NumPy would not require any dotting in this case; for more complex functions, it provides `np.vectorize` to help, but it is not a panacea.
