---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
Here we show how to use {numref}`Function {number} <function-horner>` to evaluate a polynomial. It's not a part of core Julia, so you need to download and install this text's package once, and load it for each new Julia session. The download is done by the following lines.

```{code-cell}
:tags: [remove-output]
#import Pkg
#Pkg.add("FNCBook");
```

```{index} ! Julia; using
```
Once installed, any package can be loaded with the `using` command, as follows.

```{tip}
Many Julia functions, including the ones in this text, are in packages that must be loaded via `using` or `import` in each session. Sometimes a `using` statement can take a few seconds or even minutes to execute, if packages have been installed or updated. 
```

```{code-cell}
#using FundamentalsNumericalComputation
using FNCFunctions
FNC = FNCFunctions
```


For convenience, this package also imports many other packages used throughout the book and makes them available as though you had run a `using` command for each of them. 


:::
```{tip}
If you are not sure where a particular function is defined, you can run `methods` on the function name to find all its definitions.
```

Returning to `horner`, let us define a vector of the coefficients of $p(x)=(x-1)^3=x^3-3x^2+3x-1$, in ascending degree order.

```{code-cell}
c = [-1, 3, -3, 1]
```

```{index} ! Julia; FNC, ! Julia; namespace
```


In order to avoid clashes between similarly named functions, Julia has boxed all the book functions into a **namespace** called `FNC`. We use this namespace whenever we invoke one of the functions.


:::
```{tip}
You must use the module name when a package is loaded by `import`, but when loaded via `using`, some functions may be available with no prefix.
```

```{code-cell}
FNC.horner(c, 1.6)
```

The above is the value of $p(1.6)$.

While the namespace does lead to a little extra typing, a nice side effect of using this paradigm is that if you type `FNC.` (including the period) and hit the <kbd>Tab</kbd> key, you will see a list of all the functions known in that namespace.

The multi-line string at the start of {numref}`Function {number} <function-horner>` is documentation, which we can access using `?FNC.horner`.
