---
jupytext:
  cell_metadata_filter: -all
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: '0.8'
    jupytext_version: 1.5.0
kernelspec:
  display_name: Julia faststart
  language: julia
  name: julia-fast
---

# Horner's rule

Here we show how to use the function `horner` to evaluate a polynomial. It's not a part of core Julia, so you need to download and install this text's package once, and load it for each new Julia session. The download is done with the following lines--but the second one is commented out here to avoid dumping a lot of irrelevant output into this notebook.

```{code-cell}
import Pkg
#Pkg.add(Pkg.PackageSpec(url="https://github.com/fncbook/fnc.git"));
```

Once installed, any package (including a standard library that comes with Julia) can be loaded with the `using` command, as follows:[^compile]

[^compile]: Sometimes a `using` or `import` statement can take a few seconds or even minutes to execute. This book depends on a large package that must be compiled on first use and after every update.

```{code-cell}
using FundamentalsNumericalComputation
```

For your convenience, this package also loads and imports all the packages used throughout the book.

Now let us define a vector of the coefficients of $p(x)=(x-1)^3=x^3-3x^2+3x-1$, in ascending degree order.

```{code-cell}
c = [-1,3,-3,1]
```

In order to avoid clashes between similarly named functions, Julia has sandboxed all the book functions into a **namespace** called `FNC`. We use this namespace whenever we invoke one of the functions.

```{code-cell}
FNC.horner(c,1.6)
```

The above is the value of $p(1.6)$.

While the namespace does lead to a little extra typing, a nice side effect of using this paradigm is that if you type `FNC.` (including the period) and hit the Tab key, you will see a list of all the functions known in that namespace.
