---
kernelspec:
  display_name: Python 3
  language: python
  name: python3
numbering:
  headings: false
---
```{code-cell}
:tags: [remove-cell]
exec(open("../../../python/FNC_init.py").read())
```
[**Demo %s**](#demo-algorithms-horner)


Here we show how to use `horner` to evaluate a polynomial. First, we have to ensure that the book's package is imported.
  
```{code-cell} ipython3
import fncbook as FNC
```

Here is the help string for the function:

```{code-cell} ipython3
help(FNC.horner)
```

We now define a vector of the coefficients of $p(x)=(x−1)^3=x^3−3x^2+3x−1$, in descending degree order. Note that the textbook's functions are all in a namespace called `FNC`, to help distinguish them from other Python commands and modules.

```{code-cell} ipython3
c = array([1, -3, 3, -1])
print(FNC.horner(c, 1.6))
```

The above is the value of $p(1.6)$, up to a rounding error.

