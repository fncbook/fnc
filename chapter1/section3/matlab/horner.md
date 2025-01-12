---
kernelspec:
  display_name: MATLAB
  language: matlab
  name: jupyter_matlab_kernel
numbering:
  headings: false
---
```{code-cell}
:tags: [remove-cell]
cd  /Users/driscoll/Documents/GitHub/fnc/matlab
FNC_init
```
[**Demo %s**](#demo-algorithms-horner)

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
