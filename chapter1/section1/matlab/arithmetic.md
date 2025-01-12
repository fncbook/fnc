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
[**Demo %s**](#demo-float-arithmetic)


There is no double-precision number between $1$ and $1+\epsilon_\text{mach}$. Thus the following difference is zero despite its appearance.

```{code-cell}
( 1 + eps / 2 ) - 1
```

However, the spacing between floats in $[1/2,1)$ is $\macheps/2$, so both $1-\macheps/2$ and its negative are represented exactly:

```{code-cell}
1 - (1 - eps / 2)
```

This is now the expected result. But we have found a rather shocking breakdown of the associative law of addition!
