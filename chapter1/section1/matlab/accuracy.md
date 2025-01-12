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
[**Demo %s**](#demo-float-accuracy)

:::{tip} Getting started in MATLAB
:class: dropdown
See @section-setup-matlab for instructions on how to install functions for MATLAB for this book.
:::

Recall the grade-school approximation to the number $\pi$.

```{index} MATLAB; format
```

```{tip}
:class: dropdown
The number of digits displayed is controlled by `format`, but the underlying values are not affected by it.
```

```{code-cell}
format long
p = 22/7
```
Not all the digits displayed for `p` are the same as those of $\pi$. 

```{tip}
:class: dropdown
The value of `pi` is predefined.
```

The absolute and relative accuracies of the approximation are as follows.

```{code-cell}
abs_accuracy = abs(p - pi)
rel_accuracy = abs(p - pi) / pi
```

Here we calculate the number of accurate digits in `p`.
```{tip}
:class: dropdown
The `log` function is for the natural log. For other common bases, use `log10` or `log2`.
```

```{code-cell}
format short
accurate_digits = -log10(rel_accuracy)
```

