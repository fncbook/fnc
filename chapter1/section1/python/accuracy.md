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
[**Demo %s**](#demo-float-accuracy)

```{tip} Getting started with Python
:class: dropdown
:open:
See @section-setup-python for guidance on how to set up Python for the demos in this book.
```

Recall the grade-school approximation to the number $\pi$.

```{code-cell} ipython3
p = 22/7
print(p)
```
Not all the digits displayed for `p` are the same as those of $\pi$. 
```{tip}
:class: dropdown
The value of `pi` is predefined in the `numpy` package.
```

```{code-cell} ipython3
print(pi)
```

The absolute and relative accuracies of the approximation are as follows:
```{tip}
:class: dropdown
We often use [Python f-strings](https://docs.python.org/3/tutorial/inputoutput.html#tut-f-strings) to format numerical output. 
```


```{code-cell} ipython3
print(f"absolute accuracy: {abs(p - pi)}")
```

```{code-cell} ipython3
rel_acc = abs(p - pi) / pi
print("relative accuracy: {rel_acc:.4e}")
```

Here we calculate the number of accurate digits in `p`:
```{tip}
:class: dropdown
The `log` function is for the natural log. For other common bases, use `log10` or `log2`.
```


```{code-cell} ipython3
print(f"accurate digits: {-log10(rel_acc):.1f}")
```
