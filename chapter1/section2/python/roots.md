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
[**Demo %s**](#demo-condition-roots)


The polynomial $p(x) = \frac{1}{3}(x-1)(x-1-\epsilon)$ has roots $1$ and $1+\epsilon$. For small values of $\epsilon$, the roots are ill-conditioned. 

```{tip}
:class: dropdown
The statement `x, y = 10, 20` makes individual assignments to both `x` and `y`.
```


```{code-cell}
ep = 1e-6   
a, b, c = 1/3, (-2 - ep) / 3, (1 + ep) / 3   # coefficients of p
```

Here are the roots as computed by the quadratic formula.

```{code-cell}
d = sqrt(b**2 - 4*a*c)
r1 = (-b - d) / (2*a)
r2 = (-b + d) / (2*a)
print(r1, r2)
```

The display of `r2` suggests that the last five digits or so are inaccurate. The relative error in the value is 

```{code-cell}
print(abs(r1 - 1) / abs(1))
print(abs(r2 - (1 + ep)) / abs(1 + ep))
```

The condition number of each root is 
$$
\kappa(r_i) = \frac{|r_i|}{|r_1-r_2|} \approx \frac{1}{\epsilon}. 
$$
Thus, relative error in the data at the level of roundoff can grow in the result to be roughly
```{code-cell}
print(finfo(float).eps / ep)
```

This matches the observation pretty well.
