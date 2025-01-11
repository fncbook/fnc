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
[**Demo %s**](#demo-stability-roots)



Our first step is to construct a polynomial with six known roots.

```{code-cell} ipython3
r = [-2, -1, 1, 1, 3, 6]
p = poly(r)
print(p)
```

Now we use a standard numerical method for finding those roots, pretending that we don't know them already. This corresponds to $\tilde{y}$ in {numref}`Definition {number} <definition-stability-backward>`.

```{code-cell} ipython3
r_computed = sort(roots(p))
print(r_computed)
```

Here are the relative errors in each of the computed roots.

```{code-cell} ipython3
print(abs(r - r_computed) / r)
```

It seems that the forward error is acceptably close to machine epsilon for double precision in all cases except the double root at $x=1$. This is not a surprise, though, given the poor conditioning at such roots.

Let's consider the backward error. The data in the rootfinding problem is the polynomial coefficients. We can apply poly to find the coefficients of the polynomial (that is, the data) whose roots were actually computed by the numerical algorithm. This corresponds to $\tilde{x}$ in {numref}`Definition {number} <definition-stability-backward>`. 

```{code-cell} ipython3
p_computed = poly(r_computed)
print(p_computed)
```

We find that in a relative sense, these coefficients are very close to those of the original, exact polynomial:

```{code-cell} ipython3
print(abs(p - p_computed) / p)
```

In summary, even though there are some computed roots relatively far from their correct values, they are nevertheless the roots of a polynomial that is very close to the original.
