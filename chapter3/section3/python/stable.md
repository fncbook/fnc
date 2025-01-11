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
[**Demo %s**](#demo-qr-stable)

We'll repeat the experiment of {numref}`Demo {number} <demo-normaleqns-instab>`, which exposed instability in the normal equations. 

```{code-cell}
t = linspace(0, 3, 400)
A = array([ [sin(t)**2, cos((1+1e-7)*t)**2, 1] for t in t ])
x = array([1, 2, 1])
b = A @ x
```

The error in the solution by {numref}`Function {number} <function-lsqrfact>` is similar to the bound predicted by the condition number.

```{code-cell}
print(f"observed error: {norm(FNC.lsqrfact(A, b) - x) / norm(x):.3e}")
print(f"conditioning bound: {cond(A) * finfo(float).eps:.3e}")
```
