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
[**Demo %s**](#demo-shooting-mems)

We revisit {numref}`Demo {number} <demo-shooting-naive>` but let {numref}`Function {number} <function-shoot>` do the heavy lifting.

```{code-cell}
lamb = 0.6
phi = lambda r, w, dwdr: lamb / w**2 - dwdr / r
a, b = finfo(float).eps, 1
```

We specify the given and unknown endpoint values.

```{code-cell}
ga = lambda w, dw : dw       # w'=0 at left
gb = lambda w, dw : w - 1    # w=1 at right
```

In this setting, we need to provide initial guesses for $w(a)$ and $w'(a)$.

```{code-cell}
init = array([0.8, 0])
r, w, dw_dx = FNC.shoot(phi, a, b, ga, gb, init)
plot(r, w)
title("Shooting solution")
xlabel("$r$"),  ylabel("$w(r)$");
```

The value of $w$ at $r=1$, meant to be exactly one, was computed to be

```{code-cell}
print(f"w at right end is {w[-1]}")
```

The accuracy is consistent with the error tolerance used for the IVP solution by `shoot`. The initial value $w(0)$ that gave this solution is

```{code-cell}
print(f"w at left end is {w[0]}")
```
