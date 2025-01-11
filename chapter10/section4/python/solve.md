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
[**Demo %s**](#demo-linear-solve)


```{code-cell}
exact = lambda x: exp( sin(x) )
```

The problem is presented above in our standard form, so we can identify the coefficient functions in the ODE. Each should be coded as a function.

```{code-cell}
p = lambda x: -cos(x)
q = sin
r = lambda x: 0 * x    # must be a function 
```

We solve the BVP and compare the result to the exact solution.

```{code-cell}
x, u = FNC.bvplin(p, q, r, [0, pi/2], 1, exp(1), 25)
```

```{code-cell}
subplot(2, 1, 1)
plot(x, u)
ylabel("solution"),  title("Solution of the BVP")

subplot(2, 1, 2)
plot(x, exact(x) - u, "-o")
ylabel("error");
```
