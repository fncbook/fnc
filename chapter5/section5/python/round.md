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
[**Demo %s**](#demo-fdconverge-round)

Let $f(x)=e^{-1.3x}$. We apply finite-difference formulas of first, second, and fourth order to estimate $f'(0)=-1.3$.

```{code-cell}
f = lambda x: exp(-1.3 * x)
exact = -1.3

h_ = array([1 / 10**(n+1) for n in range(12)])
FD = zeros((len(h_), 3))
for (i, h) in enumerate(h_):
    nodes = h * linspace(-2, 2, 5)
    vals = f(nodes)
    FD[i, 0] = dot(array([0, 0, -1, 1, 0]) / h, vals)
    FD[i, 1] = dot(array([0, -1/2, 0, 1/2, 0]) / h, vals)
    FD[i, 2] = dot(array([1/12, -2/3, 0, 2/3, -1/12]) / h, vals)

results = PrettyTable()
results.add_column("h", h_)
results.add_column("FD1", FD[:, 0])
results.add_column("FD2", FD[:, 1])
results.add_column("FD4", FD[:, 2])
print(results)
```

They all seem to be converging to $-1.3$. The convergence plot reveals some interesting structure to the errors, though.

```{code-cell}
loglog(h_, abs(FD[:, 0] + 1.3), "-o", label="FD1")
loglog(h_, abs(FD[:, 1] + 1.3), "-o", label="FD2")
loglog(h_, abs(FD[:, 2] + 1.3), "-o", label="FD4")
gca().invert_xaxis()
plot(h_, 0.1 * 2 ** (-52) / h_, "--", color="k", label="$O(h^{-1})$")
xlabel("$h$")
ylabel("total error")
title("FD error with roundoff")
legend();
```

Again the graph is made so that $h$ decreases from left to right. The errors are dominated at first by truncation error, which decreases most rapidly for the fourth-order formula. However, increasing roundoff error eventually equals and then dominates the truncation error as $h$ continues to decrease. As the order of accuracy increases, the crossover point moves to the left (greater efficiency) and down (greater accuracy).
