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
[**Demo %s**](#demo-flops-lufact)


We'll test the conclusion of $O(n^3)$ flops experimentally using the `lu` function imported from `scipi.linalg`.

```{code-cell}
from scipy.linalg import lu
N = arange(200, 2600, 200)
t = zeros(len(N))
for i, n in enumerate(N):
    A = random.randn(n,n)  
    start = timer()
    for j in range(5): lu(A)
    t[i] = timer() - start
```

We plot the timings on a log-log graph and compare it to $O(n^3)$. The result could vary significantly from machine to machine, but in theory the data should start to parallel the line as $n\to\infty$.

```{code-cell}
loglog(N, t, "-o", label="obseved")
loglog(N, t[-1] * (N / N[-1])**3, "--", label="$O(n^3)$")
legend();
xlabel("$n$");
ylabel("elapsed time (sec)");
title("Timing of LU factorizations");
```
