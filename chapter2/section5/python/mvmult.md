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
[**Demo %s**](#demo-flops-mvmult)


Here is a straightforward implementation of matrix-vector multiplication.

```{code-cell} 
n = 6
A = random.rand(n, n)
x = ones(n)
y = zeros(n)
for i in range(n):
    for j in range(n):
        y[i] += A[i, j] * x[j]   # 2 flops
```

Each of the loops implies a summation of flops. The total flop count for this algorithm is

$$
\sum_{i=1}^n \sum_{j=1}^n 2 = \sum_{i=1}^n 2n = 2n^2.
$$

Since the matrix $\mathbf{A}$ has $n^2$ elements, all of which have to be involved in the product, it seems unlikely that we could get a flop count that is smaller than $O(n^2)$ in general.

Let's run an experiment with the built-in matrix-vector multiplication. We assume that flops dominate the computation time and thus measure elapsed time. 

```{code-cell} 
N = 400 * arange(1, 11)
t = []
print("  n           t")
for i, n in enumerate(N):
    A = random.randn(n, n)  
    x = random.randn(n)
    start = timer()
    for j in range(50): A @ x
    t.append(timer() - start)
    print(f"{n:5}   {t[-1]:10.3e}")
```

The reason for doing multiple repetitions at each value of $n$ above is to avoid having times so short that the resolution of the timer is a factor.

Looking at the timings just for $n=2000$ and $n=4000$, they have ratio:

```{code-cell} 
print(t[9] / t[4])
```

If the run time is dominated by flops, then we expect this ratio to be 

$$
\frac{2(4000)^2}{2(2000)^2}=4.
$$
