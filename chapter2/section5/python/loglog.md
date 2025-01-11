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
[**Demo %s**](#demo-flops-loglog)

Let's repeat the experiment of the previous example for more, and larger, values of $n$.

```{code-cell} 
N = arange(400, 6200, 200)
t = zeros(len(N))
for i, n in enumerate(N):
    A = random.randn(n,n)  
    x = random.randn(n)
    start = timer()
    for j in range(20): A@x
    t[i] = timer() - start
```

Plotting the time as a function of $n$ on log-log scales is equivalent to plotting the logs of the variables, but is formatted more neatly. 

```{code-cell} 
fig, ax = subplots()
ax.loglog(N, t, "-o", label="observed")
ylabel("elapsed time (sec)");
xlabel("$n$");
title("Timing of matrix-vector multiplications");
```

You can see that while the full story is complicated, the graph is trending to a straight line of positive slope. For comparison, we can plot a line that represents $O(n^2)$ growth exactly. (All such lines have slope equal to 2.)

```{code-cell} 
ax.loglog(N, t[-1] * (N/N[-1])**2, "--", label="$O(n^2)$")
ax.legend();  fig
```
