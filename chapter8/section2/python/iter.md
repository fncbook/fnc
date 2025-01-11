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
[**Demo %s**](#demo-power-iter)

We will experiment with the power iteration on a 5×5 matrix with prescribed eigenvalues and dominant eigenvalue at 1.

```{code-cell}
ev = [1, -0.75, 0.6, -0.4, 0]
A = triu(ones([5, 5]), 1) + diag(ev)    # triangular matrix, eigs on diagonal
```

We run the power iteration 60 times. The first output should be a sequence of estimates converging to the dominant eigenvalue—which, in this case, we set up to be 1.

```{code-cell}
beta, x = FNC.poweriter(A, 60)
print(beta)
```

We check for linear convergence using a log-linear plot of the error.

```{code-cell}
err = 1 - beta
semilogy(arange(60), abs(err), "-o")
ylim(1e-10, 1)
xlabel("$k$")
ylabel("$|\\lambda_1 - \\beta_k|$")
title("Convergence of power iteration");
```

The asymptotic trend seems to be a straight line, consistent with linear convergence. To estimate the convergence rate, we look at the ratio of two consecutive errors in the linear part of the convergence curve. The ratio of the first two eigenvalues should match the observed rate.

```{code-cell}
print(f"theory: {ev[1] / ev[0]:.5f}")
print(f"observed: {err[40] / err[39]:.5f}")
```

Note that the error is supposed to change sign on each iteration. The effect of these alternating signs is that estimates oscillate around the exact value.

```{code-cell}
print(beta[26:30])
```

In practical situations, we don't know the exact eigenvalue that the algorithm is supposed to find. In that case we would base errors on the final $\beta$ that was found, as in the following plot.

```{code-cell}
err = beta[-1] - beta
semilogy(arange(60), abs(err), "-o")
ylim(1e-10, 1), xlabel("$k$")
ylabel("$|\\lambda_1 - \\beta_k|$")
title("Convergence of power iteration");
```

The results are very similar until the last few iterations, when the limited accuracy of the reference value begins to show. That is, while it is a good estimate of $\lambda_1$, it is less good as an estimate of the error in nearby estimates.

