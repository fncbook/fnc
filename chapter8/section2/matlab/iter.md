---
kernelspec:
  display_name: MATLAB
  language: matlab
  name: jupyter_matlab_kernel
numbering:
  headings: false
---
```{code-cell}
:tags: [remove-cell]
cd  /Users/driscoll/Documents/GitHub/fnc/matlab
FNC_init
```
[**Demo %s**](#demo-power-iter)

We will experiment with the power iteration on a 5×5 matrix with prescribed eigenvalues and dominant eigenvalue at 1.

```{code-cell}
ev = [1, -0.75, 0.6, -0.4, 0];
% Make a triangular matrix with eigenvalues on the diagonal.
A = triu(ones(5, 5), 1) + diag(ev);
```

We run the power iteration 60 times. The best estimate of the dominant eigenvalue is the last entry of the first output.

```{code-cell}
[beta, x] = poweriter(A, 60);
format long
beta(1:12)
```

We check for linear convergence using a log-linear plot of the error.

```{code-cell}
err = 1 - beta;
clf,  semilogy(abs(err), '.-')
title('Convergence of power iteration')
xlabel('k'),  ylabel(('|\lambda_1 - \beta_k|'));
```

The asymptotic trend seems to be a straight line, consistent with linear convergence. To estimate the convergence rate, we look at the ratio of two consecutive errors in the linear part of the convergence curve. The ratio of the first two eigenvalues should match the observed rate.

```{code-cell}
theory = ev(2) / ev(1)
observed = err(40) / err(39)
```

Note that the error is supposed to change sign on each iteration. The effect of these alternating signs is that estimates oscillate around the exact value.

```{code-cell}
beta(26:29)
```

In practical situations, we don't know the exact eigenvalue that the algorithm is supposed to find. In that case we would base errors on the final $\beta$ that was found, as in the following plot.

```{code-cell}
err = beta(end) - beta(1:end-1);
semilogy(abs(err), '.-')
title('Convergence of power iteration')
xlabel('k'),  ylabel(('|\beta_{60} - \beta_k|'));
```

The results are very similar until the last few iterations, when the limited accuracy of the reference value begins to show. That is, while it is a good estimate of $\lambda_1$, it is less good as an estimate of the error in nearby estimates.
