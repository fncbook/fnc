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
[**Demo %s**](#demo-inviter-conv)

We set up a $5\times 5$ triangular matrix with prescribed eigenvalues on its diagonal.

```{code-cell}
ev = [1, -0.75, 0.6, -0.4, 0];
A = triu(ones(5, 5), 1) + diag(ev);
```

We run inverse iteration with the shift $s=0.7$. The result should converge to the eigenvalue closest to 0.7, which we know to be 0.6 here.

```{code-cell}
s = 0.7;
[beta, x] = inviter(A, s, 30);
format short
beta(1:10)
```

The convergence is again linear.

```{code-cell}
err = abs(0.6 - beta);
semilogy(abs(err),'.-')
title('Convergence of inverse iteration')
xlabel('k'), ylabel(('|\lambda_j - \beta_k|'));
```

Let's reorder the eigenvalues to enforce {eq}`shiftorder`.
```{tip}
:class: dropdown
The second output of `sort` returns the index permutation needed to sort the given vector.
```

```{code-cell}
[~, idx] = sort(abs(ev - s));
ev = ev(idx)
```

Now it is easy to compare the theoretical and observed linear convergence rates.

```{code-cell}
theoretical_rate = (ev(1) - s) / (ev(2) - s)
observed_rate = err(26) / err(25)
```
