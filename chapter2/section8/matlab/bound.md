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
[**Demo %s**](#demo-condition-bound)


```{index} ! MATLAB; cond
```
MATLAB has a function `cond` to compute matrix condition numbers. By default, the 2-norm is used. As an example, the family of *Hilbert matrices* is famously badly conditioned. Here is the $6\times 6$ case.

```{code-cell}
A = hilb(6)
kappa = cond(A)
```

Because $\kappa\approx 10^8$, it's possible to lose nearly 8 digits of accuracy in the process of passing from $\mathbf{A}$ and $\mathbf{b}$ to $\mathbf{x}$. That fact is independent of the algorithm; it's inevitable once the data are expressed in finite precision. 

Let's engineer a linear system problem to observe the effect of a perturbation. We will make sure we know the exact answer.

```{code-cell}
x = (1:6)';
b = A * x;
```

Now we perturb the system matrix and vector randomly by $10^{-10}$ in norm.

```{code-cell} 
dA = randn(size(A));  dA = 1e-10 * (dA / norm(dA));
db = randn(size(b));  db = 1e-10 * (db / norm(db));
```

We solve the perturbed problem using pivoted LU and see how the solution was changed.

```{code-cell}
new_x = ((A + dA) \ (b + db));
dx = new_x - x;
```

Here is the relative error in the solution.

```{code-cell}
relative_error = norm(dx) / norm(x)
```

And here are upper bounds predicted using the condition number of the original matrix.

```{code-cell}
upper_bound_b = (kappa * norm(db) / norm(b))
upper_bound_A = (kappa * norm(dA) / norm(A))
```

Even if we didn't make any manual perturbations to the data, machine roundoff does so at the relative level of $\macheps$.

```{code-cell}
dx = A\b - x;
relative_error = norm(dx) / norm(x)
rounding_bound = kappa * eps
```

Larger Hilbert matrices are even more poorly conditioned:

```{code-cell}
A = hilb(14);
kappa = cond(A)
```

Note that $\kappa$ exceeds $1/\macheps$. In principle we therefore may end up with an answer that has relative error greater than 100%.

```{code-cell}
rounding_bound = kappa * eps
```

Let's put that prediction to the test.

```{code-cell}
x = (1:14)';  b = A * x;
dx = A\b - x;
relative_error = norm(dx) / norm(x)
```

As anticipated, the solution has zero accurate digits in the 2-norm.
