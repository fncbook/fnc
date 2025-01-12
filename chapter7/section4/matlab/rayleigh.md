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
[**Demo %s**](#demo-symm-eig-rayleigh)

We will use a symmetric matrix with a known EVD and eigenvalues equal to the integers from 1 to 20.

```{code-cell}
n = 20;
lambda = 1:n;
D = diag(lambda);
[V, ~] = qr(randn(n, n));    % get a random orthogonal V
A = V * D * V';
```

The Rayleigh quotient is a scalar-valued function of a vector.

```{code-cell}
R = @(x) (x' * A * x) / (x' * x);
```

The Rayleigh quotient evaluated at an eigenvector gives the corresponding eigenvalue.

```{code-cell}
format long
R(V(:, 7))
```

If the input to he Rayleigh quotient is within a small $\delta$ of an eigenvector, its output is within $O(\delta^2)$ of the corresponding eigenvalue. In this experiment, we observe that each additional digit of accuracy in an approximate eigenvector gives two more digits to the eigenvalue estimate coming from the Rayleigh quotient.

```{code-cell}
delta = 1 ./ 10 .^ (1:5)';
dif = zeros(size(delta));
for k = 1:length(delta)
    e = randn(n, 1);
    e = delta(k) * e / norm(e);
    x = V(:, 6) + e;
    dif(k) = R(x) - lambda(6);
end
disp(table(delta, dif, variablenames=["perturbation size", "R(x) - lambda"]))
```
