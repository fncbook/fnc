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
[**Demo %s**](#demo-linear-converge)



```{code-cell}
lambda = 10;
exact = @(x) sinh(lambda * x) / sinh(lambda) - 1;
```

The following functions define the ODE.

```{code-cell}
p = @(x) zeros(size(x));            
q = @(x) -lambda^2 * ones(size(x));
r = @(x) lambda^2 * ones(size(x));
```

We compare the computed solution to the exact one for increasing $n$.

```{code-cell}
p = @(x) zeros(size(x));            
q = @(x) -lambda^2 * ones(size(x));
r = @(x) lambda^2 * ones(size(x));
n = 2 * round(10.^(1:0.25:3)');
err = zeros(size(n));
for k = 1:length(n)
    [x, u] = bvplin(p, q, r, 0, 1, -1, 0, n(k));
    err(k) = norm(exact(x) - u, Inf);
end
disp(table(n, err, variableNames = ["n", "inf-norm error"]))
```

Each factor of 10 in $n$ reduces error by a factor of 100, which is indicative of second-order convergence.

```{code-cell}
clf,  loglog(n, err, 'o-')
hold on, loglog(n, n.^(-2), 'k--')
xlabel('n'),  ylabel('max error')
title('Convergence for a linear BVP') 
legend('obs. error', '2nd order')
```
