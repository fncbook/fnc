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
[**Demo %s**](#demo-improper-intinf)

```{code-cell}
:tags: [hide-input]
f = @(x) 1 ./ (1 + x.^2);
tol = 1 ./ 10.^(5:0.5:14);
err = zeros(length(tol), 2);
len = zeros(length(tol), 2);
for k = 1:length(tol)
    [I1, x1] = intadapt(f, -2/tol(k), 2/tol(k), tol(k));
    [I2, x2] = intinf(f, tol(k));
    err(k, :) = abs(pi - [I1, I2]);
    len(k, :) = [length(x1), length(x2)];
end
clf,  loglog(len, err, 'o-')   
n = [100, 10000];
hold on,  loglog(n, 1000 * n.^(-4), 'k--')  % 4th order error
legend("direct", "double exponential", "4th order", location="southwest")
title(("Comparison of integration methods"));
```

Both methods are roughly fourth-order due to Simpson's formula in the underlying adaptive integration method. At equal numbers of evaluation nodes, however, the double exponential method is consistently 2–3 orders of magnitude more accurate.
