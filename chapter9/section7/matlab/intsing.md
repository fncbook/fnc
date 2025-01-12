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
[**Demo %s**](#demo-improper-intsing)


```{code-cell}
:tags: [hide-input]
f = @(x) 1 ./ (10 * sqrt(x));
tol = 1 ./ 10.^(5:0.5:14);
err = zeros(length(tol), 2);
len = zeros(length(tol), 2);
for k = 1:length(tol)
    [I1, x1] = intadapt(f, (tol(k)/20)^2, 1, tol(k));
    [I2, x2] = intsing(f, tol(k));
    err(k, :) = abs(0.2 - [I1, I2]);
    len(k, :) = [length(x1), length(x2)];
end
clf,  loglog(len, err, 'o-')   
n = [30, 3000];
hold on,  loglog(n, 30 * n.^(-4), 'k--')  % 4th order error
legend("direct", "double exponential", "4th order", location="southwest")
title(("Comparison of integration methods"));
```

As in {numref}`Demo {number} <demo-improper-intinf>`, the double exponential method is more accurate than direct integration by a few orders of magnitude. Equivalently, the same accuracy can be reached with many fewer nodes.
