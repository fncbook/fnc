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
[**Demo %s**](#demo-stability-rungefix)

Here again is the function from {numref}`Demo {number} <demo-stability-runge>` that provoked the Runge phenomenon when using equispaced nodes.

```{code-cell} 
f = @(x) 1 ./ (x.^2 + 16);
```

```{code-cell} 
:tags: [hide-input]
clf
x = linspace(-1, 1, 1601)';
n = [4, 10, 16, 40];
for k = 1:length(n) 
    theta = linspace(0, pi, n(k) + 1)';
    t = -cos(theta);
    p = polyinterp(t, f(t));
    semilogy( x, abs(f(x) - p(x)) );  hold on
end
title('Error for degrees 4, 10, 16, 40')   
xlabel('x'), ylabel('|f(x)-p(x)|')   

```

By degree 16 the error is uniformly within machine epsilon, and, importantly, it stays there as $n$ increases. Note that as predicted by the error indicator function, the error is uniform over the interval at each value of $n$.
