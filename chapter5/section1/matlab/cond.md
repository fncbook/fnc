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
[**Demo %s**](#demo-interp-cond)

In {numref}`Demo %s <demo-interpolation-global>` and {numref}`Demo %s <demo-interpolation-pwise>` we saw a big difference between polynomial interpolation and piecewise polynomial interpolation of some arbitrarily chosen data. The same effects can be seen clearly in the cardinal functions, which are closely tied to the condition numbers.

```{code-cell}
n = 18;
t = linspace(-1, 1, n+1)';
y = [zeros(9, 1); 1; zeros(n - 9, 1)];    % 10th cardinal function
clf, scatter(t, y)
hold on
x = linspace(-1, 1, 400)';
plot(x, interp1(t, y, x, 'spline'))
title('Piecewise cubic cardinal function') 
xlabel('x'), ylabel('p(x)') 
```

The piecewise cubic cardinal function is nowhere greater than one in absolute value. This happens to be true for all the cardinal functions, ensuring a good condition number for any interpolation with these functions. But the story for global polynomials is very different.

```{code-cell}
clf, scatter(t, y)
c = polyfit(t, y, n);
hold on, plot(x, polyval(c, x))
title('Polynomial cardinal function')
xlabel('x'), ylabel(('p(x)'));
```

From the figure we can see that the condition number for polynomial interpolation on these nodes is at least 500.
