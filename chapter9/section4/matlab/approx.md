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
[**Demo %s**](#demo-orthogonal-approx)

Let's approximate $e^x$ over the interval $[−1,1]$. We can sample it at, say, 15 points, and find the best-fitting straight line to that data.

```{code-cell}
clf;  fplot(@exp, [-1, 1], displayname="function")
t = linspace(-1, 1, 15)';
y = exp(t);
V = [t.^0, t];
c = V \ y;
p = @(t) c(1) + c(2)*t;

hold on,  fplot(p, [-1, 1], displayname="LS fit at 15 points")
title('Least-squares fit to samples of exp(x)')    
xlabel('x'),  ylabel('f(x)')    
legend(location="northwest")    
```

There's nothing special about 15 points. Choosing more doesn't change the result much.

```{code-cell}
t = linspace(-1, 1, 150)';
y = exp(t);
V = [t.^0, t];
c = V \ y;
p = @(t) c(1) + c(2)*t;
fplot(p, [-1, 1], displayname="LS fit at 150 points")
```

This situation is unlike interpolation, where the degree of the interpolant increases with the number of nodes. Here, the linear fit is apparently approaching a limit that we may think of as a continuous least-squares fit.

```{code-cell}
n = (40:60:400)';
slope = zeros(size(n));
intercept = zeros(size(n));

for k = 1:length(n)
    t = linspace(-1, 1, n(k))';
    V = [t.^0, t];
    c = V \ exp(t);
    intercept(k) = c(1);
    slope(k) = c(2);
end
disp(table(n, intercept, slope))
```
