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
[**Demo %s**](#demo-secant-iqi)

Here we look for a root of $x+\cos(10x)$ that is close to 1.

```{code-cell}
f = @(x) x + cos(10 * x);
interval = [0.5, 1.5];
clf, fplot(f, interval)
set(gca, 'ygrid', 'on'), axis(axis)   
title('Objective function')    
xlabel('x'), ylabel('y')    
r = fzero(f, 1)
```

We choose three values to get the iteration started.

```{code-cell}
x = [0.8, 1.2, 1]';
y = f(x);
hold on, scatter(x, y)
title('Three initial points')    
```

If we were using forward interpolation, we would ask for the polynomial interpolant of $y$ as a function of $x$. But that parabola has no real roots.

```{code-cell}
c = polyfit(x, y, 2);    % coefficients of interpolant
q = @(x) polyval(c, x);
fplot(q, interval, '--')
title('Parabola model')     
```

To do inverse interpolation, we swap the roles of $x$ and $y$ in the interpolation.

```{tip}
:class: dropdown
By giving two functions in the `fplot` call, we get the parametric plot $(q(y),y)$ as a function of $y$.
```

```{code-cell}
cla, fplot(f, interval)
scatter(x, y)     
c = polyfit(y, x, 2);    % coefficients of interpolating polynomial
q = @(y) polyval(c, y);
fplot(q, @(y) y, ylim,'--')    % plot x=q(y), y=y
title('Sideways parabola')    
```

We seek the value of $x$ that makes $y$ zero. This means evaluating $q$ at zero.

```{code-cell}
x = [x; q(0)];
y = [y; f(x(end))]
```

We repeat the process a few more times.

```{code-cell}
for k = 4:8
    c = polyfit(y(k-2:k), x(k-2:k), 2);
    x(k+1) = polyval(c, 0);
    y(k+1) = f(x(k+1));
end
disp('final residual:')
y(end)
```

Here is the sequence of errors.

```{code-cell}
format short e
err = x - r
```

The convergence is probably superlinear:

```{code-cell}
logerr = log(abs(err));
ratios = logerr(2:end) ./ logerr(1:end-1)
```

