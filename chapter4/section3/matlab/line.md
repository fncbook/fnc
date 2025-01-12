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
[**Demo %s**](#demo-newton-line)


Suppose we want to find a root of the function

```{code-cell}
f = @(x) x .* exp(x) - 2;
clf, fplot(f, [0, 1.5])
xlabel('x'), ylabel('y')    
set(gca, 'ygrid', 'on')  
title('Objective function')    

```

From the graph, it is clear that there is a root near $x=1$. So we call that our initial guess, $x_1$.

```{code-cell}
x1 = 1;
y1 = f(x1)
hold on, scatter(x1, y1)
```

Next, we can compute the tangent line at the point $\bigl(x_1,f(x_1)\bigr)$, using the derivative.

```{code-cell}
dfdx = @(x) exp(x) .* (x + 1);
slope1 = dfdx(x1);
tangent1 = @(x) y1 + slope1 * (x - x1);
fplot(tangent1, [0, 1.5],'--')
title('Function and tangent line')    

```

In lieu of finding the root of $f$ itself, we settle for finding the root of the tangent line approximation, which is trivial. Call this $x_2$, our next approximation to the root.

```{code-cell}
x2 = x1 - y1 / slope1
scatter(x2, 0)
title('Root of the tangent')    
```

```{code-cell}
y2 = f(x2)
```

The residual (i.e., value of $f$) is smaller than before, but not zero. So we repeat the process with a new tangent line based on the latest point on the curve.

```{code-cell}
cla, fplot(f, [0.8, 0.9])
scatter(x2, y2)
slope2 = dfdx(x2);
tangent2 = @(x) y2 + slope2 * (x - x2);
fplot(tangent2, [0.8, 0.9], '--')
x3 = x2 - y2 / slope2;
scatter(x3, 0)
title('Next iteration')    
```

```{code-cell}
y3 = f(x3)
```

Judging by the residual, we appear to be getting closer to the true root each time.
