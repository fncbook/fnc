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
[**Demo %s**](#demo-secant-line)



We return to finding a root of the equation $x e^x=2$.

```{code-cell}
f = @(x) x .* exp(x) - 2;
clf, fplot(f, [0.25, 1.25])
set(gca, 'ygrid', 'on')  
xlabel('x'), ylabel('y')    
title('Objective function')    
```

From the graph, it's clear that there is a root near $x=1$. To be more precise, there is a root in the interval $[0.5,1]$. So let us take the endpoints of that interval as _two_ initial approximations.

```{code-cell}
x1 = 1;    y1 = f(x1);
x2 = 0.5;  y2 = f(x2);
hold on, scatter([x1, x2], [y1, y2])
title('Two initial values')    
```

Instead of constructing the tangent line by evaluating the derivative, we can construct a linear model function by drawing the line between the two points $\bigl(x_1,f(x_1)\bigr)$ and $\bigl(x_2,f(x_2)\bigr)$. This is called a _secant line_.

```{code-cell}
slope2 = (y2 - y1) / (x2 - x1);
secant2 = @(x) y2 + slope2 * (x - x2);
```

As before, the next root estimate in the iteration is the root of this linear model.

```{code-cell}
fplot(secant2,[0.25, 1.25],'k--')
x3 = x2 - y2 / slope2;
y3 = f(x3)
scatter(x3, 0)
title('Next value') 
```

For the next linear model, we use the line through the two most recent points. The next iterate is the root of that secant line, and so on.

```{code-cell}
slope2 = (y3 - y2) / (x3 - x2);
x4 = x3 - y3 / slope2;
y4 = f(x4)
```
