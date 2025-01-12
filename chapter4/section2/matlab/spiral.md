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
[**Demo %s**](#demo-fp-spiral)

Let's convert the roots of a quadratic polynomial $f(x)$ to a fixed point problem.

```{code-cell}
f = @(x) x.^2 - 4*x + 3.5;
r = roots([1, -4, 3.5])
```

We define $g(x)=x-p(x)$. 

```{code-cell}
g = @(x) x - f(x);
```

Intersections of $y=g(x)$ with the line $y=x$ are fixed points of $g$ and thus roots of $f$. (Only one is shown in the chosen plot range.)

```{code-cell}
clf
fplot(g, [2, 3])
hold on,  plot([2, 3], [2, 3], 'k')
title('Finding a fixed point'),  axis equal  
xlabel('x'),  ylabel('y')  
```

If we evaluate $g(2.1)$, we get a value of almost 2.6, so this is not a fixed point.

```{code-cell}
x = 2.1;
y = g(x)
```

However, $y=g(x)$ is considerably closer to the fixed point at around 2.7 than $x$ is. Suppose then that we adopt $y$ as our new $x$ value. Changing the $x$ coordinate in this way is the same as following a horizontal line over to the graph of $y=x$.

```{code-cell}
plot([x, y], [y, y], '-')
x = y;
```

Now we can compute a new value for $y$. We leave $x$ alone here, so we travel along a vertical line to the graph of $g$.

```{code-cell}
y = g(x)
plot([x, x],[x, y], '-')
```

You see that we are in a position to repeat these steps as often as we like. Let's apply them a few times and see the result.

```{code-cell}
for k = 1:5
    plot([x, y], [y, y], '-')
    x = y;       % y --> new x
    y = g(x);    % g(x) --> new y
    plot([x, x], [x, y], '-')  
end
```

The process spirals in beautifully toward the fixed point we seek. Our last estimate has almost 4 accurate digits.

```{code-cell} 
abs(y - r(1)) / r(1)
```

Now let's try to find the other fixed point $\approx 1.29$ in the same way. We'll use 1.3 as a starting approximation.

```{code-cell}
cla
fplot(g, [1, 2])
hold on, plot([1, 2], [1, 2], 'k')
ylim([1, 2])
x = 1.3;  y = g(x);
for k = 1:5
    plot([x, y], [y, y], '-'),  
    x = y;       % y --> new x
    y = g(x);    % g(x) --> new y
    plot([x, x], [x, y], '-') 
end
title('No convergence') 
```

This time, the iteration is pushing us _away from_ the correct answer.
