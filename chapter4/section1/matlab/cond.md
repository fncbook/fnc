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
[**Demo %s**](#demo-roots-cond)

Consider first the function

```{code-cell}
f  = @(x) (x - 1) .* (x - 2);
```

At the root $r=1$, we have $f'(r)=-1$. If the values of $f$ were perturbed at every point by a small amount of noise, we can imagine finding the root of the function drawn with a thick ribbon, giving a range of potential roots.


```{code-cell}
clf
interval = [0.8, 1.2];
fplot(f, interval)
grid on, hold on
fplot(@(x) f(x) + 0.02, interval, 'k')
fplot(@(x) f(x) - 0.02, interval, 'k')
axis equal 
xlabel('x'), ylabel('f(x)')
title(('Well-conditioned root'));
```

The possible values for a perturbed root all lie within the interval where the ribbon intersects the $x$-axis. The width of that zone is about the same as the vertical thickness of the ribbon.

By contrast, consider the function

```{code-cell}
f = @(x) (x - 1) .* (x - 1.01);
```

Now $f'(1)=-0.01$, and the graph of $f$ will be much shallower near $x=1$. Look at the effect this has on our thick rendering:

```{code-cell}
axis(axis), cla
fplot(f, interval)
fplot(@(x) f(x) + 0.02, interval, 'k')
fplot(@(x) f(x) - 0.02, interval, 'k')
title('Poorly conditioned root')   
```

The vertical displacements in this picture are exactly the same as before. But the potential _horizontal_ displacement of the root is much wider. In fact, if we perturb the function entirely upward by the amount drawn here, the root disappears!
