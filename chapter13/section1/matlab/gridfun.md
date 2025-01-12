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
[**Demo %s**](#demo-tensorprod-gridfun)

Here is the grid from {numref}`Example {number} <example-tensorprod-smallgrid>`.

```{code-cell}
m = 4;   
x = linspace(0, 2, m+1);
n = 2;   
y = linspace(1, 3, n+1);
```

For a given $f(x,y)$ we can find $\operatorname{mtx}(f)$ by using a comprehension syntax.

```{code-cell}
[mtx, X, Y] = tensorgrid(x, y);
f = @(x, y) cos(pi * x .* y - y);
F = mtx(f)
```

We can make a nice plot of the function by first choosing a much finer grid. However, the contour and surface plotting functions expect the *transpose* of mtx($f$).
```{tip}
:class: dropdown
To emphasize departures from a zero level, use a colormap such as `redsblues` and set the color limits to be symmetric around zero.
```

::::{warning}
The contour and surface plotting functions expect the *transpose* of the outputs of `mtx`. If you forget to do that, the $x$ and $y$ axes will be swapped.
::::

```{code-cell}
m = 80;  x = linspace(0, 2, m+1);
n = 60;  y = linspace(1, 3, n+1);
[mtx, X, Y] = tensorgrid(x, y);
F = mtx(f);

pcolor(X', Y', F')
shading interp
colormap(redsblues),  colorbar
axis equal
xlabel("x"),  ylabel("y")
```
