---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
Here is the grid from {numref}`Example {number} <example-tensorprod-smallgrid>`.

```{code-cell}
m = 4
x = range(0, 2, m+1)
n = 2
y = range(1, 3, n+1);
```

For a given $f(x,y)$ we can find $\operatorname{mtx}(f)$ by using a comprehension syntax.

```{code-cell}
f = (x, y) -> cos(π * x * y - y)
F = [f(x, y) for x in x, y in y]
```

We can make a nice plot of the function by first choosing a much finer grid. However, the contour and surface plotting functions expect the *transpose* of mtx($f$).
```{tip}
To emphasize departures from a zero level, use a colormap such as `redsblues`, and use `clims` to set balanced color differences.
```

```{code-cell}
using Plots
m, n = 80, 60
x = range(0, 2, m+1);
y = range(1, 3, n+1);
F = [f(x, y) for x in x, y in y]
contour(x, y, F';
    levels=21,  aspect_ratio=1,
    color=:redsblues,  clims=(-1, 1),
    xlabel="x",  ylabel="y" )
``` 

