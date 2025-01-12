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
[**Demo %s**](#demo-interpolation-pwise)

Let us recall the data from {numref}`Demo %s <demo-interpolation-global>`.

```{code-cell}
n = 18;
t = linspace(-1, 1, n+1);
y = t.^2 + t + 0.05 * sin(20 * t);
clf, scatter(t, y)
```

Here is an interpolant that is linear between each consecutive pair of nodes, using `interp1` from MATLAB.

```{code-cell}
x = linspace(-1, 1, 400)';
hold on, plot(x, interp1(t, y, x))
title('Piecewise linear interpolant') 
```

We may prefer a smoother interpolant that is piecewise cubic, generated using `Spline1D` from the `Dierckx` package.

```{code-cell}
cla
scatter(t, y)
plot(x, interp1(t, y, x, 'spline'))
title('Piecewise cubic interpolant')  
```
