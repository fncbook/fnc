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
[**Demo %s**](#demo-pwlin-usage)

We generate a piecewise linear interpolant of $f(x)=e^{\sin 7x}$.

```{code-cell}
f = @(x) exp(sin(7 * x));
clf
fplot(f, [0, 1], displayname="function")
xlabel("x");  ylabel(("y"));
```

First we sample the function to create the data.

```{code-cell}
t = [0, 0.075, 0.25, 0.55, 0.7, 1];    % nodes
y = f(t);                              % function values
```

Now we create a callable function that will evaluate the piecewise linear interpolant at any $x$, and then plot it.

```{code-cell}
p = plinterp(t, y);
hold on
fplot(p, [0, 1], displayname="interpolant")
scatter(t, y, displayname="values at nodes")
title("PL interpolation")
legend();
```
