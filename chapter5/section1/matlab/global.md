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
[**Demo %s**](#demo-interpolation-global)

Here are some points that we could consider to be observations of an unknown function on $[-1,1]$.

```{code-cell}
n = 5;
t = linspace(-1,1,n+1)';  
y = t.^2 + t + 0.05 * sin(20 * t);
clf, scatter(t,y)
```

```{index} ! Julia; fit
```

The polynomial interpolant, as computed using `polyfit`, looks very sensible. It's the kind of function you'd take home to meet your parents.

```{code-cell}
c = polyfit(t, y, n);     % polynomial coefficients
p = @(x) polyval(c, x);
hold on
fplot(p, [-1 1])
legend('data', 'interpolant', 'location', 'north');
```

But now consider a different set of points generated in almost exactly the same way.

```{code-cell}
n = 18;
t = linspace(-1, 1, n+1);
y = t.^2 + t + 0.05 * sin(20 * t);
clf, scatter(t, y)
```

The points themselves are unremarkable. But take a look at what happens to the polynomial interpolant.

```{code-cell}
c = polyfit(t, y, n);     % polynomial coefficients
p = @(x) polyval(c, x);
hold on, fplot(p, [-1 1])
legend('data', 'interpolant', 'location', 'north');
```

Surely there must be functions that are more intuitively representative of those points!
