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
[**Demo %s**](#demo-diffadv-vec)


```{code-cell}
m = 4;  n = 3;
x = linspace(0, 2, m+1);
y = linspace(-3, 0, n+1);

f = @(x, y) cos(0.75*pi * x .* y - 0.5*pi * y);
[mtx, X, Y, vec, unvec] = tensorgrid(x, y);
F = mtx(f);
disp("function on a 4x3 grid:")
disp(F)
```

```{code-cell}
disp("vec(F):")
disp(vec(F))
```

The `unvec` operation is the inverse of vec.

```{code-cell}
disp("unvec(vec(F)):")
disp(unvec(vec(F)))
```
