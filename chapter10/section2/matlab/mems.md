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
[**Demo %s**](#demo-shooting-mems)

We revisit {numref}`Demo {number} <demo-shooting-naive>` but let {numref}`Function {number} <function-shoot>` do the heavy lifting.

```{code-cell}
lambda = 0.6;
phi = @(r, w, dwdr) lambda ./ w.^2 - dwdr ./ r;   
a = eps;  b = 1;    % avoid r=0 in denominator
```

We specify the given and unknown endpoint values.

```{code-cell}
ga = @(u, du) du;
gb = @(u, du) u - 1;
```

```{code-cell}
init = [0.8; 0];    % initial guess for u(a) and u'(a)
[r, w, dwdx] = shoot(phi, a, b, ga, gb, init, 1e-5);
clf,  plot(r, w)
title('Correct solution')
xlabel('r'),  ylabel('w(r)')
```

The value of $w$ at $r=1$, meant to be exactly one, was computed to be

```{code-cell}
format long
w(end)
```

The accuracy is consistent with the error tolerance used for the IVP solution. The initial value $w(0)$ that gave this solution is

```{code-cell}
w(1)
```
