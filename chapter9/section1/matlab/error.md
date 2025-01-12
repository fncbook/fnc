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
[**Demo %s**](#demo-polynomial-error)


```{code-cell}
t =  [ 1, 1.6, 1.9, 2.7, 3 ];
n = length(t) - 1;
Phi = @(x) prod(x - t);

clf,  fplot(@(x) Phi(x) / 5, [1, 3])
hold on,  plot(t, 0*t, 'o')
xlabel('x'),  ylabel('\Phi(x)')   
title('Interpolation error function')   
```

The error is zero at the nodes, by the definition of interpolation. The error bound, as well as the error itself, has one local maximum between each consecutive pair of nodes.
