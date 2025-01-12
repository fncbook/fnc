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
[**Demo %s**](#demo-improper-decay)

```{code-cell}
:tags: [hide-input]
f = @(x) 1 ./ (1 + x.^2);
clf,  subplot(2, 1, 1)
fplot(f, [-4, 4]);  set(gca, 'yscale', 'log') 
xlabel('x'),  ylabel('f(x)'),  ylim([1e-20, 1])  
title('Original integrand')   

x = @(t) sinh( pi * sinh(t) / 2 );
chain = @(t) pi/2 * cosh(t) .* cosh( pi * sinh(t) / 2 );
integrand = @(t) f(x(t)) .* chain(t);
subplot(2, 1, 2)
fplot(integrand, [-4, 4]);  set(gca, 'yscale', 'log') 
xlabel('t'), ylabel('f(x(t))'),  ylim([1e-20, 1])  
title('Transformed integrand')  
```

This graph suggests that we capture all of the integrand values that are larger than machine epsilon by integrating in $t$ from $-4$ to $4$.
