---
kernelspec:
  display_name: Python 3
  language: python
  name: python3
numbering:
  headings: false
---
```{code-cell}
:tags: [remove-cell]
exec(open("../../../python/FNC_init.py").read())
```
[**Demo %s**](#demo-improper-decay)

```{code-cell}
f = lambda x: 1 / (1 + x**2)
x = linspace(-4, 4, 500)
subplot(2, 1, 1)
plot(x, f(x)),  yscale('log')
xlabel('x'),  ylabel('f(x)'),  ylim([1e-16, 1])  
title('Original integrand')   

xi = lambda t: sinh( pi * sinh(t) / 2 )
dxi_dt = lambda t: pi/2 * cosh(t) * cosh( pi * sinh(t) / 2 )
integrand = lambda t: f(xi(t)) * dxi_dt(t)
subplot(2, 1, 2)
plot(x, integrand(x)),  yscale('log')
xlabel('t'),  ylabel('f(x(t))'),  ylim([1e-16, 1])  
title('Transformed integrand')   
```

This graph suggests that we capture all of the integrand values that are larger than machine epsilon by integrating in $t$ from $-4$ to $4$.
