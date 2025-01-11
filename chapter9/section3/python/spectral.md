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
[**Demo %s**](#demo-stability-spectral)

On the left, we use a log-log scale, which makes second-order algebraic convergence $O(n^{-4})$ a straight line. On the right, we use a log-linear scale, which makes spectral convergence $O(K^{-n})$ linear.

```{code-cell} 
:tags: [hide-input]
n = arange(20, 420, 20)
algebraic = 100 / n**4
spectral = 10 * 0.85**n

subplot(2, 1, 1)
loglog(n, algebraic, 'o-', label="algebraic")
loglog(n, spectral, 'o-', label="spectral")
xlabel('n'),  ylabel('error') 
title('log–log'), ylim([1e-16, 1]);
legend() 

subplot(2, 1, 2)
semilogy(n, algebraic, 'o-', label="algebraic")
semilogy(n, spectral, 'o-', label="spectral")
xlabel('n'),  ylabel('error')
title('log–linear') ,  ylim([1e-16, 1]);  
legend();
```
