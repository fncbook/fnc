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
[**Demo %s**](#demo-stability-spectral)

On the left, we use a log-log scale, which makes second-order algebraic convergence $O(n^{-4})$ a straight line. On the right, we use a log-linear scale, which makes spectral convergence $O(K^{-n})$ linear.

```{code-cell} 
:tags: [hide-input]
n = (20:20:400)';
algebraic = 100 ./ n.^4;
spectral = 10 * 0.85.^n;
clf, subplot(2, 1, 1)
loglog(n, algebraic, 'o-', displayname="algebraic")
hold on;  loglog(n, spectral, 'o-', displayname="spectral")
xlabel('n'),  ylabel('error')   
title('log–log')   
axis tight,  ylim([1e-16, 1]);  legend(location="southwest")   

subplot(2, 1, 2)
semilogy(n, algebraic, 'o-', displayname="algebraic")
hold on;  semilogy(n, spectral, 'o-', displayname="spectral")
xlabel('n'), ylabel('error'),  ylim([1e-16, 1])   
title('log–linear')   
axis tight,  ylim([1e-16, 1]);  legend(location="southwest")   
```
