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
[**Demo %s**](#demo-stability-errfun)

We plot $|\Phi(x)|$ over the interval $[-1,1]$ with equispaced nodes for different values of $n$. 

```{code-cell} 
:tags: [hide-input]
clf
x = linspace(-1, 1, 1601)';
Phi = zeros(size(x));
for n = 10:10:50
    t = linspace(-1, 1, n+1)';
    for k = 1:length(x)
        Phi(k) = prod(x(k) - t);
    end
    semilogy(x, abs(Phi)),  hold on
end
title('Error indicator on equispaced nodes')    
xlabel('x'),  ylabel('|\Phi(x)|')   
```

Each time $\Phi$ passes through zero at an interpolation node, the value on the log scale should go to $-\infty$, which explains the numerous cusps on the curves.
