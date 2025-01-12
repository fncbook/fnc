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
[**Demo %s**](#demo-stability-errcheb)

Now we look at the error indicator function $\Phi$ for Chebyshev node sets.

```{code-cell} 
:tags: [hide-input]
clf
x = linspace(-1, 1, 1601)';
Phi = zeros(size(x));
for n = 10:10:50
    theta = linspace(0, pi, n+1)';
    t = -cos(theta);                    
    for k = 1:length(x)
        Phi(k) = prod(x(k) - t);
    end
    semilogy(x, abs(Phi));  hold on
end
axis tight, title('Effect of Chebyshev nodes')    
xlabel('x'), ylabel('|\Phi(x)|')   
ylim([1e-18, 1e-2])   
```

In contrast to the equispaced case, $|\Phi|$ decreases exponentially with $n$ almost uniformly across the interval.
