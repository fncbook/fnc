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
[**Demo %s**](#demo-stability-runge)

This function has infinitely many continuous derivatives on the entire real line and looks easy to approximate over $[-1,1]$.

```{code-cell} 
f = @(x) 1 ./ (x.^2 + 16);
clf,  fplot(f, [-1, 1])
xlabel('x'),  ylabel('f(x)')    
title('Test function')    
```

We start by doing equispaced polynomial interpolation for some small values of $n$.

```{code-cell} 
:tags: [hide-input]
x = linspace(-1, 1, 1601)';
n = (4:4:12)';
for k = 1:length(n)
    t = linspace(-1, 1, n(k) + 1)';        % equally spaced nodes
    p = polyinterp(t, f(t));
    semilogy(x, abs(f(x) - p(x)));  hold on
end
title('Error for degrees 4, 8, 12')   
xlabel('x'), ylabel('|f(x) - p(x)|')   
```

The convergence so far appears rather good, though not uniformly so. However, notice what happens as we continue to increase the degree.

```{code-cell} 
:tags: [hide-input]
n = 12 + 15 * (1:3);
clf
for k = 1:length(n)
    t = linspace(-1, 1, n(k) + 1)';        % equally spaced nodes
    p = polyinterp(t, f(t));
    semilogy(x, abs(f(x) - p(x)));  hold on
end
title('Error for degrees 27, 42, 57')   
xlabel('x'), ylabel('|f(x) - p(x)|')   
```

The convergence in the middle can't get any better than machine precision relative to the function values. So maintaining the growing gap between the center and the ends pushes the error curves upward exponentially fast at the ends, wrecking the convergence.
