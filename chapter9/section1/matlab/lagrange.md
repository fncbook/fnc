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
[**Demo %s**](#demo-polynomial-lagrange)

Here is a vector of nodes.

```{code-cell}
t = [ 1, 1.5, 2, 2.25, 2.75, 3 ];
n = 5;  k = 2;
not_k = [0:k-1 k+1:n];   % all except the kth node
```

Let's apply the definition of the cardinal Lagrange polynomial for $k=2$. First we define a polynomial $q$ that is zero at all the nodes except $i=k$. Then $\ell_2$ is found by normalizing $q$ by $q(t_k)$.
```{tip}
:class: dropdown
Whenever we index into the node vector `t`, we have to add 1 since the mathematical index starts at zero.
```

```{code-cell}
q = @(x) prod(x - t(not_k + 1));
ell_k = @(x) q(x) ./ q(t(k + 1));

```

A plot confirms the cardinal property of the result.

```{code-cell}
clf
fplot(ell_k, [1, 3])
hold on, grid on
plot(t(not_k + 1), 0 * t(not_k + 1), 'o')
plot(t(k + 1), 1, 'o')
xlabel('x'),  ylabel('\ell_2(x)')    
title('Lagrange cardinal function')   
```

Observe that $\ell_k$ is _not_ between zero and one everywhere, unlike a hat function.
