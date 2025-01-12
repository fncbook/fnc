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
[**Demo %s**](#demo-subspace-unstable)

First we define a triangular matrix with known eigenvalues, and a random vector $b$.

```{code-cell}
lambda = 10 + (1:100);
A = diag(lambda) + triu(rand(100), 1); 
b = rand(100, 1);
```

Next we build up the first ten Krylov matrices iteratively, using renormalization after each matrix-vector multiplication.

```{code-cell}
Km = b;
for m = 1:29      
    v = A * Km(:, m);
    Km(:, m+1) = v / norm(v);
end
```

Now we solve least-squares problems for Krylov matrices of increasing dimension, recording the residual in each case.

```{code-cell}
warning off  
resid = zeros(30, 1);
for m = 1:30  
    z = (A * Km(:, 1:m)) \ b;
    x = Km(:, 1:m) * z;
    resid(m) = norm(b - A * x);
end
```

The linear system approximations show smooth linear convergence at first, but the convergence stagnates after only a few digits have been found.

```{code-cell}
clf
semilogy(resid, '.-')
xlabel('m'),  ylabel('|| b-Ax_m ||')
set(gca,'ytick',10.^(-6:2:0))
axis tight, title('Residual for linear systems') 
```
