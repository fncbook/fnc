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
[**Demo %s**](#demo-pivoting-builtin)

With the syntax `A \ b`, the matrix `A` is PLU-factored, followed by two triangular solves.

```{code-cell}
A = randn(500, 500);    % 500x500 with normal random entries
tic; for k=1:50; A \ rand(500, 1); end; toc
```

In {numref}`section-linsys-efficiency` we showed that the factorization is by far the most costly part of the solution process. A factorization object allows us to do that costly step only once per unique matrix. 

```{code-cell}
[L, U, p] = lu(A, 'vector');    % keep factorization result
tic
for k=1:50
    b = rand(500, 1);
    U \ (L \ b(p));
end
toc
```
