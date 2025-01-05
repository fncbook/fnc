---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
With the syntax `A \ b`, the matrix `A` is PLU-factored, followed by two triangular solves.

```{code-cell}
A = randn(500, 500)   # 500x500 with normal random entries
A \ rand(500)          # force compilation
@elapsed for k=1:50; A \ rand(500); end
```

In {numref}`section-linsys-efficiency` we showed that the factorization is by far the most costly part of the solution process. A factorization object allows us to do that costly step only once per unique matrix. 

```{code-cell}
factored = lu(A)     # store factorization result
factored \ rand(500)   # force compilation
@elapsed for k=1:50; factored \ rand(500); end
```
