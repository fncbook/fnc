---
jupytext:
  cell_metadata_filter: -all
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.10.3
kernelspec:
  display_name: Julia 1.7.1
  language: julia
  name: julia-fast
---
# Krylov methods in linear algebra

```{index} Luke Skywalker, Return of the Jedi
```

:::{epigraph}
I warn you not to underestimate my powers.

— Luke Skywalker, *Return of the Jedi* 
:::

What are the implications of the $O(n^3)$ work requirements for solving linear systems? Suppose tomorrow your computer became a thousand times faster. (Historically this has taken about 15 years in the real world.) Assuming you are willing to wait just as long today as you were yesterday, the size of the linear system you can solve has gone up only by a factor of 10. Nice, but not nearly the jump that you got in hardware power. In fact, there is an odd paradox: faster computers make faster algorithms *more* important, not less, because they demand that you work at larger values of $n$, where asymptotic differences are large.

In practice the only reasonable way to deal with large matrices (at this writing, $n>10^4$ or so) is if they are sparse, or can be approximated sparsely. But LU factorization of a sparse matrix does not necessarily lead to sparse factors, particularly when row pivoting is required. The algorithm can be improved to be more sparse-aware, but we will not go into the details.

Instead, we will replace LU factorization with an iterative algorithm. Unlike the LU factorization, iteration gives useful intermediate and continually improving results before the exact solution is found, allowing us to stop well before the nominal exact termination. More importantly, though, these iterations, based on an idea called *Krylov subspaces*, allow us to fully exploit sparsity.

Krylov subspace methods have two other advantages that are subtle but critically relevant to applications. One is that they allow us to do linear algebra *even without having the relevant matrix*. This may sound undesirable or even impossible, but it exploits the connection between matrix-vector multiplication and a linear transformation. The other major advantage of Krylov subspace iterations is that they can exploit approximate inverses when they are available. These two features are among the most powerful ideas behind scientific computation today.

