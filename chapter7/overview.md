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
# Matrix analysis

```{index} Yoda, The Empire Strikes Back
```

:::{epigraph}
Judge me by my size, do you?

— Yoda, *The Empire Strikes Back* 
:::

In previous chapters, we have seen how matrices that represent square or overdetermined linear systems of equations can be manipulated into LU and QR factorizations. But matrices have other factorizations that are more intrinsic to their nature as mathematical linear transformations. The most fundamental of these are the eigenvalue and singular value decompositions.

These decompositions can be used to solve linear and least-squares systems, but they have greater value in how they represent the matrix itself. They lead to critical and quantitative insights about the structure of the underlying transformation and suggest ways to approximate it efficiently. In this chapter, we will look at both of these fundamental decompositions and hint at just a few of their computational applications.

