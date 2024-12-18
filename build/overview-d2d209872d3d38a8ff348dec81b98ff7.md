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
# Introduction

```{index} Yoda, The Empire Strikes Back
```

```{epigraph}
You must unlearn what you have learned. 

— Yoda, _The Empire Strikes Back_
```

Our first step is to discretize the real numbers—specifically, to replace them with a finite surrogate set of numbers. This step keeps the time and storage requirements for operating with each number at constant levels, but virtually every data set and arithmetic operation is perturbed slightly away from its idealized mathematical value. We can easily keep the individual roundoff errors very small, so small that simple random accumulation is unlikely to bother us. However, some problems are extremely sensitive to these perturbations, a trait we quantify using a *condition number*. Problems with large condition numbers are difficult to solve accurately using finite precision. Furthermore, even when the condition number of a problem is not large, some algorithms for solving it allow errors to grow enormously. We call these algorithms *unstable*. In this chapter we discuss these ideas in simple settings before moving on to the more realistic problems in the rest of the book.

**Software**

Instructions for obtaining Julia and the codes used throughout the text can be found at

[https://github.com/fncbook/FundamentalsNumericalComputation.jl](https://github.com/fncbook/FundamentalsNumericalComputation.jl)

The installation process, which can take 5-10 minutes, only needs to be performed once.