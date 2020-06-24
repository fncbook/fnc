# FundamentalsNumericalComputation.jl

These are Julia implementations the functions that are defined in the first edition of [*Fundamentals of Numerical Computation*](https://tobydriscoll.net/fnc) by Driscoll and Braun. They are mostly compatible with the original MATLAB functions, except that input and output formats have been changed to be compatible with major packages listed below.

The package does not export its own functions, but it does export the name `FNC` as a shorthand for the full package name. It does re-export the following:

* LinearAlgebra
* Plots
* SparseArrays
* Polynomials
* NLsolve
* DifferentialEquations
* DataFrames
* Interpolations

**This is work in progress**. There will be cosmetic changes as well as possible breaking changes in future versions. In particular, the functions from Chapters 7-12 have not yet been thoroughly reviewed and tested. No warranty for correctness is implied.
