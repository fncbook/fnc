---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering: false
---
# How to use this book

## For everyone

First, set up your computing environment:

- [Julia](#section-setup-julia)
- [MATLAB](#section-setup-matlab)
- [Python](#section-setup-python)

This book is littered with computer examples. All the code used to generate the content you see is given to you. Here is a Julia example:

```{code-cell}
println("Welcome to Julia! Do you know π to 100 digits? Because I do! Look:")
setprecision(328)
BigFloat(π)
```

If you hover over the code cell above, you will see an icon in the upper right corner that allows you to copy that code to your clipboard. You will probably do this kind of thing a lot as you work though this book. Many exercises are simply variations on the examples presented in the text.

One thing you do *not* have to copy is the code for the functions that are presented in the book as implementations of key algorithms. They will all be available once you install them, as described in the links above.[^goodlook]

[^goodlook]: You didn't skip the "first, do this" part, did you? Not a great look for you.

Here is a short overview of major features of the book.

![usage video](_static/FNC-usage.mp4)

---

## For instructors

If you find errors or want to suggest improvements, please [open an issue](https://github.com/fncbook/fnc/issues/new/choose).

My publisher, SIAM, is cautious about how this kind of text affects sales of the print book. If you are an instructor, please consider [supporting us directly](https://buymeacoffee.com/tobydriscoll) or buying a print copy in [MATLAB](https://epubs.siam.org/doi/10.1137/1.9781611975086) or [Julia](https://epubs.siam.org/doi/10.1137/1.9781611977011) so that we can show them this kind of resource can be financially viable.

The table of contents is the same as for the print editions, and most of the content is the same as well. A few of the MATLAB functions in later chapters have been updated from the original for better clarity. Exercises are largely the same, but there have been updates, changes, and additions to those in print, so always check the online version for exercise numbers and details.

Chapters 1–6 can be used as a single-semester introduction to numerical methods, and Chapters 7–13 can be used for a follow-up course. Many sections have downstream connections, but if you are looking for sections that can be skipped with minimal impact on other sections, we suggest the following:

- [Section 3.4](#section-leastsq-house) on Householder QR factorization
- [Section 4.7](#section-nonlineqn-nlsq) on nonlinear least squares
- [Section 5.3](#section-localapprox-splines) on cubic splines
- [Section 5.7](#section-localapprox-adaptive) on adaptive integration
- [Section 6.8](#section-ivp-zerostability) on zero-stability of ODE methods
- [Section 7.5](#section-matrixanaly-dimreduce) on dimension reduction
- [Section 8.7](#section-krylov-matrixfree) on matrix-free linear algebra
- [Section 8.8](#section-krylov-precod) on preconditioning
- [Section 9.4](#section-globalapprox-orthogonal) on orthogonal polynomials
- [Section 9.7](#section-globalapprox-improper) on improper integrals
- [Section 10.6](#section-bvp-galerkin) on the Galerkin method
- [Section 12.4](#section-advection-wave) on the wave equation
