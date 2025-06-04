---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering: false
---
# How to use this book

## Everyone

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

If you hover over the code cell, you will see an icon in the upper right corner that allows you to copy that code to your clipboard. You will probably do this kind of thing a lot, if you work though this book.

One thing you do *not* have to copy is the code for the functions that are presented in the book as implementations of key algorithms. They will all be available once you install them, as described in the links above. (You didn't skip the "first, do this" part, did you? Not a great look for you.)

Here is a short overview of major features of the book.

![usage video](_static/FNC-usage.mp4)


## Instructors