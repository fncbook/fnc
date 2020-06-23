# Introduction

```{epigraph}
You must unlearn what you have learned. 

-- Yoda, _The Empire Strikes Back_
```

Our first step is to discretize the real numbers---specifically, to replace them with a finite surrogate set of numbers. This step keeps the time and storage requirements for operating with each number at constant levels, but (almost) every data set and arithmetic operation is perturbed slightly away from its idealized mathematical value. We can easily keep the individual roundoff errors very small, so small that simple random accumulation is unlikely to bother us. However, some problems are extremely sensitive to these perturbations, a trait we quantify using a *condition number*. Problems with large condition numbers are difficult to solve accurately. Furthermore, even when the condition number of a problem is not large, some algorithms for solving it allow errors to grow enormously. We call these algorithms *unstable*. In this chapter we discuss these ideas in simple settings before moving on to the more realistic problems in the rest of the book.

## Important terms

```{glossary}
algorithm
  Set of instructions for transforming data into a result.

backward error
  Change in the input required to produce the result found by an inexact algorithm.

condition number
  Ratio of the size of change in the output of a function to the size of change in the infinitesimal input that produced it.

double precision
  Typical standard in floating-point representation, using 64 bits to achieve about 16 decimal significant digits of precision.

floating point numbers
  A finite set that substitutes for the real numbers in machine calculations. Denoted by $\mathbb{F}$.

ill-conditioned
  Exhibiting a large condition number, indicating high sensitivity of a result to changes in the data.

machine epsilon
  Distance from 1 to the next-largest floating point number. Also called unit roundoff or machine precision, though the usages are not completely consistent across different references.

subtractive cancellation
  Growth in relative error that occurs when two numbers are added/subtracted to get a result that is much smaller in magnitude than the operands.

unstable
  Allowing perturbations of the data to have much larger effects on the results than can be explained by the problem's condition number.
```

## Important Julia commands and keywords

```{glossary}
`function`
  Keyword used to start a function definition.

`NaN`
  Not a Number, the result of an undefined arithmetic operation.
```
