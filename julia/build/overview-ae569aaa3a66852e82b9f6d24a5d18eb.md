(chapter-nonlineqn)=

# 4. Roots of nonlinear equations

```{index} C3PO, A New Hope
```

```{epigraph}
He says "I found her," and keeps repeating, "She's here."

— C3PO, *Star Wars: A New Hope*
```

In this chapter we extend from linear algebra to deal with *nonlinear* algebraic problems. This kind of problem arises when there is a parameter or variable that can be changed in order to satisfy a constraint or achieve some goal. We start with scalar functions of a single variable, then generalize to $n$ variables and $n$ nonlinear equations. Finally, we generalize the problem of linear least squares to situations with more nonlinear constraints to satisfy than there are variables. In every case the strategy used is one of the cornerstones of numerical computing: *replace a problem you can't solve with an approximate one that you can.* In the context of nonlinear algebraic problems, the particular tactic is to set up and solve a sequence of linear problems of the types covered in the two previous chapters.
