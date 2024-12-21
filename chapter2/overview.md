(chapter-linsys)=
# Linear systems of equations

```{index} Han Solo, Star Wars: A New Hope
```

```{epigraph}
It's all a lot of simple tricks and nonsense.

— Han Solo, *Star Wars: A New Hope*
```

One of the most frequently encountered tasks in scientific computation is the solution of the linear system of equations $\mathbf{A} \mathbf{x}=\mathbf{b}$ for a given square matrix $\mathbf{A}$ and vector $\mathbf{b}$.  This problem can be solved in a finite number of steps, using an algorithm equivalent to Gaussian elimination. Describing the algorithm is mostly an exercise in organizing some linear algebra.

Analyzing the algorithm requires new tools. Because the computations will take place in floating point, we must first discuss a system for measuring the "size" of a perturbation to a vector or matrix data. Once that is understood, we find that the conditioning of the square linear system problem is quite straightforward to describe. Finally, we will see that the algorithm may change when certain things are known about the matrix $\mathbf{A}$.

