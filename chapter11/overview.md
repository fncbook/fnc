# 11. Diffusion equations

```{index} Wedge Antilles, A New Hope
```

```{epigraph}
That's impossible, even for a computer.

— Wedge Antilles, *Star Wars: A New Hope* 
```

To this point we have considered only ordinary differential equations, those having only one independent variable. In the rest of this book, we introduce the huge topic of solving partial differential equations (PDEs). We begin by pairing time with space.

As we have seen with initial- and boundary-value problems, the crucial difference between time and space is that information can flow only forward in time. PDEs that include both time and space variables are sometimes referred to as initial-boundary-value problems (IBVPs), because they bear characteristics of both types. As with IVP methods, our IBVP methods advance a solution in time. Like BVP methods, the boundary values must be respected, and all values in space are represented simultaneously.

Finally, while we usually refer to "time" and "space," the independent variables can be more abstract. (One such example is introduced in the next section.) However, they are all time-like or space-like. Mathematical analysis presented in courses on PDEs can be used to determine definitively what type of PDE one is dealing with.

