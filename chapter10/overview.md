(chapter-bvp)=
# Boundary-value problems

```{index} C3PO, A New Hope
```

```{epigraph}
Don't shoot! Don't shoot!

— C3PO, *Star Wars: A New Hope* 
```

In @chapter-ivp we examined how to solve initial-value problems (IVP) for ordinary differential equations (ODEs). In an IVP the supplemental conditions give complete information about the state of the system at one value of the independent variable. However, not all ODE problems come in this form. Instead we might have only partial information about the solution at two different points. Such ODE problems are called *boundary-value problems* (BVP).

We begin with a numerical method that attempts to treat a BVP using IVP methods. It's unstable, so we turn to numerical methods that find the solution everywhere simultaneously. One of these directly uses a discrete form of the ODE to express the computational problem, while another first restates the ODE in terms of integration. Both methods lead to algebraic systems of equations to be solved for the entire solution all at once, using the methods of @chapter-linsys and @chapter-nonlineqn.
