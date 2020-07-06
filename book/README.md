# Fundamentals of Numerical Computation

These are source files used to generate the online text [*Fundamentals of Numerical Computation*](https://fncbook.github.io/fnc) by T. A. Driscoll and R. J. Braun.

## Installing locally

Each chapter includes notebooks that were run with Julia 1.4. They should work with any later 1.x release as well (provided dependencies honor the standard versioning scheme). To get the book functions as well as all demos, first [install Julia](https://julialang.org/downloads/). Then execute

```julia
] add https://github.com/fncbook/fnc
```

After this, hit the backspace key to go back to the main Julia prompt. This process should only need to be done once per installation. In order to use the functions you must enter

```julia
using FundamentalsNumericalComputation
```

at the Julia prompt in each session of Julia. After this, all the text's functions can be accessed with the prefix `FNC`. E.g., `FNC.lufact`, `FNC.rk23`, etc.

### Startup speed

There are two sources of noticeable slowness:

1. Any `using` statement requires a compilation step whenever a package or its dependencies have been updated. As of July 2020, this can consume several minutes for `FNC` due mainly to its dependence on ['Plots'](http://docs.juliaplots.org/latest/) and [`DifferentialEquations`](https://docs.sciml.ai/latest/index.html). Even on subsequent invocations, `using` this package can take 10-30 seconds for loading.
2. The first plot created in a Julia session may take up to a minute, which is standard for `Plots`.

Both of these lag issues are well known to Julia developers and are the target of future improvements. In the meantime power users might investigate using [`PackageCompiler`](https://julialang.github.io/PackageCompiler.jl/dev/) to get vastly better performance on both issues.

## Tooling

The book was created using [jupyter-book](https://jupyterbook.org) and [VS Code](https://code.visualstudio.com/) with the [Julia extension](https://github.com/julia-vscode/julia-vscode). These tools are already excellent and continue to evolve rapidly.

## License

The content is licensed under CC Attribution-ShareAlike 4.0, while the code is under an MIT license. All material is copyrighted.
