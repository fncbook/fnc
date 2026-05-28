---
numbering: false
---
(section-setup-julia)=

# Setting up Julia for this book

Julia, and all the packages this book depends on, is free and open-source. Much of the functionality outside the core is distributed via packages that need to be installed once per system.

## Installing Julia

To install Julia, go to the [Julia website](https://julialang.org/downloads/) and follow the instructions there to download the appropriate version for your operating system. Any 1.x version at 1.11 or greater should be fine. I strongly recommend using the `juliaup` version manager.

The most basic way to use Julia is to double-click on the app icon, or to open a Terminal/Command window and enter `julia`. Either way, you should end up with a console that has a `julia>` prompt, where you can enter commands and get results back in a read–eval–print loop (REPL). But read on...

## Coding environments

You *could* interact with Julia only by typing in at the REPL and then pasting the results into a word processor, but you can do much better. Two of the most popular ways to use Julia are:

- **[Jupyter lab](https://jupyter.org)**. This is a notebook-based interface that mixes cells having text and code, including text and graphical output. You write and execute code within your web browser, seeing results as they are computed. Code is run and saved on your own computer. If you have ever installed Python on your computer, you may have a working copy of Jupyter already.
- **[VS Code](https://code.visualstudio.com)**. This is a free, full-featured code editor that can be extended by installing free [Julia-specific tools](https://code.visualstudio.com/docs/languages/julia). This gives you a "Julia: Start REPL" command, or you can write and run Jupyter notebooks within the app. VS Code has tight integration with Microsoft's Copilot AI. This book (the version you are reading now, anyway) was written in VS Code.

## Installing the book's functions

Most Julia functionality is contained in packages. A few of these form the *standard library* that is installed with Julia itself. The rest are downloaded and installed by Julia's package manager, `Pkg`.

This book's functions are distributed in a package called `FNCFunctions`. Running all the demos and code examples requires installing a number of other independently maintained packages as well. To do all the installations in one shot, open a Julia prompt, then copy and execute the following code block there:

```julia
import Pkg
Pkg.add(
    [
    "FNCFunctions", "Arpack", "BoundaryValueDiffEq", "Dierckx", "FFTW", "GraphRecipes", "IJulia", "Images", "IncompleteLU", "IterativeSolvers", "JLD2", "LaTeXStrings", "LinearMaps", "MAT", "MatrixDepot", "NLsolve", 
    "OrdinaryDiffEq", "OrdinaryDiffEqLowOrderRK", "OrdinaryDiffEqRosenbrock", "Plots", "Polynomials", 
    "Preconditioners", "PrettyTables", "QuadGK", "SpecialFunctions", "TestImages"
    ]
)
Pkg.precompile() 
Pkg.build("IJulia")
```

This will take several minutes to finish. Once it's complete, you won't need to repeat this step. In addition to installing the packages needed to reproduce the examples, the code above also makes Julia available to run via Jupyter notebooks, if Jupyter was already installed on your computer.

If this step doesn't work, see @section-packages.

## Using packages

Once a package is installed (or is part of the standard library), it has to be loaded each time Julia is started in order to make its functions available. Julia offers both `import` and `using` as ways to load packages. The difference is that `import` loads a package into its own namespace, so that you have to refer to its functions with a prefix. For example:

``` julia-repl
julia> import Statistics

julia> mean
ERROR: UndefVarError: `mean` not defined in `Main`
Suggestion: check for spelling errors or missing imports.
Hint: a global variable of this name also exists in Statistics.

julia> Statistics.mean
mean (generic function with 6 methods)
```

If you load a package with `using` rather than `import`, then the package may also make some of its functions available in the global namespace. For example:

``` julia-repl
julia> using Statistics

julia> mean
mean (generic function with 6 methods)
```

This is convenient for functions that you will call frequently.

The package `FNCFunctions` written for this book does *not* put functions into the global namespace. This is a deliberate reminder that they are for learning purposes only and not meant as tools for serious work. The package does define an alias called `FNC`, though, so you can write, for instance, `FNC.lufact` instead of `FNCFunctions.lufact`.

(section-packages)=

## Packages used in the book

In order to avoid repeating low-information code, the book demos are run with a few packages installed and always imported:

- `FNCFunctions`
- `Printf`, `LinearAlgebra` (part of the default Julia installation)
- `Plots`, `PrettyTables`, `LaTeXStrings` (external packages)

Of the external packages, `Plots` is essential for the exercises, while the others are used to make the results look nicer in the book.

Throughout the book demos there are other external packages loaded as needed. These include:

- `Arpack`
- `BoundaryValueDiffEq`
- `Dierckx`
- `FFTW`
- `GraphRecipes`
- `Images`
- `IncompleteLU`
- `IterativeSolvers`
- `JLD2`
- `LinearMaps`
- `MatrixDepot`
- `MAT`
- `NLsolve`
- `OrdinaryDiffEq`
- `OrdinaryDiffEqLowOrderRK`
- `OrdinaryDiffEqRosenbrock`
- `Polynomials`
- `Preconditioners`
- `QuadGK`
- `SpecialFunctions`
- `TestImages`

It's possible that future versions of these packages will cause errors in some of this book's examples, which is not an uncommon frustration in software. To deal with such situations, Julia's package manager makes it possible to recreate the same versions of all packages used to produce this book. It requires just a few extra steps for you.

1. Create a folder (directory) on your computer. I'll call it `FNC-class` here to be concrete. You will need to know the entire system path to that folder; for demonstration here I will use `/path/to/folder/FNC-class`.
2. Download the file [Project.toml](/Project.toml) and the file [Manifest-v1.12.toml](/Manifest-v1.12.toml), and save both into the `FNC-class` folder with those exact file names.
3. Install version 1.12 of Julia. (Any 1.12.x should work.)
4. Open Julia in a terminal window or a blank Jupyter notebook. Enter the following:

```julia
import Pkg
Pkg.activate("/path/to/folder/FNC-class")
Pkg.instantiate()
Pkg.precompile()
```

This should download and install everything in the preserved state.
5. Each time you start Julia, you need to run the `import Pkg` and `Pkg.activate` lines above. If you forget, then the packages won't be available to load until you do.
