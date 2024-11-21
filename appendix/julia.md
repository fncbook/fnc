# Julia cheatsheet

## Ways to work in Julia

### REPL (Read-Eval-Print Loop)

This is the [default way to start and run](https://docs.julialang.org/en/v1/stdlib/REPL/) Julia. 

- Press `]` to enter package management mode; press <kbd>backspace</kbd> to leave it.
- Press <kbd>Ctrl</kbd>-<kbd>C</kbd> to interrupt execution.
- Press <kbd>Ctrl</kbd>-<kbd>D</kbd> to exit.
- Press `?` to enter help mode.
- Press Up/Down arrows to walk through the command history.
- There are other interactive utilities described [in the manual]](https://docs.julialang.org/en/v1/stdlib/InteractiveUtils/).

### Jupyter notebook

Install the `IJulia` package, then import it and run `IJulia.notebook()` at the REPL to launch a session in a web browser. 

- You can put multiple lines in a code cell. Press <kbd>Shift</kbd>-<kbd>Enter</kbd> to execute the cell.
- Click the Stop button (square icon) to interrupt execution.
- Use `include Pkg` and then `Pkg.add`, etc. for package management.

### Pluto notebook

Unlike a Jupyter notebook, a Pluto notebook is reactive in the sense that a change anywhere affects all of the code cells immediately. To use it, install and import the \texttt{Pluto} package, then use \texttt{Pluto.run()} to launch it in a web browser.

### Visual Studio Code
Install the Julia extension and run the REPL in a pane. This option is most like an IDE such as PyCharm or Spyder and has many features to assist with coding.

## Basics

|  Item/Reference     |   Example                 |
|:--------------------|:--------------------------|
| Help on a function  | `?func`                   |
| [Search help](https://docs.julialang.org/en/v1/stdlib/InteractiveUtils/#Base.Docs.apropos)         | `apropos("topic")`        |
| [Assign to variable](https://docs.julialang.org/en/v1/manual/variables/#man-variables)  | `=`                       |
| Comment             | starts with `#`           |
| [Function definition](https://docs.julialang.org/en/v1/manual/functions/#man-anonymous-functions) | `function foo(x,y)` or `foo = (x,y) -> ` |
| Other [defined symbols](https://docs.julialang.org/en/v1/base/punctuation/#Punctuation) |  `%` `÷` `∘`   
| Symbols using [LaTeX nomenclature](https://en.wikipedia.org/wiki/List_of_mathematical_symbols_by_subject) | `\beta`+<kbd>Tab</kbd>, `\le`+<kbd>Tab</kbd>, `x\_0`+<kbd>Tab</kbd>, `z\hat`+<kbd>Tab</kbd> |

## Input and output

|                     |    Example                |   Package      |
|:--------------------|:--------------------------|:---------------|
| [Integers](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/#Integers) |  `176`, `-1`, `0`  |    |
| [Floats](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/#Floating-Point-Numbers)  |  `3.14`, `1.0`, `0.`, `NaN`, `Inf`   |
| Scientific notation | `1.234e-5`                |                |
| Complex numbers     | `3+4im`, `complex(3,4)`   |                |
| Constants           | `pi`, `π`, `exp(1)`       |                |
| [Strings](https://docs.julialang.org/en/v1/manual/strings/)              | `"string"`                |                |
| Show result         | `@show`                   |                |
| Print to screen     | `print`, `println`         |                |
| [Formatted print](https://docs.julialang.org/en/v1/stdlib/Printf/)     | `@sprintf`, `@printf`      | [Printf](https://docs.julialang.org/en/v1/stdlib/Printf/)       |
| [Interpolate into string](https://docs.julialang.org/en/v1/manual/strings/#string-interpolation) | `"x is $x"`           |                |
| [Table](https://ronisbr.github.io/PrettyTables.jl/stable/)        | `pretty_table`             | [PrettyTables](https://ronisbr.github.io/PrettyTables.jl/stable/)  |

## Operators

|                     |     Example              |
|:--------------------|:--------------------------|
| Arithmetic          | `+` `-` `*` `/` `^`   |
| [Broadcast over array](https://docs.julialang.org/en/v1/manual/arrays/#Broadcasting)    | prefix with dot, or `@.` in front of all       |
| [Equality test](https://docs.julialang.org/en/v1/base/math/#Base.:==)       | `==`                       |
| Comparison          | `<` `<=` `>` `>=`  |
| Logical AND / OR / NOT   | `&`  `|`  `!`       |
| [Short-circuit logic](https://docs.julialang.org/en/v1/manual/control-flow/#Short-Circuit-Evaluation)  | `&&` `||`              |

## Errors
|  Message/Symptom          |       Possible cause            |
|:--------------------|:--------------------------|
| `BoundsError`       | Accessed a nonexistent array element  | 
| `MethodError`         | Probably used the wrong number/type of function arguments |
| [`InexactError`](https://docs.julialang.org/en/v1/manual/conversion-and-promotion/)       | Illegal type conversion, maybe by assigning to an array |
| "Cannot juxtapose string literal" |   Invalid string construction | 
| Unexpected [array changes](https://docs.julialang.org/en/v1/manual/arrays/#man-multi-dim-arrays) | Use `copy` instead of array reference  |
| Wrong [norm](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.norm)         |  Used `norm` when `opnorm` is needed   |

## Iteration

|                     |      Example           |
|:--------------------|:--------------------------|
| [Predetermined number](https://docs.julialang.org/en/v1/manual/control-flow/#man-loops)  | `for i in 1:10`         |
| [Based on condition](https://docs.julialang.org/en/v1/manual/control-flow/#man-loops)   | `while abs(x) > 1`       |
| [Comprehension](https://docs.julialang.org/en/v1/manual/arrays/#man-comprehensions)       |  `[i+j for i in 1:3, j in 1:3]` |
| [Generator](https://docs.julialang.org/en/v1/manual/arrays/#Generator-Expressions)           | `sum(k for k in 1:10)`    |

## Vectors, matrices, arrays

|                     |    Example                |   Package      |
|:--------------------|:--------------------------|:---------------|
| All ones, all zeros  | `ones(4)`, `zeros(2,5)`  |                |
| Random elements      | `rand(100)`, `randn(3,3)` |                |
| [Concatenate](https://docs.julialang.org/en/v1/manual/arrays/#man-array-concatenation)          | `[1;2;3]` vertically, `[1 2 3]` horizontally |           |
| Get dimensions       | `length`, `size`         |                |
| [Extract diagonal](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.diag)     | `diag(A)`, `diag(A,-1)`  |  LinearAlgebra   |
| [Build by diagonal](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.diagm)    | `diagm([1,2,3])`, `diagm(1=>[-1,-1])`  |   LinearAlgebra  |
| [Access element](https://docs.julialang.org/en/v1/manual/arrays/#man-array-indexing)       | `x[2]`,  `A[1,4]`      |                |
| [Access block](https://docs.julialang.org/en/v1/manual/arrays/#man-array-indexing)         | `x[1:end-1]`, `A[1:4,2:3]`  |                |
| [Access row/column](https://docs.julialang.org/en/v1/manual/arrays/#man-array-indexing)    | `A[4,:]`, `A[:,2:2:end]`   |                |
| [Range by step size](https://docs.julialang.org/en/v1/base/math/#Base.::)  | `0:4`, `1:2:9`,  `5:-1:1`   |                |
| [Range by length](https://docs.julialang.org/en/v1/base/math/#Base.range)      | `range(a,b,length=10)`   |                |
| Make a copy          | `copy`                  |                |


## Linear algebra

|                     |    Example                |   Package      |
|:--------------------|:--------------------------|:---------------|
| [Adjoint/transpose](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#Base.adjoint)    | `A'`                     |    |
| [Inner (scalar) product](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.dot) |  `dot(x,y)`            |  LinearAlgebra   |
| Outer product        | `x*v'`                   |    |
| [Solve linear system](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#man-linalg) | `A\b`                    |    |
| [Solve overdetermined system](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#man-linalg) | `A\b`             |    |
| [Identity](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.I)             | `I`, `I(5)`              |  LinearAlgebra   |
| [Factorizations](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#man-linalg-factorizations)       | `lu`, `qr`, `Cholesky`      |  LinearAlgebra   |
| [Norm](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.norm)                 | `norm`, `opnorm`         |  LinearAlgebra   |
| [Condition number](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.cond)     | `cond`                    |  LinearAlgebra   |
| [Eigenvalues](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#man-linalg)          | `eigen`, `eigvals`, `eigs`     |  LinearAlgebra   |
| [Singular values](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#man-linalg)      | `svd`, `svdvals`         |  LinearAlgebra   |

## Major problem types

|                     |    Example                |   Package      |
|:--------------------|:--------------------------|:---------------|
| [Linear system / least squares](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#man-linalg) |  `A\b`          |     |
| [Sparse matrix](https://docs.julialang.org/en/v1/stdlib/SparseArrays/#Sparse-Arrays)        | `A\b`, `gmres`, `minres`, `cg`, `eigs`  | SparseArrays, [IterativeSolvers](https://iterativesolvers.julialinearalgebra.org/stable/),  [Arpack](https://arpack.julialinearalgebra.org/stable/) |
| [Polynomial interpolation/approximation](https://juliamath.github.io/Polynomials.jl/stable/#Fitting-arbitrary-data-1)   | `fit`    | [Polynomials](https://juliamath.github.io/Polynomials.jl/stable/) |
| [Polynomial roots](https://juliamath.github.io/Polynomials.jl/stable/#Root-finding-1)   |  `roots(Polynomial(-1,1,1))`  | [Polynomials](https://juliamath.github.io/Polynomials.jl/stable/) |
| Rootfinding         | `nlsolve`                  | [NLsolve](https://github.com/JuliaNLSolvers/NLsolve.jl) |
| Integration        | `quadgk(sin,0,pi)`            | [QuadGK](https://juliamath.github.io/QuadGK.jl/stable/) |
| [Initial-value problem](https://diffeq.sciml.ai/latest/tutorials/ode_example/#ode_example)  | `ODEProblem`, `solve`  | [DifferentialEquations](https://diffeq.sciml.ai/latest/) | 
| [Boundary-value problem](https://diffeq.sciml.ai/latest/tutorials/bvp_example/#Boundary-Value-Problems) | `BVProblem`, `solve`  | [DifferentialEquations](https://diffeq.sciml.ai/latest/) |


## Graphics

|                     |     Example              |   Package      |
|:--------------------|:--------------------------|:---------------|
| [One variable](https://docs.juliaplots.org/latest/generated/gr/#gr-ref1)           | `plot`, `scatter`, `hline`, `vline`   |  [Plots](https://docs.juliaplots.org/latest/)  |
| [Two variables](https://docs.juliaplots.org/latest/generated/gr/#gr-ref22)          | `surface`, `contour`, `contourf`   |  [Plots](https://docs.juliaplots.org/latest/)  |
| [Network](https://docs.juliaplots.org/latest/graphrecipes/introduction/)               | `graphplot`            |   GraphRecipes
| [Matrix values](https://docs.juliaplots.org/latest/generated/gr/#gr-ref32)         | `spy`                 |  [Plots](https://docs.juliaplots.org/latest/)  |
| Modify existing plot  | `plot!` etc., `xlabel!`, `xlims!`, `title!`      |  [Plots](https://docs.juliaplots.org/latest/)  |
| [Layouts](https://docs.juliaplots.org/latest/layouts/#layouts)               | `layout=(2,1)`, `subplot=2`   |  [Plots](https://docs.juliaplots.org/latest/)  |
| Log scales            | `yscale=:log10`    |   [Plots](https://docs.juliaplots.org/latest/)  |
| Colors and images        | `RGB(1,0,1)`, `Gray.(img)` |     [Images](https://juliaimages.org/stable/)  |

