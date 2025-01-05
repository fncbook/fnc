---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
Here is a straightforward implementation of matrix-vector multiplication.

```{code-cell}
n = 6
A = randn(n, n)
x = rand(n)
y = zeros(n)
for i in 1:n
    for j in 1:n
        y[i] += A[i, j] * x[j]    # 1 multiply, 1 add
    end
end
```

Each of the loops implies a summation of flops. The total flop count for this algorithm is

$$ \sum_{i=1}^n \sum_{j=1}^n 2 = \sum_{i=1}^n 2n = 2n^2. $$

Since the matrix $\mathbf{A}$ has $n^2$ elements, all of which have to be involved in the product, it seems unlikely that we could get a flop count that is smaller than $O(n^2)$ in general.

```{index} ! Julia; push\!, ! Julia; for
```

Let's run an experiment with the built-in matrix-vector multiplication. Note that Julia is unusual in that loops have a variable scope separate from its enclosing code. Thus, `for n in n` below means that inside the loop, the name `n` will take on each one of the values that were previously assigned to the vector `n`.
```{tip}
The `push!` function attaches a new value to the end of a vector.
```

```{code-cell}
n = 1000:1000:5000
t = []
for n in n
    A = randn(n, n)  
    x = randn(n)
    time = @elapsed for j in 1:80; A * x; end
    push!(t, time)
end
```

The reason for doing multiple repetitions at each value of $n$ in the loop above is to avoid having times so short that the resolution of the timer is significant.

```{code-cell}
@pt :header = ["size", "time (sec.)"] [n t]
```

```{index} Julia; Boolean indexing
```

Looking at the timings just for $n=2000$ and $n=4000$, they have ratio
```{tip}
The expression `n.==4000` here produces a vector of Boolean (true/false) values the same size as `n`. This result is used to index within `t`, accessing only the value for which the comparison is true.
```

```{code-cell}
@show t[n.==4000] ./ t[n.==2000];
```

If the run time is dominated by flops, then we expect this ratio to be

$$
\frac{2(4000)^2}{2(2000)^2}=4.
$$
