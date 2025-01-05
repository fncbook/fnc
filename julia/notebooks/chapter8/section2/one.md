---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
Here we choose a random 5×5 matrix and a random 5-vector.

```{code-cell}
A = rand(1.0:9.0, 5, 5)
A = A ./ sum(A, dims=1)
x = randn(5)
```

Applying matrix-vector multiplication once doesn't do anything recognizable.

```{code-cell}
y = A * x
```

Repeating the multiplication still doesn't do anything obvious.

```{code-cell}
z = A * y
```

But if we keep repeating the matrix-vector multiplication, something remarkable happens: $\mathbf{A} \mathbf{x} \approx \mathbf{x}$.

```{code-cell}
for j in 1:8
    x = A * x
end
[x A * x]
```

This phenomenon seems to occur regardless of the starting vector.

```{code-cell}
x = randn(5)
for j in 1:8
    x = A * x
end
[x A * x]
```
