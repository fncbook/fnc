---
kernelspec:
  display_name: MATLAB
  language: matlab
  name: jupyter_matlab_kernel
numbering:
  headings: false
---
```{code-cell}
:tags: [remove-cell]
cd  /Users/driscoll/Documents/GitHub/fnc/matlab
FNC_init
```
[**Demo %s**](#demo-power-one)

Here we choose a magic 5×5 matrix and a random 5-vector.

```{code-cell}
A = magic(5) / 65;
x = randn(5, 1);
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
for j = 1:8
    x = A * x;
end
[x, A * x]
```

This phenomenon seems to occur regardless of the starting vector.

```{code-cell}
x = randn(5, 1);
for j = 1:8
    x = A * x;
end
[x, A * x]
```
