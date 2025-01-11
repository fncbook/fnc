---
kernelspec:
  display_name: Python 3
  language: python
  name: python3
numbering:
  headings: false
---
```{code-cell}
:tags: [remove-cell]
exec(open("../../../python/FNC_init.py").read())
```
[**Demo %s**](#demo-house-qr)


We will use Householder reflections to produce a QR factorization of a matrix.

```{code-cell}
A = 1.0 + floor(9 * random.rand(6,4))
m, n = A.shape
print(A)
```

Our first step is to introduce zeros below the diagonal in column 1 by using {eq}`hhvector` and {eq}`hhreflect`. 

```{code-cell}
z = A[:, 0]
v = z - norm(z) * hstack([1, zeros(m-1)])
P_1 = eye(m) - (2 / dot(v, v)) * outer(v, v)   # reflector
```

We check that this reflector introduces zeros as it should:

```{code-cell}
print(P_1 @ z)
```

Now we replace $\mathbf{A}$ by $\mathbf{P}_1\mathbf{A}$.

```{code-cell}
A = P_1 @ A
print(A)
```

We are set to put zeros into column 2. We must not use row 1 in any way, lest it destroy the zeros we just introduced. So we leave it out of the next reflector.

```{code-cell}
z = A[1:, 1]
v = z - norm(z) * hstack([1, zeros(m-2)])
P_2 = eye(m-1) - (2 / dot(v, v)) * outer(v, v) 
```

We now apply this reflector to rows 2 and below only.

```{code-cell}
A[1:, 1:] = P_2 @ A[1:, 1:]
print(A)
```

We need to iterate the process for the last two columns.

```{code-cell}
for j in [2, 3]:
    z = A[j:, j]
    v = z - norm(z) * hstack([1, zeros(m-j-1)])
    P = eye(m-j) - (2 / dot(v, v)) * outer(v, v)
    A[j:, j:] = P @ A[j:, j:]
```

We have now reduced the original to an upper triangular matrix using four orthogonal Householder reflections:

```{index} Python; triu
```

```{code-cell}
R = triu(A)
print(R)
```
