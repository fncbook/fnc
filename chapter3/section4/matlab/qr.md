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
[**Demo %s**](#demo-house-qr)


We will use Householder reflections to produce a QR factorization of a matrix.

```{code-cell}
A = magic(6);
A = A(:, 1:4);
[m, n] = size(A)
```

```{index} ! MATLAB; eye
```

Our first step is to introduce zeros below the diagonal in column 1 by using {eq}`hhvector` and {eq}`hhreflect`. 

```{code-cell}
z = A(:, 1);
v = z - norm(z) * eye(m,1);
P_1 = eye(m) - 2 / (v' * v) * (v * v'); 
```

We check that this reflector introduces zeros as it should:

```{code-cell}
P_1 * z
```

Now we replace $\mathbf{A}$ by $\mathbf{P}_1\mathbf{A}$.

```{code-cell}
A = P_1 * A
```

We are set to put zeros into column 2. We must not use row 1 in any way, lest it destroy the zeros we just introduced. So we leave it out of the next reflector.

```{code-cell}
z = A(2:m, 2);
v = z - norm(z) * eye(m-1, 1);
P_2 = eye(m-1) - 2 / (v' * v) * (v * v');
```

We now apply this reflector to rows 2 and below only.

```{code-cell}
A(2:m, 2:n) = P_2 * A(2:m, 2:n)
```

We need to iterate the process for the last two columns.

```{code-cell}
for j = 3:n
    z = A(j:m,j);
    k = m-j+1;
    v = z - norm(z) * eye(k, 1);
    P = eye(k) - 2 / (v' * v) * (v * v');
    A(j:m, j:n) = P * A(j:m, j:n);
end
```

We have now reduced the original to an upper triangular matrix using four orthogonal Householder reflections:

```{code-cell}
R = A
```
