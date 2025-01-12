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
[**Demo %s**](#demo-matrixfree-deblur)

We repeat the earlier process to blur an original image $\mathbf{X}$ to get $\mathbf{Z}$.

```{code-cell}
:tags: [hide-input]
load mandrill
[m, n] = size(X);
v = [1/4, 1/2, 1/4];
B = spdiags(v, -1:1, m, m);
C = spdiags(v, -1:1, n, n);
blur = @(X) B^12 * X * C^12;
```

```{code-cell}
Z = blur(X);
clf,  imshow(Z, [0, 255])
title(("Blurred image"));
```

Now we imagine that $\mathbf{X}$ is unknown and that we want to recover it from $\mathbf{Z}$. We first need functions that translate between vector and matrix representations.

```{code-cell}
vec = @(X) reshape(X,m*n,1);
unvec = @(x) reshape(x,m,n);
T = @(x) vec( blur(unvec(x)) );
```
The blurring operators are symmetric, so we apply `minres` to the composite blurring transformation `T`.

```{code-cell}
y = gmres(T, vec(Z), 50, 1e-5);
Y = unvec(y);

subplot(121)
imshow(X, [0, 255])
title("Original")
subplot(122)
imshow(Y, [0, 255])
title(("Deblurred"));
```
