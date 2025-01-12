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
[**Demo %s**](#demo-matrixfree-blur)

We use a readily available test image.

```{code-cell}
load mandrill
[m, n] = size(X);
clf
imshow(X, [0, 255])
title('Original image')    % ignore this 
```

We define the one-dimensional tridiagonal blurring matrices.

```{code-cell}
v = [1/4, 1/2, 1/4];
B = spdiags(v, -1:1, m, m);
C = spdiags(v, -1:1, n, n);
```

Finally, we show the results of using $k=12$ repetitions of the blur in each direction.

```{code-cell}
blur = @(X) B^12 * X * C^12;
imshow(blur(X), [0, 255])
title(('Blurred image'));
```
