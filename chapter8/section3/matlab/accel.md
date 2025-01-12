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
[**Demo %s**](#demo-inviter-accel)

```{code-cell}
ev = [1, -0.75, 0.6, -0.4, 0];
A = triu(ones(5, 5), 1) + diag(ev);
```

We begin with a shift $s=0.7$, which is closest to the eigenvalue 0.6.

```{code-cell}
s = 0.7;
x = ones(5, 1);
y = (A - s * eye(5)) \ x;
beta = x(1) / y(1) + s
```

Note that the result is not yet any closer to the targeted 0.6. But we proceed (without being too picky about normalization here).

```{code-cell}
s = beta;
x = y / y(1);
y = (A - s * eye(5)) \ x;
beta = x(1) / y(1) + s
```

Still not much apparent progress. However, in just a few more iterations the results are dramatically better.

```{code-cell}
format long
for k = 1:4
    s = beta;
    x = y / y(1);
    y = (A - s * eye(5)) \ x;
    beta = x(1) / y(1) + s
end
```
