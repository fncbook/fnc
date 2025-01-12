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
[**Demo %s**](#demo-norms-vector)

```{index} ! MATLAB; norm
```

```{code-cell}
x = [2; -3; 1; -1];
twonorm = norm(x)    % or norm(x, 2)
```

```{code-cell}
infnorm = norm(x, Inf)
```

```{code-cell}
onenorm = norm(x, 1)
```
