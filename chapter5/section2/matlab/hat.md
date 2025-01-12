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
[**Demo %s**](#demo-pwlin-hat)

Let's define a set of four nodes (i.e., $n=3$ in our formulas).

```{index} ! Julia; annotate!
```

```{code-cell}
t = [0, 0.55, 0.7, 1];
```

We plot the hat functions $H_0,\ldots,H_3$.

```{code-cell}
clf
for k = 0:3
    subplot(4, 1, k+1)
    Hk = hatfun(t, k);
    fplot(Hk, [0, 1])
    hold on
    scatter(t, Hk(t))
    text(t(k+1), 0.6, sprintf("H_%d", k))
end
```
