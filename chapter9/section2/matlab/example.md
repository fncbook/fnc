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
[**Demo %s**](#demo-barycentric-example)

```{code-cell}
f = @(x) sin( exp(2 * x) );
clf,  fplot(f, [0, 1], displayname="function")
xlabel('x'),  ylabel('f(x)')   
legend(location="southwest");
```

We start with 4 equally spaced nodes ($n=3$).

```{code-cell}
t = linspace(0, 1, 4)'; 
y = f(t);
p = polyinterp(t, y);
hold on,  fplot(p, [0, 1], displayname="interpolant on 4 nodes")
scatter(t, y, 'k', displayname="nodes")
```

The curves always intersect at the interpolation nodes. For $n=6$, the interpolant is noticeably better.

```{code-cell}
cla,  fplot(f, [0, 1], displayname="function")
t = linspace(0, 1, 7)'; 
y = f(t);
p = polyinterp(t, y);
hold on,  fplot(p, [0, 1], displayname="interpolant on 7 nodes")
scatter(t, y, 'k', displayname="nodes")
```
