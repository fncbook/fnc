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
[**Demo %s**](#demo-int-antideriv)

The antiderivative of $e^x$ is, of course, itself. That makes evaluation of $\int_0^1 e^x\,dx$ by the Fundamental Theorem trivial.

```{code-cell}
format long
exact = exp(1) - 1
```

```{index} ! MATLAB; integral
```

MATLAB has numerical integrator `integral` that estimates the value without finding the antiderivative first. As you can see here, it can be as accurate as floating-point precision allows.

```{code-cell}
integral(@(x) exp(x), 0, 1)
```

The numerical approach is also far more robust. For example, $e^{\,\sin x}$ has no useful antiderivative. But numerically, it's no more difficult.

```{code-cell}
integral(@(x) exp(sin(x)), 0, 1)
```

When you look at the graphs of these functions, what's remarkable is that one of these areas is basic calculus while the other is almost impenetrable analytically. From a numerical standpoint, they are practically the same problem.

```{code-cell}
:tags: [hide-input]
x = linspace(0, 1, 201)';
subplot(2,1,1), fill([x; 1; 0], [exp(x); 0;0 ], [1, 0.9, 0.9])
title('exp(x)')
ylabel('f(x)')
subplot(2, 1, 2), fill([x; 1; 0], [exp(sin(x)); 0; 0], [1, 0.9, 0.9])
title('exp(sin(x))')
xlabel('x'), ylabel(('f(x)'));
```
