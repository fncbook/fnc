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
[**Demo %s**](#demo-rootproblem-bessel)


```{code-cell}
J3 = @(x) besselj(3,x);
fplot(J3, [0, 20])
grid on
xlabel('x'), ylabel('J_3(x)')  
title('Bessel function') 
```
From the graph we see roots near 6, 10, 13, 16, and 19. We use `nlsolve` from the `NLsolve` package to find these roots accurately. It uses vector variables, so we have to code accordingly.
```{tip}
:class: dropdown
Type `\omega` followed by <kbd>Tab</kbd> to get the character `ω`.
The argument `ftol=1e-14` below is called a **keyword argument**. Here it sets a goal for the maximum value of $|f(x)|$.
```

```{code-cell}
omega = [];
for guess = [6, 10, 13, 16, 19]
    omega = [omega; fzero(J3, guess)];
end
omega
```

```{code-cell}
table(omega, J3(omega), 'VariableNames', {'root estimate', 'function value'})
```

```{code-cell}
hold on
scatter(omega, J3(omega))
title('Bessel roots')    
```

If instead we seek values at which $J_3(x)=0.2$, then we must find roots of the function $J_3(x)-0.2$.

```{code-cell}
omega = [];
for guess = [3, 6, 10, 13]
    f = @(x) J3(x) - 0.2;
    omega = [omega; fzero(f, guess)];
end
scatter(omega, J3(omega), '<')
```
