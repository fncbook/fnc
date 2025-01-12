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
[**Demo %s**](#demo-stability-quadbad)


```{index} ! MATLAB; scientific notation
```

We apply the quadratic formula to find the roots of a quadratic via {eq}`quadunstable`. 

```{tip}
:class: dropdown
A number in scientific notation is entered as `1.23e4` rather than as `1.23 * 10^4`.
```

```{code-cell}
a = 1;  b = -(1e6 + 1e-6);  c = 1;
x1 = (-b + sqrt(b^2 - 4*a*c)) / (2*a)
x2 = (-b - sqrt(b^2 - 4*a*c)) / (2*a)
```

The first value is correct to all stored digits, but the second has fewer than six accurate digits:

```{code-cell}
error = abs(x2 - 1e-6) / 1e-6 
accurate_digits = -log10(error)
```

 The instability is easily explained. Since $a=c=1$, we treat them as exact numbers. First, we compute the condition numbers with respect to $b$ for each elementary step in finding the "good" root:

| Calculation | Result | $\kappa$ |
|:------------|:-------|:---------|
|$u_1 = b^2$  | $1.000000000002000\times 10^{12}$ |  2 |
|$u_2 = u_1 - 4$ | $9.999999999980000\times 10^{11}$  | $\approx 1.00$ |
|$u_3 = \sqrt{u_2}$ | $999999.9999990000$ | 1/2 |
|$u_4 = u_3 - b$ | $2000000$ | $\approx 0.500$ |
|$u_5 = u_4/2$ | $1000000$  | 1 |

Using {eq}`condition-chain`, the chain rule for condition numbers, the conditioning of the entire chain is the product of the individual steps, so there is essentially no growth of relative error here. However, if we use the quadratic formula for the "bad" root, the next-to-last step becomes $u_4=(-u_3) - b$, and now  $\kappa=|u_3|/|u_4|\approx 5\times 10^{11}$. So we can expect to lose 11 digits of accuracy, which is what we observed. The key issue is the subtractive cancellation in this one step.
