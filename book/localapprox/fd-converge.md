# Convergence of finite differences

```{index} finite differences
```

All of the finite difference formulas in the previous section based on equally spaced nodes converge as the node spacing $h$ decreases to zero. However, note that to discretize a function over an interval $[a,b]$, we use $h=(b-a)/n$, which implies $n=(b-a)/h=O(h^{-1})$. As $h\to 0$, the total number of nodes needed grows without bound. So we would like to make $h$ as large as possible while still achieving some acceptable accuracy.

```{index} truncation error; of a finite difference formula
```

To measure this kind of performance, we introduce the {term}`truncation error` of a finite difference formula, defined as the difference between the two sides of {eq}`fdxformula`:

```{math}
:label: truncFD
\tau_f(h) = f'(0) - \frac{1}{h} \sum_{k=-p}^{q} a_k f(kh),
```

where we have used translation invariance to set $x=0$ for simplicity. In order to make this expression useful, we expand it in a Taylor series about $h=0$ and ignore all but the leading terms.[^trunc]

[^trunc]: The term *truncation error* is derived from the idea that the finite difference formula, being finite, has to truncate the series representation and thus cannot be exactly correct for all functions.

````{proof:example}
The finite difference formula {eq}`forwardFD11` implies

```{math}
:label: fd1trunc
\begin{split}
  \tau_f(h) &= f'(0) - \frac{ f(h)-f(0)}{h} \\
  &=f'(0) - h^{-1} \left[ \bigl( f(0) + h f'(0) + \tfrac{1}{2}h^2f''(0)+ \cdots \bigr) - f(0) \right] \\
  & = -\frac{1}{2}h f''(0) + O(h^2).
  \end{split}
```
````

````{sidebar} Demo
:class: demo
{doc}`demos/fdconverge-1`
````

The most important conclusion of {eq}`fd1trunc` is that $\tau_f(h)=O(h)$. The dependence on $h$ raised to the first power leads us to call this a **first-order accurate** formula. In a first-order formula, cutting $h$ in half should reduce the error in the $f'$ estimate by about half as well, when $h$ is small enough.

## Higher-order accuracy

```{index} order of accuracy; of a finite difference formula
```

As a rule, including more function values in a finite difference formula (i.e., increasing $p$ and $q$ in {eq}`fdxformula`) gives a truncation error that depends on a higher power of $h$ and thus vanishes more quickly as $h\to 0$. The power of $h$ in the leading term of the truncation error is known as the {term}`order of accuracy`.

````{proof:example}
We compute the truncation error of {eq}`centerFD12`:
  
```{math}
:label: fd2trunc
\begin{split}
  \tau_f(h) &= f'(0) - \frac{ f(h)-f(-h)}{2h}\\
  &= f'(0) - (2h)^{-1} \left[ \bigl( f(0) + h f'(0) + \tfrac{1}{2}h^2f''(0)+ \tfrac{1}{6}h^3f'''(0)+ O(h^4) \bigr) \right.\\
  &\qquad - \left.  \bigl( f(0) - h f'(0) + \tfrac{1}{2}h^2f''(0) - \tfrac{1}{6}h^3f'''(0)+O(h^4) \bigr) \right] \\
  &= (2h)^{-1} \left[ \tfrac{1}{3}h^3f'''(0) + O(h^4) \right] = O(h^2).
\end{split}
```

Because the lowest-order term is proportional to $h^2$, we say that the method has order of accuracy equal to 2.
````

````{sidebar} Demo
:class: demo
{doc}`demos/fdconverge-2`
````

Equation {eq}`fd2trunc` implies that {eq}`centerFD12` is second-order accurate. Reducing $h$ by a factor of two should cut the error in the $f'$ estimate by a factor of roughly four. This convergence is much more rapid than for the first-order formula {eq}`forwardFD11`.

## Stability

The truncation error $\tau_f(h)$ of a finite difference formula is dominated by a leading term $O(h^m)$ for an integer $m$. This error decreases as $h\to 0$. However, we have not yet accounted for the effects of roundoff error. To keep matters as simple as possible, let's consider the forward difference

```{math}
\delta(h) = \frac{f(x+h)-f(x)}{h}.
```

As $h\to 0$, the numerator approaches zero even though the values $f(x+h)$ and $f(x)$ are not necessarily near zero. This is the recipe for subtractive cancellation error! In fact, finite difference formulas are inherently ill-conditioned as $h\to 0$ when evaluated in floating point arithmetic.

To be precise, recall that the condition number for the problem of computing $f(x+h)-f(x)$ is

```{math}
\kappa(h) = \frac{ \max\{|f(x+h)|,|f(x)|\} }{ |f(x+h)-f(x) | },
```

implying a relative error of size $\kappa(h) \epsilon_\text{mach}$ in its computation. Hence the numerical value we actually compute for $\delta$ is

```{math}
\begin{split}
\tilde{\delta}(h) &= \frac{f(x+h)-f(x)}{h}\, (1+\kappa(h)\epsilon_\text{mach}) \\
&= \delta(h) + \frac{ \max\{|f(x+h)|,|f(x)|\} }{ |f(x+h)-f(x) | }\cdot \frac{f(x+h)-f(x)}{h} \cdot \epsilon_\text{mach}\\
&= \delta(h) \pm  \frac{ \max\{|f(x+h)|,|f(x)|\} }{ h}\,\epsilon_\text{mach}.
\end{split}
```

As $h\to 0$, $f(x+h)=f(x)+O(h)$, so we can simplify the maximization in the expression to get

```{math}
\bigl| \tilde{\delta}(h) - \delta(h) \bigr| = |f(x)| \epsilon_\text{mach} h^{-1} + \epsilon_\text{mach} \cdot O(1).
```

Hence as $h \to 0$, the roundoff error is $O(h^{-1})$, which grows without bound. Combining truncation error with roundoff error leads to

```{math}
:label: FDround
\bigl|  f'(x) - \tilde{\delta}(h) \bigr| \le \bigl| \tau_f(h) \bigr| + \bigl|f(x) \bigr|\, \epsilon_\text{mach} h^{-1}.
```

```{index} subtractive cancellation
```

Equation {eq}`FDround` quantifies the contributions of both truncation and roundoff errors. While the truncation error $\tau$ vanishes as $h$ decreases, the roundoff error actually *increases* thanks to the subtractive cancellation. At some value of $h$ the two error contributions will be of roughly equal size. This occurs roughly when

```{math}
\bigl|f(x)\bigr|\, \epsilon_\text{mach} h^{-1} \approx C h, \quad \text{or} \quad h \approx K \sqrt{\epsilon_\text{mach}}
```

for a constant $K$ that depends on $x$ and $f$, but not $h$. For a method of truncation order $m$, the details of the subtractive cancellation are a bit different, but the conclusion generalizes to

```{math}
:label: FDtruncbalance
\epsilon_\text{mach} h^{-1} \approx C h^m, \quad \text{or} \quad h \approx K \epsilon_\text{mach}^{1/(m+1)}.
```

Finally, at the optimal $h$ from {eq}`FDtruncbalance`, both the truncation and roundoff errors are

```{math}
:label: FDoptimal
O(h^m) = O\left( \epsilon_\text{mach}^{ m/(m+1) } \right).
```

````{sidebar} Demo
:class: demo
{doc}`demos/fdconverge-round`
````

Hence for a first order formula ($m=1$), we can expect only $O\left(\sqrt{\epsilon_\text{mach}}\right)$ error, or about half of the number of accurate machine digits. As $m$ increases, we get ever closer to using the full accuracy available.

```{margin}
Higher-order finite difference methods are both more efficient and less vulnerable to roundoff than low-order methods.
```

The observations in {doc}`demos/fdconverge-round` match the analysis above quite well, with the optimal $h$ becoming larger and the optimal error getting smaller as the order increases. Thus, higher-order finite difference methods are both more efficient and less vulnerable to roundoff than low-order methods.

## Exercises

1. ⌨ Write a code to evaluate the centered 2nd-order finite difference approximation to $f'(\pi/7)$ for $f(x)=\cos(x)$ and $h=2^{-1},2^{-2},\ldots,2^{-7}$. On a log--log graph, plot the error as a function of  $h$ and compare it graphically to second-order convergence.

    ````{only} solutions
    ````

2. ✍ Derive the first two nonzero terms of the Taylor series at $h=0$ of the truncation error $\tau_{f}(h)$ for the formula {eq}`backwardFD11`.

    ````{only} solutions
    % $$\frac{1}{2}h f''(0) - \frac{1}{6}h^2 f'''(0)$$
    ````

3. ✍ Calculate the first nonzero term in the Taylor series of the truncation error $\tau_{f}(h)$ for the finite difference formula defined by the second row of~{ref}`table-FDforward`.

    ````{only} solutions
    % In[1]:= Series[f'[0]-(\[Minus]3/2 f [0]+2 f [h]\[Minus]1/2 f [2h])/h,{h,0,3}]
    % Out[1]= 1/3 (f^(3))[0] h^2+1/4 (f^(4))[0] h^3+O[h]^4
    ````

4. ✍ Calculate the first nonzero term in the Taylor series of the truncation error $\tau_{f}(h)$ for the finite difference formula defined by the third row of~{ref}`table-FDforward`.

    ````{only} solutions
    ````

5. ✍ Using natural extensions of our definitions, show that the formula {eq}`centerFD22` is second-order accurate. 

    ````{only} solutions
    ````

    (problem-fd-muc)=
6. ✍  A different way to derive finite difference formulas is the **method of undetermined coefficients**. Starting from {eq}`fdxformula`,

    ```{math}
    f'(x) \approx \frac{1}{h}\sum_{k=-p}^q a_k f(x+kh),
    ```

    let each $f(x+k h)$ be expanded in a series around $h=0$. When the coefficients of powers of $h$ are collected, one obtains

    ```{math}
    \frac{1}{h} \sum_{k=-p}^q a_k f(x+kh) = \frac{b_0}{h} + b_1 f'(x) + b_2 f''(x)h + \cdots,
    ```

    where

    ```{math}
    b_i = \sum_{k=-p}^q k^i a_k.
    ```

    In order to make the result as close as possible to $f'(x)$, we impose the conditions 

    ```{math}
    b_0 = 0,\, b_1=1,\, b_2=0,\, b_3=0,\,\ldots,\,b_{p+q}=0.
    ```

    This provides a system of linear equations for the weights.

    **(a)** For $p=q=2$, write out the system of equations for $a_{-2}$, $a_{-1}$, $a_0$, $a_1$, $a_2$.

    **(b)** Verify that the coefficients from the appropriate row of {numref}`table-FDcenter` satisfy the equations you wrote down in part (a).

    **(c)** Derive the finite difference formula for $p=1$, $q=2$ using the method of undetermined coefficients.

    ````{only} solutions
    ````
