---
numbering:
  enumerator: 1.2.%s
---
(section-intro-conditioning)=
# Problems and conditioning

Let's think a bit about what must be the easiest math problem you've dealt with in quite some time: adding 1 to a number. Formally, we describe this problem as a function $f(x)=x+1$, where $x$ is any real number. 

On a computer, $x$ will be represented by its floating-point counterpart, $\fl(x)$. Given the property {eq}`fpbound`, we have $\fl(x)=x(1+\epsilon)$ for some $\epsilon$ satisfying $|\epsilon| < \macheps/2$.  There is no error in representing the value 1. 

Let's suppose that we are fortunate and that the addition proceeds exactly, with no additional errors. Then the machine result is just

```{math}
:label: plus1fp
  {y} = x(1+\epsilon)+1.
```

We can derive the relative error in this result:

```{math}
:label: plus1relerr
  \frac{ |{y}-f(x)| }{ |f(x)| } = \frac{ |(x+\epsilon x+1) - (x+1)| }{ |x+1| }
  = \frac{ |\epsilon x| }{ |x+1| } .
```

This error could be quite large if the denominator is small. In fact, we can make the relative error as large as we please by taking $x$ very close to $-1$. This is essentially what happened in {numref}`Demo %s <demo-float-arithmetic>`.

```{index} ! subtractive cancellation
```

You may have encountered this situation before when using significant digits for scientific calculations. Suppose we round all results to five decimal digits, and we add $-1.0012$ to $1.0000$. The result is $-0.0012$, or $-1.2\times 10^{-3}$ in scientific notation. Notice that even though both operands are specified to five digits, it makes no sense to write more than two digits in the answer because there is no information in the problem beyond their decimal places. 

This phenomenon is known as **subtractive cancellation**, or loss of significance. We may say that three digits were "lost" in the mapping from $-1.0012$ to $-0.0012$. There's no way the loss could be avoided, *regardless of the algorithm*, once we decided to round off everything to a fixed number of digits.

:::{prf:observation} Subtractive cancellation
Subtractive cancellation is a loss of accuracy that occurs when two numbers add or subtract to give a result that is much smaller in magnitude. It is one of the most common mechanisms introducing dramatic growth of errors in floating-point computation.
:::

In double precision, all  the values are represented to about 16 significant decimal digits of precision, but it's understood that subtractive cancellation may render some of those digits essentially meaningless.

## Condition numbers

Now we consider problems more generally. As above, we represent a problem as a function $f$ that maps a real data value $x$ to a real result $f(x)$. We abbreviate this situation by the notation $f:\real \mapsto \real$, where $\real$ represents the real number set.

When the problem $f$ is approximated in floating point on a computer, the data $x$ is represented as a floating-point value $\tilde{x}=\fl(x)$. Ignoring all other sources of error, we define the quantitative measure

```{math}
:label: condition1
   \frac{ \vphantom{\dfrac{\bigl|}{\bigl|}}\dfrac{|f(x)-f(\tilde{x})|}{|f(x)|} }{%
     \vphantom{\dfrac{\bigl|}{\bigl|}}\dfrac{|x-\tilde{x}|}{|x|} },
```

which is the ratio of the relative changes in result and data. We make this expression more convenient if we recall that floating-point arithmetic gives $\tilde{x}=x(1+\epsilon)$ for some value $|\epsilon|\le \macheps/2$. Hence

```{math}
:label: condition2
   \dfrac{\left|f(x)-f(x+\epsilon x)\right| } {|\epsilon f(x)|}.
```

Finally, we idealize what happens in a perfect computer by taking a limit as $\macheps\to 0$.

```{index} ! condition number; of a scalar function
```

(definition-intro-condno)=
::::{prf:definition} Condition number (scalar function)
The relative **condition number** of a scalar function $f(x)$ is 

```{math}
:label: condition
   \kappa_f(x) = \lim_{\epsilon\to 0} \dfrac{ |f(x)-f(x(1+\epsilon))| }{ |\epsilon f(x)| }.
```
::::

The condition number is a ratio of the relative error of the output to the relative error of the input. It depends only on the problem and the data, not the computer or the algorithm.

Assuming that $f$ has at least one continuous derivative, we can simplify the expression {eq}`condition` through some straightforward manipulations:

```{math}
:label: conditionderiv
\begin{split}
   \kappa_f(x) &= \lim_{\epsilon\to 0} \left| \dfrac{ f(x+\epsilon x) - f(x) }{ \epsilon f(x)}
   \right| \\
  &= \lim_{\epsilon\to 0}  \left| \dfrac{ f(x+\epsilon x) - f(x) }{ \epsilon x} \cdot  \frac{x}{f(x)}  \right| \\
   &= \left| \dfrac{ x f'(x)} {f(x)}  \right|.
\end{split}
```

In retrospect, it should come as no surprise that the change in values of $f(x)$ due to small changes in $x$ involves the derivative of $f$. In fact, if we were making measurements of changes in absolute rather than relative terms, the condition number would be simply $|f'(x)|$.

````{prf:example}
Let's return to our "add 1" problem and generalize it slightly to $f(x)=x-c$ for constant $c$. We compute, using {eq}`conditionderiv`,

```{math}
:label: condadd
\kappa_f(x)=\left| \frac{(x)(1)}{x-c} \right| = \left| \frac{x}{x-c}\right|.
```

```{index} subtractive cancellation
```

The result is the relative change {eq}`plus1relerr` normalized by the size of the perturbation $\epsilon$. The condition number is large when $|x|\gg |x-c|$. Considering that $c$ can be negative, this result applies to both addition and subtraction. Furthermore, the situation is symmetric in $x$ and $c$; that is, if we perturbed $c$ and not $x$, the result would be $|c|/|x-c|$.
````

````{prf:example}
Another elementary operation is to multiply by a constant: $f(x)=cx$ for nonzero $c$. We compute
  
```{math}
:label: condmult
\kappa_f(x) = \left| \dfrac{ x f'(x)} {f(x)}  \right| = \left| \frac{(x)(c)}{cx} \right| = 1.
```

We conclude that multiplication by a real number leads to the same  relative error in the result as in the data. In other words,  multiplication does not have the potential for cancellation error that addition does.
````

```{index} condition number; of elementary functions
```

Condition numbers of the major elementary functions are given in {numref}`Table {number} <table-condition-functions>`.

```{list-table} Relative condition numbers of elementary functions
:name: table-condition-functions
:header-rows: 1

* - Function 
  - Condition number
* - $f(x) = x + c$ 
  - $\kappa_f(x) = \dfrac{\lvert x \rvert}{\lvert x+c\rvert}$ 
* - $f(x) = cx$
  - $\kappa_f(x) = 1$ 
* - $f(x) = x^p$ 
  - $\kappa_f(x) = \lvert p \rvert$
* - $f(x) = e^x$
  - $\kappa_f(x) = \lvert x \rvert$
* - $f(x) = \sin(x)$ 
  - $\kappa_f(x) = \lvert x\cot(x) \rvert$ 
* - $f(x) = \cos(x)$
  - $\kappa_f(x) = \lvert x\tan(x) \rvert$
* - $f(x) = \log(x)$
  - $\kappa_f(x) = \dfrac{1}{\lvert \log(x) \rvert}$
```

As you are asked to show in [Exercise 4](problem-conditioning-chain), when two functions $f$ and $g$ are combined in a chain as $h(x)=f(g(x))$, the composite condition number is

```{math}
:label: condition-chain
\kappa_h(x) = \kappa_f(g(x)) \cdot \kappa_g(x).
```

## Estimating errors

Refer back to the definition of $\kappa_f$ as a limit in {eq}`condition`. Approximately speaking, if $|\epsilon|$ is small, we expect 

```{math}
\left| \dfrac{ f(x+\epsilon x) - f(x) }{ f(x)}  \right| \approx \kappa_f(x)\, |\epsilon|.
```

That is, whenever the data $x$ is perturbed by a small amount, we expect that relative perturbation to be magnified by a factor of $\kappa_f(x)$ in the result. 

::::{prf:observation}
If $\kappa_f \approx 10^d$, then we expect to lose up to $d$ decimal digits of accuracy in computing $f(x)$ from $x$.
::::

Large condition numbers signal when errors cannot be expected to remain comparable in size to roundoff error. We call a problem poorly conditioned or **ill-conditioned** when $\kappa_f(x)$ is large, although there is no fixed threshold for the term.

If $\kappa_f \approx 1/\macheps$, then we can expect the result to have a relative error of as much as 100% simply by expressing the data $x$ in finite precision. Such a function is essentially not computable at this machine epsilon.

````{prf:example}
Consider the problem $f(x)= \cos(x)$. By the table above, $\kappa_f(x) = |x \tan x|$. There are two different ways in which $\kappa$ might become large:

- If $|x|$ is very large, then perturbations that are small relative to $x$ may still be large compared to $1$. Because $|f(x)|\le 1$ for all $x$, this implies that the perturbation will be large relative to the result, too.
- The condition number grows without bound as $x$ approaches an odd integer multiple of $\pi/2$, where $f(x)=0$. A perturbation which is small relative to a nonzero $x$ may not be small relative to $f(x)$ in such a case.
````

You may have noticed that for some functions, such as the square root, the condition number can be less than 1. This means that relative changes get *smaller* in the passage from input to output. However, every result in floating-point arithmetic is still subject to rounding error at the relative level of $\macheps$. In practice, $\kappa_f<1$ is no different from $\kappa_f=1$.

## Polynomial roots

Most problems have multiple input and output values. These introduce complications into the formal definition of the condition number. Rather than worry over those details here, we can still look at variations in only one output value with respect to one data value at a time.

(example-quad-root-cond)=
````{prf:example}
Consider the problem of finding the roots of a quadratic polynomial; that is, the values of $t$ for which $at^2+bt+c=0$. Here the data are the coefficients $a$, $b$, and $c$ that define the polynomial, and the solution to the problem are the two (maybe complex-valued) roots $r_1$ and $r_2$. Formally, we might write $f([a,b,c])=[r_1,r_2]$ using vector notation.

Let's pick one root $r_1$ and consider what happens to it as we vary just the leading coefficient $a$. This suggests a scalar function $f(a)=r_1$. Starting from $ar_1^2 + br_1 + c = 0$, we differentiate implicitly with respect to $a$ while holding $b$ and $c$ fixed:  

```{math}
r_1^2 + 2a r_1 \left(\frac{dr_1}{da}\right) + b \,\frac{dr_1}{da} = 0.
```

Solving for the derivative, we obtain
  
```{math}
:label: rootderiv1
\frac{dr_1}{da} = \frac{-r_1^2}{2a r_1 + b}. 
```

Hence the condition number for the problem $f(a)=r_1$ is

```{math}
:label: rootcond1
\kappa_f(a)  = \left|\frac{a}{r_1} \cdot \frac{dr_1}{da} \right| = \left| \frac{a r_1}{ 2a r_1 + b} \right|
= \left| \frac{r_1}{ r_1-r_2} \right|,
```

where in the last step we used the quadratic formula: 

:::{math}
|2ar_1 + b | = \left| \sqrt{b^2-4ac} \, \right| = |a(r_1 - r_2)|.
:::

Based on {eq}`rootcond1`, we can expect poor conditioning in the rootfinding problem if and only if $|r_1| \gg |r_1-r_2|$. Similar conclusions apply for $r_2$ and for variations with respect to the coefficients $b$ and $c$.
````

The calculation in {numref}`Example %s <example-quad-root-cond>` generalizes to polynomials of higher degree.

:::{prf:observation}
Roots of polynomials are ill-conditioned with respect to changes in the polynomial coefficients when they are much closer to each other than to the origin.
:::

The condition number of a root can be arbitrarily large. In the extreme case of a repeated root, the condition number is formally infinite, which implies that the ratio of changes in the root to changes in the coefficients cannot be bounded. 

(demo-condition-roots)=
```````{prf:example}
``````{tab-set}
`````{tab-item} Julia
:sync: julia
````{dropdown} Conditioning of polynomial roots
:open:
```{include} julia/roots.ipynb
```
````
`````

`````{tab-item} MATLAB
:sync: matlab
````{dropdown} Conditioning of polynomial roots
:open:
```{include} matlab/roots.ipynb
```
````
`````

`````{tab-item} Python
:sync: python
````{dropdown} Conditioning of polynomial roots
:open:
```{include} python/roots.ipynb
```
````
`````
``````
```````
    

## Exercises

1. ✍ Use {eq}`conditionderiv` to derive the relative condition numbers of the following functions appearing in {numref}`table-condition-functions`.

    **(a)** $f(x) = x^p,\quad$
    **(b)** $f(x) = \log(x),\quad$
    **(c)** $f(x) = \cos(x),\quad$
    **(d)** $f(x) = e^x$.

2. ✍ Use the chain rule {eq}`condition-chain` to find the relative condition number of the given function. Then check your result by applying {eq}`conditionderiv` directly.

    **(a)** $f(x) = \sqrt{x+5},\quad$
    **(b)** $f(x) = \cos(2\pi x),\quad$
    **(c)** $f(x) = e^{-x^2}.$

3. ✍ Calculate the relative condition number of each function, and identify all values of $x$ at which $\kappa_{f}(x)\to\infty$ (including limits as $x\to\pm\infty$).

    **(a)** $f(x) = \tanh(x),\quad$
    **(b)** $f(x) = \dfrac{e^x-1}{x},\quad$
    **(c)** $f(x) = \dfrac{1-\cos(x)}{x}.$

    (problem-conditioning-chain)=
4. ✍ Suppose that $f$ and $g$ are real-valued functions that have relative condition numbers $\kappa_f$ and $\kappa_g$, respectively. Define a new function $h(x)=f\bigl(g(x)\bigr)$. Show that for $x$ in the domain of $h$, the relative condition number of $h$ satisfies {eq}`condition-chain`.

5. ✍ Suppose that $f$ is a function with relative condition number $\kappa_f$, and that $f^{-1}$ is its inverse function. Show that the relative condition number of $f^{-1}$ satisfies
  
    ```{math}
    \kappa_{f^{-1}}(x) = \frac{1}{\kappa_f\Bigl( f^{-1}(x) \Bigr)},
    ```

    provided the denominator is nonzero.

   (problem-quadrootcond)=
6. ✍  Referring to the derivation of {eq}`rootcond1`, derive an expression for the relative condition number of a root of $ax^2+bx+c=0$ due to perturbations in $b$ only. 

7. The polynomial $x^2-2x+1$ has a double root at 1. Let $r_1(\epsilon)$ and $r_2(\epsilon)$ be the roots of the perturbed polynomial $x^2-(2+\epsilon)x+1$.
  
    **(a)** ✍/⌨ Using a computer or calculator, make a table with rows for $\epsilon = 10^{-2}$, $10^{-4}$, $10^{-6}$, $\ldots$, $10^{-12}$ and columns for $\epsilon$, $r_1(\epsilon)$, $r_2(\epsilon)$, $|r_1(\epsilon)-1|$, and $|r_2(\epsilon)-1|$.

    **(b)** ✍ Show that the observations of part (a) satisfy
    
    $$\max\{\, |r_1(\epsilon)-1|, |r_2(\epsilon)-1| \,\} \approx C \epsilon^q$$
    
    for some $0<q<1$. (This supports the conclusion that $\kappa=\infty$ at the double root.) 
  
8. ✍ Generalize {eq}`rootcond1` to finding a root of the $n$th degree polynomial $p(x) = a_nx^n + \cdots + a_1 x + a_0$, and show that the relative condition number of a root $r$ with respect to perturbations only in $a_k$ is
  
    ```{math}
    \kappa_r(a_k) = \left| \frac{a_k r^{k-1}}{p'(r)} \right|.
    ```
