---
numbering:
  enumerator: 9.7.%s
---
(section-globalapprox-improper)=
# Improper integrals

```{index} improper integral
```

When the interval of integration or the integrand itself is unbounded, we say an integral is *improper*. Improper integrals present particular challenges to numerical computation.

## Infinite interval

When the integration domain is $(-\infty,\infty)$, the integrand has to decay as $x \to \pm \infty$ in order for the improper integral to be finite. This fact brings up the possibility of truncating the domain:

:::{math}
:label: infiniteinterval
\int_{-\infty}^{\infty} f(x)\, dx \approx  \int_{-M}^{M} f(x)\, dx.
:::

This integral can be discretized finitely by, say, the trapezoid formula, or an adaptive integrator. However, this approach can be inefficient.

(example-improper-slowdecay)=
::::{prf:example}
Consider the integral of $f(x)=1/(1+x^2)$,

$$
  \int_{-\infty}^\infty \frac{1}{1+x^{2}}\, dx = \pi.
$$

In this case we can easily estimate the effect of truncation on the result. For large $M$, 

$$
\int_M^\infty f\,dx \approx \int_M^\infty x^{-2}\,dx =  M^{-1}.
$$

The same estimate applies to the integral over $(-\infty,-M)$. To get eight digits of accuracy, for instance, we need to truncate with $M > 2 \times 10^8$. 
::::

## Double exponential transformation

In order to do better than direct truncation, we want to encourage the function to decay faster. In practice this means a change of variable, $x(t)$. If $|x(t)|$ grows rapidly as $|t| \to \infty$, then $f(x(t))$ will decay more rapidly in $t$ than in $x$. 

One way to accomplish this feat is to use

:::{math}
:label: DEquadtrans1
  x(t) = \sinh\left(  \sinh t \right).
:::

Noting the asymptotic behavior as $t \rightarrow \pm\infty$ that

:::{math}
:label: sinh-asymp
  \left| \sinh(t) \right| \sim \frac{1}{2} e^{ |t| },
:::

we find that in the same limits,

$$
x(t) \approx \pm \frac{1}{2} \exp\left( \frac{1}{2} e^{ |t| } \right).
$$

```{index} ! double exponential transformation
```
Thus, {eq}`DEquadtrans1` is often referred to as a **double exponential** transformation.

By the chain rule,

:::{math}
:label: DEquadchain1
\begin{split}
\int_{-\infty}^\infty f(x)\, dx &= \int_{-\infty}^\infty f(x(t))\frac{dx}{dt}\, dt \\
  &= \int_{-\infty}^\infty f(x(t))\, \cosh\left( \sinh t \right)  \cosh t  \, dt.
\end{split}
:::

The exponential terms introduced by the chain rule grow double exponentially, but the more rapid decay of $f$ in the new variable more than makes up for this.

(demo-improper-decay)=
::::{prf:example} Decay by transformation
Consider again $f(x)=1/(1+x^2)$ from {numref}`Example %s <example-improper-slowdecay>`, with $x(t)$ given by {eq}`DEquadtrans1`. As $t\to\infty$, 

$$
f(x(t)) \approx x^{-2} \approx 4 \exp\left( -e^{t} \right).
$$

The chain rule terms in {eq}`DEquadchain1` become

$$
\cosh\left( \sinh t \right)  \cosh t 
\approx \frac{1}{2} \exp\left( \frac{1}{2} e^t  \right)  \cdot \frac{1}{2} e^t,
$$ 

yielding a product that is roughly

$$
2 \exp\left( -\frac{1}{2} e^t  \right).
$$

The total integrand in {eq}`DEquadchain1` therefore has double exponential decay in $t$, essentially because of the squaring of $x$ in the denominator of $f$. The same result holds as $t\to-\infty$.

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-improper-decay-julia
:::
````

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-improper-decay-matlab
:::
````

````{tab-item} Python
:sync: python
:::{embed} #demo-improper-decay-python
:::
````
`````
::::

{numref}`Function {number} <function-intinf>` implements double exponential integration by applying the adaptive integrator {numref}`Function {number} <function-intadapt>` to {eq}`DEquadchain1`. It truncates the interval to $-M\le t \le M$ by increasing $M$ until the integrand is too small to matter relative to the error tolerance.

(function-intinf)=
``````{prf:algorithm} intinf
`````{tab-set} 
````{tab-item} Julia
:sync: julia
:::{embed} #function-intinf-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #function-intinf-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #function-intinf-python
:::
````
`````
``````

(demo-improper-intinf)=
::::{prf:example} Infinite interval
We compare direct truncation in $x$ to the double exponential method of {numref}`Function {number} <function-intinf>` for $f(x)=1/(1+x^2)$.


`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-improper-intinf-julia
:::
````

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-improper-intinf-matlab
:::
````

````{tab-item} Python
:sync: python
:::{embed} #demo-improper-intinf-python
:::
````
`````
::::

## Integrand singularity

If $f$ asymptotically approaches infinity as $x$ approaches an integration endpoint, its exact integral may or may not be finite. If $f$ is integrable, then the part of the integration interval near the singularity needs to be more finely resolved than the rest of it. 

Let's consider 

:::{math}
:label: intsing
\int_0^1 f(x)\,dx,
:::

where $f$ and/or a derivative of $f$ is unbounded at the left endpoint, zero. The change of variable

:::{math}
:label: DEquadtrans2
  x(t) = \frac{2}{1+\exp(2 \sinh t)}
:::

satisfies $x(0)=1$ and $x\to 0^+$ as $t\to \infty$, thereby transforming the integration interval to $t\in(0,\infty)$ and placing the singularity at infinity. The chain rule implies

:::{math}
:label: DEquadchain2
\begin{split}
  \int_{0}^1 f(x)\, dx &= \int_{0}^\infty f(x(t)) \frac{dx}{dt}\, dt \\
  &= \int_{0}^\infty f(x(t)) \frac{\cosh t}{\cosh(\sinh t)^2}  \,  dt.
\end{split}
:::

Now the growth of $f$ and $\cosh t$ together are counteracted by the double exponential denominator, allowing easy truncation of {eq}`DEquadchain2`. This variable transformation is paired with adaptive integration in {numref}`Function {number} <function-intsing>`.

(function-intsing)=
``````{prf:algorithm} intsing
`````{tab-set} 
````{tab-item} Julia
:sync: julia
:::{embed} #function-intsing-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #function-intsing-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #function-intsing-python
:::
````
`````
``````

(demo-improper-intsing)=
::::{prf:example} Singularity at an endpoint
Let's use {numref}`Function {number} <function-intsing>` to compute

$$
\int_0^{0.01} \frac{1}{\sqrt{x}}\, dx.
$$

Since the integration interval is not $[0,1]$, we must first use the change of variable $s=100t$, yielding

$$
\int_0^{1} \frac{1}{10\sqrt{s}}\, ds = 0.2.
$$

In order to use {numref}`Function {number} <function-intadapt>`, we must truncate on the left to avoid evaluation at zero, where $f$ is infinite. Since the integral from $0$ to $\delta$ is $20\sqrt{\delta}$, we use $\delta=(\epsilon/20)^2$ to achieve error tolerance $\epsilon$.
`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-improper-intsing-julia
:::
````

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-improper-intsing-matlab
:::
````

````{tab-item} Python
:sync: python
:::{embed} #demo-improper-intsing-python
:::
````
`````
::::

Double exponential integration is an effective general-purpose technique for improper integrals that usually outperforms interval truncation in the original variable. There are specialized methods tailored to specific singularity types that can best it, but those require more analytical work to use properly.

## Exercises


``````{exercise}
⌨ Use {numref}`Function {number} <function-intinf>` to estimate the given integral with error tolerances $10^{-3},10^{-6},10^{-9},10^{-12}$. For each result, show the actual error and the number of nodes used.

**(a)** $\displaystyle\int_{-\infty}^\infty \dfrac{1}{1+x^2+x^4}\, dx = \dfrac{\pi}{\sqrt{3}}$

**(b)** $\displaystyle\int_{-\infty}^\infty e^{-x^2}\cos(x)\, dx = e^{-1/4}\sqrt{\pi}$

**(c)** $\displaystyle\int_{-\infty}^\infty (1+x^2)^{-2/3}\, dx = \dfrac{\sqrt{\pi}\,\Gamma(1/6)}{\Gamma(2/3)}$  (use `gamma()` for $\Gamma()$)

``````


``````{exercise}
⌨ Use {numref}`Function {number} <function-intsing>` to estimate the given integral, possibly after rewriting the integral into the form {eq}`intsing` with a left-endpoint singularity. Use error tolerances $10^{-3},10^{-6},10^{-9},10^{-12}$, and for each result, show the actual error and the number of nodes used.

**(a)** $\displaystyle\int_{0}^1 (\log x)^2\, dx = 2$

**(b)** $\displaystyle\int_{0}^{\pi/4} \sqrt{\tan(x)}\, dx = \dfrac{\pi}{\sqrt{2}}$

**(c)** $\displaystyle\int_{0}^1 \frac{1}{\sqrt{1-x^2}}\, dx = \dfrac{\pi}{2}$
``````


``````{exercise}
For integration on a semi-infinite interval such as $x\in [0,\infty)$, another double exponential transformation is useful: $x(t)=\exp\left( \sinh t \right)$.

**(a)** ✍ Show that $t\in(-\infty,\infty)$ is mapped to $x\in (0,\infty)$. 

**(b)** ✍ Derive an analog of {eq}`DEquadchain1` for the chain rule on $\int_0^\infty f(x)\,dx$. 

**(c)** ✍ Show that truncation of $t$ to $[-M,M]$ will truncate $x$ to $[1/\mu,\mu]$ for some positive $\mu$. 

**(d)** ⌨ Write a function `intsemi(f,tol)` for the semi-infinite integration problem. Test it on the integral

$$  
\displaystyle\int_0^\infty \frac{e^{-x}}{\sqrt{x}}\,dx = \sqrt{\pi}.
$$


``````
