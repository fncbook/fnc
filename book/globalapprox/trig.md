# Trigonometric interpolation

Up to this point all of our global approximating functions have been polynomials. While they are versatile and easy to work with, they are not always the best choice.

Suppose we want to approximate a function $f$ that is periodic, with one period represented by the standard interval $[-1,1]$. Mathematically, $f(x+2)=f(x)$ for all real $x$. We could use polynomials to interpolate or project $f$. However, it seems reasonable to replace polynomials by functions that are also periodic, i.e., trigonometric functions.

```{index} interpolation; by trigonometric functions
```
Doing so leads to {term}`trigonometric interpolation`. As it happens, trigonometric interpolation allows us to return to equally spaced nodes without any problems. We therefore define $N=2n+1$ equally spaced nodes inside the interval $[-1,1]$ by

```{margin}
Trigonometric interpolation works on equispaced nodes without stability problems.
```

:::{math}
  :label: trignodes
  t_k = \frac{2k}{N}, \quad k=-n,\ldots,n.
:::

The formulas in this section require some minor but important adjustments if $N$ is even instead. We have modified our standard indexing scheme here to make the symmetry within $[-1,1]$ about $x=0$ more transparent. Note that the endpoints $\pm1$ are *not* among the nodes.

As usual, we have sample values $y_{-n},\ldots,y_n$, perhaps representing values of a function $f(x)$ at the nodes.  We also now assume that the sample values can be extended periodically forever in both directions, so that $y_{k+mN}=y_k$ for any integer $m$.

## Cardinal functions

```{index} cardinal functions
```
We can explicitly state the cardinal function basis for equispaced trigonometric interpolation. It starts with

:::{math}
  :label: trigcardinal
  \tau(x) = \frac{2}{N} \left( \frac{1}{2} + \cos \pi x + \cos 2\pi x
    + \cdots \cos n\pi x\right) = \frac{\sin(N\pi x/2)}{N\sin(\pi x/2)}.
:::

You can directly check (see [this exercise](problem-triginterp-checktau)) that this is 2-periodic, that $\tau(t_k)=0$ for $k\neq 0$, and that $\tau(t_0)=\tau(0)=1$ in the limiting sense.

Because any shift of a periodic function is also periodic, the cardinal basis for trigonometric interpolation is defined by $\tau_k(x) = \tau(x-t_k)$. Because the function $\tau_{-n},\ldots,\tau_n$ form a cardinal basis, the coefficients of the interpolant are just the sampled function values:

:::{math}
  :label: trigcardinalinterp
  p(x)=\sum_{k=-n}^n y_k\tau_k(x).
:::

{numref}`Function {number}<function-triginterp>` is an implementation of trigonometric interpolation based on {eq}`trigcardinalinterp`. The function accepts an $N$-vector of equally spaced nodes. (The case of even $N$ is included in the code.) Note too that evaluation of $\tau(0)=1$ from {eq}`trigcardinal` properly requires an application of L'H\^opital's rule, but we patch it after the fact by replacing "NaN" values.

(function-triginterp)=
````{proof:function} triginterp
**Trigonometric interpolation**
```{code-block} julia
:lineno-start: 1
```
````

::::{prf:example} Julia demo
:class: demo
:label: demos-trig-interp
{doc}`demos/trig-interp`
::::

## Fast Fourier transform

Although the cardinal form of the interpolant is useful and stable, there is also a lot of importance attached to the equivalent forms

:::{math}
  :label: triginterp
  p(x) = \sum_{k=-n}^n c_k e^{ik\pi x} = \frac{a_0}{2} + \sum_{k=1}^n  a_k \cos(k\pi x) + b_k \sin(k\pi x) ,
:::

where the complex constants $c_k$, or alternatively the real constants $a_k$ and $b_k$, are determined by interpolation conditions. The connection between real and complex versions are Euler's formula,

:::{math}
  :label: eulerformula
  e^{i\theta} = \cos(\theta) + i\sin(\theta),
:::

and the resultant identities

:::{math}
  :label: eulersincos
  \cos \theta = \frac{e^{i \theta}+e^{-i\theta}}{2}, \qquad \sin \theta = \frac{e^{i \theta}-e^{-i\theta}}{2i}.
:::

While working with an all-real formulation seems natural when the data are real, the complex-valued version leads to more elegant formulas and standard usages, so we adopt it exclusively.

The $N=2n+1$ unknown coefficients are determined by interpolation nodes at the $N$ nodes within $[-1,1]$. By evaluating the complex exponential functions at these nodes, we get the $N\times N$ linear system

$$
   \begin{bmatrix}
     e^{-i n\pi t_{-n}}   & \cdots & 1      & e^{i \pi t_{-n}}   & \cdots & e^{i n\pi t_{-n}}   \\
     e^{-i n\pi t_{-n+1}} & \cdots & 1      & e^{i \pi t_{-n+1}} & \cdots & e^{i n\pi t_{-n+1}} \\
     \vdots             &        & \vdots & \vdots            &        & \vdots              \\[1mm]
     e^{-i n\pi t_{0}}    & \cdots & 1      & e^{i \pi t_{0}}    & \cdots & e^{i n\pi t_{0}}    \\
     \vdots             &        & \vdots & \vdots             &        & \vdots              \\[1mm]
     e^{-i n\pi t_{n}}    & \cdots & 1      & e^{i \pi t_{n}}    & \cdots & e^{i n\pi t_{n}}    \\
   \end{bmatrix}
   \mathbf{c} = \mathbf{y},
$$

```{index} FFT
```
to be solved for the coefficients $\mathbf{c}$. One of the most important (though not entirely original) algorithmic observations of the 20th century was that the linear system for the interpolation coefficients can be solved in just $O(N\log N)$ operations by an algorithm now known as the **fast Fourier transform** or FFT.

The `FFTW` package provides a function `fft` to perform this transform, but its conventions are a little different from ours. Instead of nodes in $(-1,1)$, it expects the nodes to be defined in $[0,2)$, and it returns the complex coefficients

$$
\begin{bmatrix}
  c_0, & c_1, & \cdots & c_n, & c_{-n}, & \cdots & c_{-1}
\end{bmatrix}.
$$

::::{prf:example} Julia demo
:class: demo
:label: demos-trig-interp
{doc}`demos/trig-fft`
::::

The theoretical and computational aspects of Fourier analysis are vast and far-reaching. We have given only the briefest of introductions.

\begin{exercises}
  \input{globalfuncapprox/exercises/Trig}
\end{exercises}
