# Trigonometric interpolation
\label{sec:triginterp}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Up to this point all of our global approximating functions have been polynomials. While they are versatile and easy to work with, they are not always the best choice.

Suppose we want to approximate a function $f$ that is periodic, with one period represented by the standard interval $[-1,1]$. Mathematically, $f(x+2)=f(x)$ for all real $x$. We could use polynomials to interpolate or project $f$. However, it seems reasonable to replace polynomials by functions that are also periodic, i.e., trigonometric functions.

Doing so leads to 
```{index} interpolation!by trigonometric functions
```
 {term}`trigonometric interpolation`. As it happens, \texthighlight{trigequi}{trigonometric interpolation allows us to return to equally spaced nodes without any problems.} We therefore define $N=2n+1$ equally spaced nodes inside the interval $[-1,1]$ by

$$
  :label: trignodes
  t_k = \frac{2k}{N}, \quad k=-n,\ldots,n.
$$

The formulas in this section require some minor but important adjustments if $N$ is even instead. We have modified our standard indexing scheme here to make the symmetry within $[-1,1]$ about $x=0$ more transparent. Note that the endpoints $\pm1$ are *not* among the nodes.

As usual, we have sample values $y_{-n},\ldots,y_n$, perhaps representing values of a function $f(x)$ at the nodes.  We also now assume that the sample values can be extended periodically forever in both directions, so that $y_{k+mN}=y_k$ for any integer $m$.

\fxnote{Figure for the periodic interpolation setup?}
%
::::{prf:example}
%  \inputexample{globalfuncapprox}{periodicfunction}
%::::

%
%See Figure~\ref{fig:trigcardinal}.
%\begin{figure}[tbh]
%  \centering
%  \psfrag{x}[t][t]{$x$}
%  \includegraphics[width=\textwidth]{trigcardinal}
%  \caption{Periodic interpolation. The top shows sample values and
%    their underlying function within the interval $[-1,1]$ for $n=5$
%    ($N=11$ nodes), and how both can be extended periodically to the
%    entire real line. The bottom shows a cardinal trigonometric
%    interpolant.}
%  \label{fig:trigcardinal}
%\end{figure}


## Cardinal functions

We can explicitly state the 
```{index} cardinal functions
```
 cardinal function basis for equispaced trigonometric interpolation. It starts with

$$
  :label: trigcardinal
  \tau(x) = \frac{2}{N} \left( \frac{1}{2} + \cos \pi x + \cos 2\pi x
    + \cdots \cos n\pi x\right) = \frac{\sin(N\pi x/2)}{N\sin(\pi x/2)}.
$$

%(Proof of this summation is a cute exercise in complex analysis.)
You can directly check (see {ref}`prob-triginterp-checktau`) that this is 2-periodic, that $\tau(t_k)=0$ for $k\neq 0$, and that $\tau(t_0)=\tau(0)=1$ in the limiting sense.

Because any shift of a periodic function is also periodic, the cardinal basis for trigonometric interpolation is defined by $\tau_k(x) = \tau(x-t_k)$. Because the function $\tau_{-n},\ldots,\tau_n$ form a cardinal basis, the coefficients of the interpolant are just the sampled function values:

$$
  :label: trigcardinalinterp
  p(x)=\sum_{k=-n}^n y_k\tau_k(x).
$$


{numref}`Function {number}<function-triginterp>` is an implementation of trigonometric interpolation based on {eq}`trigcardinalinterp`. The function accepts an $N$-vector of equally spaced nodes. (The case of even $N$ is included in the code.) Note too that evaluation of $\tau(0)=1$ from {eq}`trigcardinal` properly requires an application of L'H\^opital's rule, but we patch it after the fact by replacing "NaN" values.

(function-triginterp)=
````{proof:function} triginterp
**Trigonometric interpolation**

```{code-block} julia
:lineno-start: 1
```
````

::::{prf:example}
  \inputexample{globalfuncapprox}{periodicinterp}
::::



%## Gibbs phenomenon

%As always, spectral convergence is predicated on having infinitely many continuous derivatives. At the other extreme is a function with a jump discontinuity. Trigonometric interpolation across a jump leads to a lack of convergence altogether, a fact famously known as the **Gibbs phenomenon**.
%
::::{prf:example}
%  \inputexample{globalfuncapprox}{gibbs}
%::::


%{ref}`example-gibbs` demonstrates that at a discontinuity, the trigonometric interpolant overshoots the correct values and then has a slow oscillatory decay. As $N$ grows, the widths of the overshoots go to zero, but the heights of the overshoots converge not to zero but to about $0.28$, or $14\%$ of the step height. A similar phenomenon affects the interpolation of discontinuous nonperiodic functions by polynomials.

## Fast Fourier transform

Although the cardinal form of the interpolant is useful and stable, there is also a lot of importance attached to the equivalent forms

$$
  :label: triginterp
  p(x) = \sum_{k=-n}^n c_k e^{ik\pi x} = \frac{a_0}{2} + \sum_{k=1}^n \Bigl[ a_k \cos(k\pi x) + b_k \sin(k\pi x) \Bigr],
$$

where the constants $c_k$, or alternatively the constants $a_k$ and $b_k$, are determined by interpolation conditions. The connection between real and complex versions are Euler's formula,

$$
  :label: eulerformula
  e^{i\theta} = \cos(\theta) + i\sin(\theta),
$$

and the consequential

$$
  :label: eulersincos
  \cos \theta = \frac{e^{i \theta}+e^{-i\theta}}{2}, \qquad \sin \theta = \frac{e^{i \theta}-e^{-i\theta}}{2i}.
$$

While working with an all-real formulation seems natural when the data are real, the complex-valued version leads to more elegant formulas and standard usages, so we adopt it exclusively.

%and its generalization,
%
$$
%  :label: eulerk
%  \bigl(e^{i\theta}\bigr)^k = e^{ik\theta} = \cos(k\theta) + i \sin(k\theta).
%$$

%Note that the complex form can be used even when the function $f$ is real, though that implies \fxnote*{Make an exercise of this!}{some constraints} on the complex-valued coefficients $c_k$.
%Thinking in terms of $z^k$ for $z=\exp(i\theta)$,
%equation {eq}`eulerk` explains why {eq}`triginterp` is often referred to as a **trigonometric polynomial**. Its periodicity is an automatic consequence of its construction.

The $N=2n+1$ unknown coefficients are determined by interpolation nodes at the $N$ nodes within $[-1,1]$. By evaluating the complex exponential functions at these nodes, we get the $N\times N$ linear system
\[
   \begin{bmatrix}
     e^{-i n\pi t_{-n}}   & \cdots & 1      & e^{i \pi t_{-n}}   & \cdots & e^{i n\pi t_{-n}}   \\
     e^{-i n\pi t_{-n+1}} & \cdots & 1      & e^{i \pi t_{-n+1}} & \cdots & e^{i n\pi t_{-n+1}} \\
     \vdots             &        & \vdots & \vdots            &        & \vdots              \\[1mm]
     e^{-i n\pi t_{0}}    & \cdots & 1      & e^{i \pi t_{0}}    & \cdots & e^{i n\pi t_{0}}    \\
     \vdots             &        & \vdots & \vdots             &        & \vdots              \\[1mm]
     e^{-i n\pi t_{n}}    & \cdots & 1      & e^{i \pi t_{n}}    & \cdots & e^{i n\pi t_{n}}    \\
   \end{bmatrix}
   \mathbf{c} = \mathbf{y},
\]
to be solved for the coefficients $\mathbf{c}$. One of the most important (though not entirely original) observations of the 20th century was that \texthighlight{fft}{the linear system for the interpolation coefficients can be solved in just $O(N\log N)$ operations by an algorithm now known as the 
```{index} FFT
```
 **fast Fourier transform** or FFT.} \matlab has a function {term}`fft` to perform this transform, but its conventions are a little different from ours. Instead of nodes in $(-1,1)$, it expects the nodes to be defined in $[0,2)$, and it returns the complex coefficients
\[
N
\begin{bmatrix}
  c_0 & c_1 & \cdots & c_n & c_{-n} & \cdots & c_{-1}
\end{bmatrix}^T.
\]


::::{prf:example}
  \inputexample{globalfuncapprox}{fft_ex}
::::


The theoretical and computational aspects of Fourier analysis are vast and far-reaching. We have given only the briefest of introductions.

\begin{exercises}
  \input{globalfuncapprox/exercises/Trig}
\end{exercises}
\clearpage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
