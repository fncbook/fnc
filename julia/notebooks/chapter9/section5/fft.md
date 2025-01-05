---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
This function has frequency content at $2\pi$, $-2\pi$, and $\pi$. 

```{code-cell}
f(x) = 3 * cospi(2x) - cispi(x)    # cispi(x) := exp(1im * π * x)
```

To use `fft`, we set up nodes in the interval $[0,2)$. 

```{code-cell}
n = 4;  N = 2n+1;
t = [ 2j / N for j in 0:N-1 ]      # nodes in [0,2)
y = f.(t);
```

We perform Fourier analysis using `fft` and then examine the resulting coefficients.

```{code-cell}
using FFTW
c = fft(y) / N
freq = [0:n; -n:-1]
@pt :header=["k", "coefficient"] [freq round.(c, sigdigits=5)]
```

Note that $1.5 e^{2i\pi x}+1.5 e^{-2i\pi x} = 3 \cos(2\pi x)$, so this result is sensible.

Fourier's greatest contribution to mathematics was to point out that *every* periodic function is just a combination of frequencies—infinitely many of them in general, but truncated for computational use. Here we look at the magnitudes of the coefficients for $f(x) = \exp( \sin(\pi x) )$.

```{code-cell}
:tags: [hide-input]
f(x) = exp( sin(pi*x) )     # content at all frequencies
n = 9;  N = 2n+1;
t = [ 2j / N for j in 0:N-1 ]      # nodes in [0,2)
c = fft(f.(t)) / N

freq = [0:n; -n:-1]
scatter(freq, abs.(c);
    xaxis=(L"k", [-n, n]),  yaxis=(L"|c_k|", :log10), 
    title="Fourier coefficients",  legend=:none)
```

The Fourier coefficients of smooth functions decay exponentially in magnitude as a function of the frequency. This decay rate is determines the convergence of the interpolation error.
