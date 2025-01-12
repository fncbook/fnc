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
[**Demo %s**](#demo-trig-fft)

This function has frequency content at $2\pi$, $-2\pi$, and $\pi$. 

```{code-cell}
f = @(x) 3 * cos(2*pi * x) - exp(1i*pi * x);
```

To use `fft`, we set up nodes in the interval $[0,2)$. 

```{code-cell}
n = 4;
N = 2*n + 1;
t = 2 * (0:N-1)' / N;      % nodes in $[0,2)$
y = f(t);
```

We perform Fourier analysis using `fft` and then examine the resulting coefficients.

```{code-cell}
c = fft(y) / N;
freq = [0:n, -n:-1]';
format short
disp(table(freq, c, variableNames=["k", "coefficient"]))
```

Note that $1.5 e^{2i\pi x}+1.5 e^{-2i\pi x} = 3 \cos(2\pi x)$, so this result is sensible.

Fourier's greatest contribution to mathematics was to point out that *every* periodic function is just a combination of frequencies—infinitely many of them in general, but truncated for computational use. Here we look at the magnitudes of the coefficients for $f(x) = \exp( \sin(\pi x) )$.

```{code-cell}
:tags: [hide-input]
f = @(x) exp( sin(pi*x) );    % content at all frequencies
n = 9;  N = 2*n + 1;
t = 2 * (0:N-1)' / N;         % nodes in $[0,2)$
y = f(t);
c = fft(y) / N;
freq = [0:n, -n:-1]';

clf
semilogy(freq, abs(c), 'o')
xlabel('k'),  ylabel('|c_k|')   
title('Fourier coefficients')    
```

The Fourier coefficients of smooth functions decay exponentially in magnitude as a function of the frequency. This decay rate is determines the convergence of the interpolation error.
