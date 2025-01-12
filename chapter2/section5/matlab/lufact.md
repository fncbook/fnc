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
[**Demo %s**](#demo-flops-lufact)

We'll test the conclusion of $O(n^3)$ flops experimentally, using the built-in `lu` function instead of the purely instructive `lufact`.
```{tip}
:class: dropdown
The first time a function is invoked, there may be significant time needed to compile it in memory. Thus, when timing a function, run it at least once before beginning the timing.
```

```{code-cell}
n_ = (200:100:2400)';
t_ = zeros(size(n_));
for i = 1:length(n_)
    n = n_(i);
    A = randn(n, n);  
    tic    % start a timer
    for j = 1:6,  [L, U] = lu(A);  end
    t = toc;
    t_(i) = t / 6;  
end
```

We plot the timings on a log-log graph and compare it to $O(n^3)$. The result could vary significantly from machine to machine, but in theory the data should start to parallel the line as $n\to\infty$.

```{code-cell}
clf
loglog(n_,t_,'.-')
hold on, loglog(n_,t_(end)*(n_/n_(end)).^3,'--')
axis tight
xlabel('size of matrix'), ylabel('time (sec)')
title('Timing of LU factorization')
legend('lu','O(n^3)','location','southeast');
```
