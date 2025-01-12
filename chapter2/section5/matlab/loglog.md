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
[**Demo %s**](#demo-flops-loglog)


Let's repeat the previous experiment for more, and larger, values of $n$.

```{code-cell}
n_ = (400:400:6000)';
t_ = zeros(size(n_));
for i = 1:length(n_)
    n = n_(i);
    A = randn(n, n);  x = randn(n, 1);
    tic    % start a timer
    for j = 1:100      % repeat ten times
        A*x;
    end
    t = toc;          % read the timer
    t_(i) = t / 100;   % seconds per instance
end
```

Plotting the time as a function of $n$ on log-log scales is equivalent to plotting the logs of the variables.

```{code-cell}
clf    % clear any existing figure
loglog(n_, t_, '.-')
xlabel('size of matrix')
ylabel('time (sec)')
title(('Timing of matrix-vector multiplications'));
```

You can see that while the full story is complicated, the graph is trending to a straight line of positive slope. For comparison, we can plot a line that represents $O(n^2)$ growth exactly. (All such lines have slope equal to 2.)

```{code-cell}
hold on
loglog(n_, t_(1) * (n_ / n_(1)).^2, '--')
axis tight
legend('data', 'O(n^2)', 'location', 'southeast');
```
