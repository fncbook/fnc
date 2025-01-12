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
[**Demo %s**](#demo-flops-mvmult)

Here is a straightforward implementation of matrix-vector multiplication.

```{code-cell}
n = 6;
A = magic(n);
x = ones(n,1);
y = zeros(n,1);
for i = 1:n
    for j = 1:n
        y(i) = y(i) + A(i,j)*x(j);   % 2 flops
    end
end
```

Each of the loops implies a summation of flops. The total flop count for this algorithm is

$$ \sum_{i=1}^n \sum_{j=1}^n 2 = \sum_{i=1}^n 2n = 2n^2. $$

Since the matrix $\mathbf{A}$ has $n^2$ elements, all of which have to be involved in the product, it seems unlikely that we could get a flop count that is smaller than $O(n^2)$ in general.

```{index} ! MATLAB; tic and toc
```

Let's run an experiment with the built-in matrix-vector multiplication, using `tic` and `toc` to time the operation.

```{code-cell}
n_ = (400:400:4000)';
t_ = zeros(size(n_));
for i = 1:length(n_)
    n = n_(i);
    A = randn(n, n);  x = randn(n, 1);
    tic    % start a timer
    for j = 1:100      % repeat 100 times
        A*x;
    end
    t = toc;           % read the timer
    t_(i) = t / 100;   % seconds per instance
end
```

The reason for doing multiple repetitions at each value of $n$ in the loop above is to avoid having times so short that the resolution of the timer is significant.

```{code-cell}
table(n_, t_, 'variablenames', {'size', 'time'})
```

```{index} MATLAB; Boolean indexing
```

Looking at the timings just for $n=2000$ and $n=4000$, they have ratio
```{tip}
:class: dropdown
The expression `n_==4000` here produces a vector of Boolean (true/false) values the same size as `n_`. This result is used to index within `t_`, accessing only the value for which the comparison is true.
```

```{code-cell}
t_(n_==4000) / t_(n_==2000)
```

If the run time is dominated by flops, then we expect this ratio to be 

$$
\frac{2(4000)^2}{2(2000)^2}=4.
$$
