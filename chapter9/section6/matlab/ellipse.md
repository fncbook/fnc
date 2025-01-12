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
[**Demo %s**](#demo-integration-ellipse)

```{code-cell}
f = @(t) pi * sqrt( cos(pi*t).^2 + sin(pi*t).^2 / 4 );
N = (4:4:48)';
perim = zeros(size(N));
for k = 1:length(N)
    h = 2 / N(k);
    t = h * (0:N(k)-1);
    perim(k) = h * sum(f(t));
end
err = abs(perim - perim(end));    % use last value as "exact"
format long
disp(table(N, perim, err, variableNames=["number of nodes", "perimeter", "error"]))
```
The approximations gain about one digit of accuracy for each constant increment of $n$, which is consistent with spectral convergence.
