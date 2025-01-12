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
[**Demo %s**](#demo-blackscholes-unstable)

Let's try to do everything the same as in {numref}`Demo {number} <demo-blackscholes-solve>`, but extending the simulation time to $T=8$.

```{code-cell}
T = 8;
n = 1000;  tau = T / n;
t = tau*(0:n)';
lambda = tau / h^2;  mu = tau / h;
for j = 1:n
    % Fictitious value from Neumann condition.
    Vfict = 2*h + V(m,j);
    Vj = [ V(:,j); Vfict ];
    % First row is zero by the Dirichlet condition.
    for i = 2:m+1 
        diff1 = (Vj(i+1) - Vj(i-1));
        diff2 = (Vj(i+1) - 2*Vj(i) + Vj(i-1));
        V(i,j+1) = Vj(i) ...
             + (lambda*sigma^2 * x(i)^2/2) * diff2  ...
             + (r*mu * x(i))/2 * diff1 - r*tau * Vj(i);
    end   
end
clf
for j = index_times
    str = sprintf("t = %.2f", t(j));
    plot(x, V(:, j), displayname=str) 
    hold on
end
title('Black-Scholes instability')
xlabel('stock price'),  ylabel('option value')
axis tight,  grid on
legend(location="northwest")
```

```{code-cell}
:tags: [remove-output]
clf
plot(x, V(:,1))
hold on,  grid on
axis([0, 8, 0, 6])
title('Black-Scholes solution...?') 
xlabel('stock price'),  ylabel('option value')
vid = VideoWriter("figures/black-scholes-8.mp4","MPEG-4");
vid.Quality = 85;
open(vid);
for frame = 1:10:n+1
    cla, plot(x, V(:, frame))
    str = sprintf("t = %.2f", t(frame));
    text(0.4, 5.2, str);
    writeVideo(vid, frame2im(getframe(gcf)));
end
close(vid)
```


This so-called solution is nonsense!
