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
[**Demo %s**](#demo-blackscholes-solve)

We consider the Black–Scholes problem for the following parameter values:

```{code-cell}
Smax = 8;  T = 6;
K = 3;  sigma = 0.06;  r = 0.08;
```

We discretize space and time.

```{code-cell}
m = 200;  h = Smax / m;
x = h*(0:m)';
n = 1000;  tau = T / n;
t = tau*(0:n)';
lambda = tau / h^2;  mu = tau / h;
```

We set the initial condition and then march forward in time.

```{code-cell}
V = zeros(m+1,n+1);
V(:,1) = max( 0, x-K );
for j = 1:n
    % Fictitious value from Neumann condition.
    Vfict = 2*h + V(m,j);
    Vj = [ V(:,j); Vfict ];
    % First row is zero by the Dirichlet condition.
    for i = 2:m+1 
        diff1 = (Vj(i+1) - Vj(i-1));
        diff2 = (Vj(i+1) - 2*Vj(i) + Vj(i-1));
        V(i,j+1) = Vj(i) ...
            + (lambda * sigma^2* x(i)^2/2) * diff2  ...
            + (r*mu * x(i))/2 * diff1 - r*tau * Vj(i);
    end 
end
```

Here is a plot of the solution after every 250 time steps.

```{code-cell}
index_times = 1 + 250*(0:4);
show_times = t(index_times);
clf
for j = index_times
    str = sprintf("t = %.2f", t(j));
    plot(x, V(:, j), displayname=str) 
    hold on
end
title('Black-Scholes solution')
xlabel('stock price'),  ylabel('option value')
axis tight,  grid on
legend(location="northwest")
```

```{index} ! MATLAB; animation
```

Alternatively, here is an animation of the solution.

```{code-cell}
:tags: [remove-output]
clf
plot(x, V(:,1))
hold on,  grid on
axis([0, 8, 0, 6])
title('Black-Scholes solution') 
xlabel('stock price'),  ylabel('option value')
vid = VideoWriter("figures/black-scholes-6.mp4","MPEG-4");
vid.Quality = 85;
open(vid)
for frame = 1:10:n+1
    cla, plot(x, V(:, frame))
    str = sprintf("t = %.2f", t(frame));
    text(0.4, 5.2, str);
    writeVideo(vid, frame2im(getframe(gcf)));
end
close(vid)
```


The results are easy to interpret, recalling that the time variable really means *time until strike*. Say you are close to the option's strike time. If the current stock price is, say, $S=2$, then it's not likely that the stock will end up over the strike price $K=3$, and therefore the option has little value. On the other hand, if presently $S=3$, then there are good odds that the option will be exercised at the strike time, and you will need to pay a substantial portion of the stock price in order to take advantage. As the time to strike increases, there is an expectation that the stock price is more likely to rise somewhat, making the value of the option larger at each fixed $S$. 
