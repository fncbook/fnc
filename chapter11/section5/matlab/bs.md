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
[**Demo %s**](#demo-boundaries-bs)


```{code-cell}
K = 3;  sigma = 0.06;  r = 0.08;  Smax = 8;
phi = @(t, x, u, ux, uxx) sigma.^2/2 * (x.^2 .* uxx) + r * x.*ux - r*u;
ga = @(u, ux) u;
gb = @(u, ux) ux - 1;
```

```{code-cell}
init = @(x) max(0, x - K);
[x, u] = parabolic(phi, [0, Smax], 80, ga, gb, [0, 15], init);
```

```{code-cell}
:tags: [remove-input]
clf
plot(x, u(0))
hold on,  grid on
axis([0, Smax, -0.1, 8])
title("Black–Scholes equation with boundaries")
xlabel('x'),  ylabel('u(x,t)')
vid = VideoWriter("figures/boundaries-bs.mp4", "MPEG-4");
vid.Quality = 85;
open(vid);
for t = linspace(0, 15, 151)
    cla, plot(x, u(t))
    str = sprintf("t = %.1f", t);
    text(0.5, 7.1, str);
    writeVideo(vid, frame2im(getframe(gcf)));
end
close(vid)
```

Recall that $u$ is the value of the call option, and time runs backward from the strike time. The longer the horizon, the more value the option has due to anticipated growth in the stock price.
