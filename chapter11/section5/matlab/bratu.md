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
[**Demo %s**](#demo-boundaries-bratu)


```{code-cell}
phi = @(t, x, u, ux, uxx) u.^2 + uxx;
ga = @(u, ux) u;
gb = @(u, ux) ux;

init = @(x) 400 * x.^4 .* (1 - x).^2;
[x, u] = parabolic(phi, [0, 1], 60, ga, gb, [0, 0.1], init);
```

```{code-cell}
:tags: [hide-cell]
clf
plot(x, u(0))
hold on,  grid on
axis([0, 1, 0, 10])
title("Heat equation with source")
xlabel('x'),  ylabel('u(x,t)')
vid = VideoWriter("figures/boundaries-source.mp4", "MPEG-4");
vid.Quality = 85;
open(vid);
for t = linspace(0, 0.1, 101)
    cla, plot(x, u(t))
    str = sprintf("t = %.3f", t);
    text(0.05, 9.2, str);
    writeVideo(vid, frame2im(getframe(gcf)));
end
close(vid) 
```

