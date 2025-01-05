---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
---
We consider the Black–Scholes problem for the following parameter values:

```{code-cell}
Smax = 8 
T = 6
K, σ, r = (3, 0.06, 0.08);
```

We discretize space and time.

```{code-cell}
m = 200;  h = Smax / m;
x = h * (0:m)
n = 1000;  τ = T / n;
t = τ * (0:n)
λ = τ / h^2
μ = τ / h
```

We set the initial condition and then march forward in time.

```{code-cell}
V = zeros(m+1, n+1)
V[:, 1] = @. max(0, x - K)
for j in 1:n
    # Fictitious value from Neumann condition.
    Vfict = 2*h + V[m,j]
    Vj = [ V[:, j]; Vfict ]
    # First row is zero by the Dirichlet condition.
    for i in 2:m+1 
        diff1 = (Vj[i+1] - Vj[i-1])
        diff2 = (Vj[i+1] - 2Vj[i] + Vj[i-1])
        V[i,j+1] = Vj[i] +
            (λ * σ^2 * x[i]^2 / 2) * diff2 
            + (r * x[i] * μ) / 2 * diff1 
            - (r * τ) * Vj[i]
    end 
end
```

Here is a plot of the solution after every 250 time steps.

```{code-cell}
using Plots
idx = 1:250:n+1
label = reshape(["t = $t" for t in t[idx]], 1, length(idx))
plot(x, V[:, idx]; 
    label, legend=:topleft,
    xaxis=("stock price"),  yaxis=("option value"),
    title="Black–Scholes solution")
```

```{index} ! Julia; @animate
```

Alternatively, here is an animation of the solution.

```{code-cell}
anim = @animate for j in 1:10:n+1
    plot(x, V[:, j];
        xaxis=(L"S"),  yaxis=([0,6],L"v(S,t)"),
        title="Black–Scholes solution",
        dpi=150,    
        label=@sprintf("t = %.2f", t[j]))
end
mp4(anim, "figures/black-scholes-6.mp4")
```

The results are easy to interpret, recalling that the time variable really means *time until strike*. Say you are close to the option's strike time. If the current stock price is, say, $S=2$, then it's not likely that the stock will end up over the strike price $K=3$, and therefore the option has little value. On the other hand, if presently $S=3$, then there are good odds that the option will be exercised at the strike time, and you will need to pay a substantial portion of the stock price in order to take advantage. As the time to strike increases, there is an expectation that the stock price is more likely to rise somewhat, making the value of the option larger at each fixed $S$. 
