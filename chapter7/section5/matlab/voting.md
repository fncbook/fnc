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
[**Demo %s**](#demo-dimreduce-voting)

This matrix describes the votes on bills in the 111th session of the United States Senate. (The data set was obtained from [https://voteview.com].) Each row is one senator, and each column is a vote item.

```{code-cell}
load voting
```

If we visualize the votes (yellow is "yea," blue is "nay"), we can see great similarity between many rows, reflecting party unity.

```{code-cell}
clf
imagesc(A)
colormap parula
title('Votes in 111th U.S. Senate')
ylabel(('senator'),  xlabel('bill'));
```

We use {eq}`sing-val-decay` to quantify the decay rate of the values.

```{code-cell}
[U, S, V] = svd(A);
sigma = diag(S);
tau = cumsum(sigma.^2) / sum(sigma.^2);
plot(tau(1:16), 'o')
xlabel('k'),  ylabel('\tau_k')
title(('Fraction of singular value energy'));
```

The first and second singular triples contain about 58% and 17%, respectively, of the energy of the matrix. All others have far less effect, suggesting that the information is primarily two-dimensional. The first left and right singular vectors also contain interesting structure.

```{code-cell}
subplot(211), plot(U(:, 1), '.')
xlabel('senator number'), title('left singular vector')
subplot(212), plot(V(:, 1), '.')
xlabel('bill number'), title(('right singular vector'));
```

Both vectors have values greatly clustered near $\pm C$ for a constant $C$. These can be roughly interpreted as how partisan a particular senator or bill was, and for which political party. Projecting the senators' vectors into the first two $\mathbf{V}$-coordinates gives a particularly nice way to reduce them to two dimensions. Political scientists label these dimensions *partisanship* and *bipartisanship*. Here we color them by actual party affiliation (also given in the data file): red for Republican, blue for Democrat, and yellow for independent.

```{code-cell}
clf
x1 = V(:, 1)'*A';   x2 = V(:, 2)'*A'; 
scatter(x1(Dem), x2(Dem), 20, 'b'),  hold on
scatter(x1(Rep), x2(Rep), 20, 'r')
scatter(x1(Ind), x2(Ind), 20, 'k')
xlabel('partisanship'),  ylabel('bipartisanship')
title(('111th US Senate in 2D'));
```
