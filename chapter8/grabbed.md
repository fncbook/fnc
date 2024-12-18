(demo-minrescg-indefinite-julia)=
#### minrescg-indefinite





The following matrix is indefinite.

```{code-cell}
A = FNC.poisson(10) - 20I
λ = eigvals(Matrix(A))
isneg = @. λ < 0
@show sum(isneg),sum(.!isneg);
```

We can compute the relevant quantities from {numref}`Theorem {number} <theorem-minrescg-indefinite>`.

```{code-cell}
mn,mx = extrema(-λ[isneg])
κ₋ = mx/mn
mn,mx = extrema(λ[.!isneg])
κ₊ = mx/mn
ρ = (sqrt(κ₋*κ₊)-1) / (sqrt(κ₋*κ₊)+1)
```

Because the iteration number $m$ is halved in {eq}`minres-conv`, the rate constant of linear convergence is the square root of this number, which makes it even closer to 1. 

Now we apply MINRES to a linear system with this matrix, and compare the observed convergence to the upper bound from the theorem.

```{index} ! Julia; minres
```

```{code-cell}
b = rand(100)
x,hist = minres(A,b,reltol=1e-10,maxiter=51,log=true);

relres = hist[:resnorm] / norm(b)
m = 0:length(relres)-1
plot(m,relres,label="observed",leg=:left,
	xaxis=L"m",yaxis=(:log10,"relative residual"),
	title=("Convergence of MINRES") )
plot!(m,ρ.^(m/2),l=:dash,label="upper bound")
```

The upper bound turns out to be pessimistic here, especially in the later iterations. However, you can see a slow linear phase in the convergence that tracks the bound closely.
(demo-minrescg-converge-julia)=
#### minrescg-converge





We will compare MINRES and CG on some quasi-random SPD problems.  The first matrix has a condition number of 100. 

```{code-cell}
n = 5000
density = 0.001
A = FNC.sprandsym(n,density,1/100)
x = (1:n)/n
b = A*x;
```

```{index} ! Julia; cg
```

Now we apply both methods and compare the convergence of the system residuals, using implementations imported from `IterativeSolvers`.

```{code-cell}
plt = plot(title="Convergence of MINRES and CG",
    xaxis=("Krylov dimension"),yaxis=(:log10,"relative residual norm"))
for method in [minres,cg]
    x̃,history = method(A,b,reltol=1e-6,maxiter=1000,log=true);
    relres = history[:resnorm] / norm(b)
    plot!(0:length(relres)-1,relres,label="$method")
    err = round( norm( x̃ - x ) / norm(x), sigdigits=4 )
    println("$method error: $err")
end
plt
```

There is little difference between the two methods here. Next, we increase the condition number of the matrix by a factor of 25. The rule of thumb predicts that the number of iterations required should increase by a factor of about 5.

```{code-cell}
A = FNC.sprandsym(n,density,1/2500)
b = A*x;
```

```{code-cell}
:tags: [hide-input]

plt = plot(title="Convergence of MINRES and CG",
    xaxis=("Krylov dimension"),yaxis=(:log10,"relative residual norm"))
for method in [minres,cg]
    x̃,history = method(A,b,reltol=1e-6,maxiter=1000,log=true);
    relres = history[:resnorm] / norm(b)
    plot!(0:length(relres)-1,relres,label="$method")
    err = round( norm( x̃ - x ) / norm(x), sigdigits=4 )
    println("$method error: $err")
end
plt
```

Both methods have an early superlinear phase that allow them to finish slightly sooner than the factor of 5 predicted: {numref}`Theorem {number} <theorem-minrescg-converge>` is an upper bound, not necessarily an approximation. Both methods ultimately achieve the same reduction in the residual; MINRES stops earlier, but with a slightly larger error.

(demo-power-one-julia)=
#### power-one





Here we choose a random 5×5 matrix and a random 5-vector.

```{code-cell}
A = rand(1.:9.,5,5)
A = A./sum(A,dims=1)
x = randn(5)
```

Applying matrix-vector multiplication once doesn't do anything recognizable.

```{code-cell}
y = A*x
```

Repeating the multiplication still doesn't do anything obvious.

```{code-cell}
z = A*y
```

But if we keep repeating the matrix-vector multiplication, something remarkable happens: $\mathbf{A} \mathbf{x} \approx \mathbf{x}$.

```{code-cell}
for j in 1:8;  x = A*x;  end
[x A*x]
```

This phenomenon seems to occur regardless of the starting vector.

```{code-cell}
x = randn(5)
for j in 1:8;  x = A*x;  end
[x A*x]
```
(demo-power-iter-julia)=
#### power-iter

We will experiment with the power iteration on a 5×5 matrix with prescribed eigenvalues and dominant eigenvalue at 1.

```{code-cell}
λ = [1,-0.75,0.6,-0.4,0]
# Make a triangular matrix with eigenvalues on the diagonal.
A = triu(ones(5,5),1) + diagm(λ)
```

We run the power iteration 60 times. The best estimate of the dominant eigenvalue is the last entry of the first output.

```{code-cell}
β,x = FNC.poweriter(A,60)
eigval = β[end]
```

We check for linear convergence using a log-linear plot of the error.

```{code-cell}
err = @. 1 - β
plot(0:59,abs.(err),m=:o,title="Convergence of power iteration",
    xlabel=L"k",yaxis=(L"|\lambda_1-\beta_k|",:log10,[1e-10,1]))
```

The asymptotic trend seems to be a straight line, consistent with linear convergence. To estimate the convergence rate, we look at the ratio of two consecutive errors in the linear part of the convergence curve. The ratio of the first two eigenvalues should match the observed rate.  

```{code-cell}
@show theory = λ[2]/λ[1];
@show observed = err[40]/err[39];
```

Note that the error is supposed to change sign on each iteration. The effect of these alternating signs is that estimates oscillate around the exact value.

```{code-cell}
β[26:30]
```

In practical situations, we don't know the exact eigenvalue that the algorithm is supposed to find. In that case we would base errors on the final $\beta$ that was found, as in the following plot.

```{code-cell}
err = @. β[end] - β[1:end-1]
plot(0:58,abs.(err),m=:o,title="Convergence of power iteration",
    xlabel=L"k",yaxis=(L"|\beta_{60}-\beta_k|",:log10,[1e-10,1]))
```

The results are very similar until the last few iterations, when the limited accuracy of the reference value begins to show. That is, while it is a good estimate of $\lambda_1$, it is less good as an estimate of the error in nearby estimates.
(demo-precond-diagonal-julia)=
#### precond-diagonal

Here is an SPD matrix that arises from solving partial differential equations.

```{code-cell}
A = matrixdepot("wathen",60)
n = size(A,1)
@show n,nnz(A);
```

```{index} ! Julia; DiagonalPreconditioner
```

There is an easy way to use the diagonal elements of $\mathbf{A}$, or any other vector, as a diagonal preconditioner.

```{code-cell}
b = ones(n)
M = DiagonalPreconditioner(diag(A));
```

We now compare CG with and without the preconditioner.

```{code-cell}
:tags: [hide-input]
plain(b) = cg(A,b,maxiter=200,reltol=1e-4,log=true)
time_plain = @elapsed x,hist1 = plain(b)
prec(b) = cg(A,b,Pl=M,maxiter=200,reltol=1e-4,log=true)
time_prec = @elapsed x,hist2 = prec(b)
@show time_plain,time_prec

rr = hist1[:resnorm]
plot(0:length(rr)-1,rr/rr[1],yscale=:log10,label="plain")
rr = hist2[:resnorm]
plot!(0:length(rr)-1,rr/rr[1],yscale=:log10,label="preconditioned")
title!("Diagonal preconditioning in CG")
```

The diagonal preconditioner cut down substantially on the number of iterations. The effect on the total time is less dramatic, but this is not a large version of the problem.
(demo-precond-gmres-julia)=
#### precond-gmres

Here is a nonsymmetric matrix arising from a probabilistic model in computational chemistry.

```{code-cell}
A = sparse(matrixdepot("Watson/chem_master1"))
n = size(A,1)
@show n,nnz(A),issymmetric(A)
```

Without a preconditioner, GMRES makes essentially no progress after 100 iterations. 

```{code-cell}
b = rand(40000)
const GMRES = IterativeSolvers.gmres
x,history = GMRES(A,b,maxiter=100,reltol=1e-5,log=true)
resnorm = history[:resnorm]
@show resnorm[end] / resnorm[1];
```

```{index} ! Julia; ilu
```

The following version of incomplete LU factorization drops all sufficiently small values (i.e., replaces them with zeros). This keeps the number of nonzeros in the factors under control. 

```{code-cell}
iLU = ilu(A,τ=0.25)
@show nnz(iLU)/nnz(A);
```

The result is almost 10 times as dense as $\mathbf{A}$ and yet still not a true factorization of it. However, it's close enough for an approximate inverse in a preconditioner. The actual preconditioning matrix is $\mathbf{M}=\mathbf{L}\mathbf{U}$, but we just supply the factorization to `gmres`.

```{code-cell}
_,history = GMRES(A,b,Pl=iLU,maxiter=100,reltol=1e-5,log=true)
history
```

The $\tau$ parameter in `ilu` balances the accuracy of the iLU factorization with the time needed to compute it and invert it. As $\tau\to 0$, more of the elements are kept, making the preconditioner more effective but slower per iteration. 

```{code-cell}
:tags: [hide-input]

plt = plot(0:40,resnorm[1:41]/resnorm[1],label="no preconditioning",
        xaxis=("iteration number"),yaxis=(:log10,"residual norm"),
        leg=:bottomright,title="Incomplete LU preconditioning")
    for τ in [2,1,0.25,0.1]
        t = @elapsed iLU = ilu(A;τ)
        t += @elapsed _,history = GMRES(A,b,Pl=iLU,maxiter=100,
                                        reltol=1e-5,log=true)
        resnorm = history[:resnorm]
        label = "τ = $τ, time = $(round(t,digits=3))"
        plot!(0:length(resnorm)-1,resnorm/resnorm[1];label)
    end
plt
```

In any given problem, it's impossible to know in advance where the right balance lies between fidelity and speed for the preconditioner.
(demo-structure-sparse-julia)=
#### structure-sparse

Here we load the adjacency matrix of a graph with 2790 nodes. Each node is a web page referring to Roswell, NM, and the edges represent links between web pages. (Credit goes to Panayiotis Tsaparas and the University of Toronto for making this data public.)

```{code-cell}
@load "roswell.jld2" A;      # file is on the book's website
```

```{index} ! Julia; nnz
```
::::{grid} 1 1 2 2

:::{grid-item}


We may define the density of $\mathbf{A}$ as the number of nonzeros divided by the total number of entries.


:::
:::{card}


Use `nnz` to count the number of nonzeros in a sparse matrix.

:::
::::

```{code-cell}
m,n = size(A)
@show density = nnz(A) / (m*n);
```

```{index} ! Julia; summarysize
```

The computer memory consumed by any variable can be discovered using `summarysize`. We can use it to compare the space needed for the sparse representation to its dense counterpart, that is, the space needed to store all the elements, whether zero or not.

```{code-cell}
F = Matrix(A)
Base.summarysize(F) / Base.summarysize(A)
```

As you can see, the storage savings are dramatic. Matrix-vector products are also much faster using the sparse form because operations with structural zeros are skipped.

```{code-cell}
x = randn(n)
A*x;   # make sure * is loaded and compiled
@elapsed for i in 1:300; A*x; end
```

```{code-cell}
F*x;
@elapsed for i in 1:300; F*x; end
```
(demo-structure-fill-julia)=
#### structure-fill

Here is the adjacency matrix of a graph representing a small-world network, featuring connections to neighbors and a small number of distant contacts.

```{code-cell}
@load "smallworld.jld2" A
graphplot(A,linealpha=0.5)
```

Because each node connects to relatively few others, the adjacency matrix is quite sparse.

```{code-cell}
spy(A,title="Nonzero locations",m=2,color=:blues)
```

By {numref}`Theorem {number} <theorem-insight-adjmat>`, the entries of $\mathbf{A}^k$ give the number of walks of length $k$ between pairs of nodes, as with "*k* degrees of separation" within a social network. As $k$ grows, the density of $\mathbf{A}^k$ also grows.

```{code-cell}
plt = plot(layout=(1,3),legend=:none,size=(600,240))
for k in 2:4
    spy!(A^k,subplot=k-1,color=:blues,
        title=latexstring("\\mathbf{A}^$k"))
end
plt
```
(demo-structure-banded-julia)=
#### structure-banded

```{index} ! Julia; spdiagm
```

The `spdiagm` function creates a sparse matrix given its diagonal elements. The main or central diagonal is numbered zero, above and to the right of that is positive, and below and to the left is negative.

```{code-cell}
n = 50;
A = spdiagm(-3=>fill(n,n-3),
             0=>ones(n),
             1=>-(1:n-1),
             5=>fill(0.1,n-5) )
Matrix(A[1:7,1:7])
```

```{index} ! Julia; sparse
```

::::{grid} 1 1 2 2

:::{grid-item}


Without pivoting, the LU factors have the same lower and upper bandwidth as the orignal matrix.


:::
:::{card}


The `sparse` function converts any matrix to sparse form. But it's usually better to construct a sparse matrix directly, as the standard form might not fit in memory.

:::
::::

```{code-cell}
L,U = FNC.lufact(A)
plot(layout=2)
spy!(sparse(L),m=2,subplot=1,title=L"\mathbf{L}",color=:blues)
spy!(sparse(U),m=2,subplot=2,title=L"\mathbf{U}",color=:blues)
```

However, if we introduce row pivoting, bandedness may be expanded or destroyed.

```{code-cell}
fact = lu(A)
plot(layout=2)
spy!(sparse(fact.L),m=2,subplot=1,title=L"\mathbf{L}",color=:blues)
spy!(sparse(fact.U),m=2,subplot=2,title=L"\mathbf{U}",color=:blues)
```
(demo-structure-linalg-julia)=
#### structure-linalg
The following generates a random sparse matrix with prescribed eigenvalues.

```{code-cell}
n = 4000
density = 4e-4
λ = @. 1 + 1/(1:n)   # exact eigenvalues
A = FNC.sprandsym(n,density,λ);
```

```{index} ! Julia; eigs
```

The `eigs` function finds a small number eigenvalues meeting some criterion. First, we ask for the 5 of largest (complex) magnitude using `which=:LM`.

```{code-cell}
λmax,V = eigs(A,nev=5,which=:LM)    # Largest Magnitude
fmt = ft_printf("%20.15f")
pretty_table([λmax λ[1:5]], header=["found","exact"], formatters=fmt)
```

Now we find the 5 closest to the value 1 in the complex plane, via `sigma=1`.

```{code-cell}
λ1,V = eigs(A,nev=5,sigma=1)    # closest to sigma
data = [λ1 λ[end:-1:end-4]]
pretty_table(data, header=["found","exact"], formatters=fmt)
```

```{index} Julia; \\
```

The time needed to solve a sparse linear system is not easy to predict unless you have some more information about the matrix. But it will typically be orders of magnitude faster than the dense version of the same problem.

```{code-cell}
x = @. 1/(1:n);  b = A*x;
```

```{code-cell}
norm(x-A\b);  # force compilation
t = @elapsed sparse_err = norm(x - A\b)
println("Time for sparse solve: $t")
```

```{code-cell}
D = Matrix(A)  # convert to regular matrix
norm(x-D\b);
t = @elapsed dense_err = norm(x - D\b)
println("Time for dense solve: $t")
```

```{code-cell}
@show sparse_err;
@show dense_err;
```
(demo-subspace-unstable-julia)=
#### subspace-unstable

First we define a triangular matrix with known eigenvalues, and a random vector $b$.

```{code-cell}
λ = @. 10 + (1:100)
A = triu(rand(100,100),1) + diagm(λ)
b = rand(100);
```

Next we build up the first ten Krylov matrices iteratively, using renormalization after each matrix-vector multiplication. 

```{code-cell}
Km = [b zeros(100,29)]
for m in 1:29      
    v = A*Km[:,m]
    Km[:,m+1] = v/norm(v)
end
```

Now we solve least-squares problems for Krylov matrices of increasing dimension, recording the residual in each case.

```{code-cell}
resid = zeros(30)
for m in 1:30  
    z = (A*Km[:,1:m])\b
    x = Km[:,1:m]*z
    resid[m] = norm(b-A*x)
end
```

The linear system approximations show smooth linear convergence at first, but the convergence stagnates after only a few digits have been found.

```{code-cell}
plot(0:29,resid,m=:o,
    xaxis=(L"m"),yaxis=(:log10,L"\| b-Ax_m \|"), 
    title="Residual for linear systems",leg=:none)
```
(demo-subspace-arnoldi-julia)=
#### subspace-arnoldi
We illustrate a few steps of the Arnoldi iteration for a small matrix.

```{code-cell}
A = rand(1.:9.,6,6)
```

The seed vector we choose here determines the first member of the orthonormal basis.

```{code-cell}
u = randn(6)
Q = u/norm(u);
```

Multiplication by $\mathbf{A}$ gives us a new vector in $\mathcal{K}_2$. 

```{code-cell}
Aq = A*Q[:,1];
```

We subtract off its projection in the previous direction. The remainder is rescaled to give us the next orthonormal column.

```{code-cell}
v = Aq - dot(Q[:,1],Aq)*Q[:,1]
Q = [Q v/norm(v)];
```

On the next pass, we have to subtract off the projections in two previous directions.

```{code-cell}
Aq = A*Q[:,2]
v = Aq - dot(Q[:,1],Aq)*Q[:,1] - dot(Q[:,2],Aq)*Q[:,2]
Q = [Q v/norm(v)];
```

At every step, $\mathbf{Q}_m$ is an ONC matrix.

```{code-cell}
@show opnorm( Q'*Q - I );
```

And $\mathbf{Q}_m$ spans the same space as the three-dimensional Krylov matrix.

```{code-cell}
K = [ u A*u A*A*u ];
@show rank( [Q K] );
```
