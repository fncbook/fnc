---
numbering:
  enumerator: 6.7.%s
---
(section-ivp-implicit)=
# Implementation of multistep methods

```{index} multistep method; implementation of
```

We now consider some of the practical issues that arise when multistep formulas are used to solve IVPs. In this section we emphasize the vector IVP, $\mathbf{u}'=\mathbf{f}(t,\mathbf{u})$, and use boldface in the difference formula {eq}`multistep` as well. 

## Explicit methods

As a concrete example, the AB4 method is defined by the formula

```{math}
:label: ab4
\mathbf{u}_{i+1} = \mathbf{u}_i + h\, ( \tfrac{55}{24}\mathbf{f}_i - \tfrac{59}{24} \mathbf{f}_{i-1} + \tfrac{37}{24}\mathbf{f}_{i-2} - \tfrac{9}{24}\mathbf{f}_{i-3}), \quad i=3,\ldots,n-1.
```

{numref}`Function {number} <function-ab4>` shows a basic implementation of AB4. 

Observe that {numref}`Function {number} <function-rk4>` is used to find the starting values $\mathbf{u}_1,\mathbf{u}_2,\mathbf{u}_3$ that are needed before the iteration formula takes over. As far as RK4 is concerned, it needs to solve  (the same step size as in the AB4 iteration). These results are then used to find $\mathbf{f}_0,\ldots,\mathbf{f}_3$ and get the main iteration started.

(function-ab4)=
``````{prf:algorithm} ab4
`````{tab-set} 
````{tab-item} Julia
:sync: julia
:::{embed} #function-ab4-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #function-ab4-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #function-ab4-python
:::
````
`````
``````

(demo-implicit-ab4)=
```````{prf:example}
``````{tab-set}
`````{tab-item} Julia
:sync: julia
````{dropdown} Convergence of Adams–Bashforth
:open:
```{include} julia/ab4.ipynb
```
````
`````

`````{tab-item} MATLAB
:sync: matlab
````{dropdown} Convergence of Adams–Bashforth
:open:
```{include} matlab/ab4.ipynb
```
````
`````

`````{tab-item} Python
:sync: python
````{dropdown} Convergence of Adams–Bashforth
:open:
```{include} python/ab4.ipynb
```
````
`````
``````
```````
    

## Implicit methods

```{index} implicit IVP solver
```

The implementation of an implicit multistep method is more difficult. Consider the second-order implicit formula AM2, also known as the trapezoid method. To advance from step $i$ to $i+1$, we need to solve

```{math}
:label: AM2solve
  \mathbf{z} - \tfrac{1}{2} h f(t_{i+1},\mathbf{z})  = \mathbf{u}_i + \tfrac{1}{2} h \mathbf{f}(t_i,\mathbf{u}_i)
```

for $\mathbf{z}$. This equation can be written as $\mathbf{g}(\mathbf{z})=\boldsymbol{0}$, so the rootfinding methods of Chapter 4 can be used. The new value $\mathbf{u}_{i+1}$ is equal to the root of this equation.  

An implementation of AM2 using {numref}`Function {number} <function-levenberg>` from {numref}`section-nonlineqn-quasinewton` is shown in {numref}`Function {number} <function-am2>`. 

(function-am2)=
``````{prf:algorithm} am2
`````{tab-set} 
````{tab-item} Julia
:sync: julia
:::{embed} #function-am2-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #function-am2-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #function-am2-python
:::
````
`````
``````

## Stiff problems

At each time step in {numref}`Function {number} <function-am2>`, or any implicit IVP solver, a rootfinding iteration of uncertain expense is needed, requiring multiple calls to evaluate the function $\mathbf{f}$. This fact makes the cost of an implicit method much greater on a per-step basis than for an explicit one. Given this drawback, you are justified to wonder whether implicit methods are ever competitive! The answer is emphatically yes, as {numref}`Demo %s <demo-implicit-stiff>` demonstrates.

```{index} stiff differential equation
```

(demo-implicit-stiff)=
```````{prf:example}
``````{tab-set}
`````{tab-item} Julia
:sync: julia
````{dropdown} Stiffness
:open:
```{include} julia/stiff.ipynb
```
````
`````

`````{tab-item} MATLAB
:sync: matlab
````{dropdown} Stiffness
:open:
```{include} matlab/stiff.ipynb
```
````
`````

`````{tab-item} Python
:sync: python
````{dropdown} Stiffness
:open:
```{include} python/stiff.ipynb
```
````
`````
``````
```````
    

Although the result of {numref}`Demo %s <demo-implicit-stiff>` may seem counter-intuitive, there is no contradiction. A fourth-order explicit formula is more accurate than a second-order implicit one, in the limit $h\to 0$. But there is another limit to consider, $t\to \infty$ with $h$ fixed, and in this one the implicit method wins. 

Problems for which implicit methods are much more efficient than explicit counterparts are called **stiff**. A complete mathematical description will wait for Chapter 11, but a sure sign of stiffness is the presence of phenomena on widely different time scales. In {numref}`Demo {number} <demo-implicit-stiff>`, for instance, there are two slow periods during which the solution changes very little, interrupted by a very fast transition in the state. An explicit method "thinks" that the step size must always be dictated by the time scale of the fast transition, whereas an implicit method can take large steps during the slow periods.

```{index} adaptivity; in IVP solver
```

## Adaptivity

As with RK methods, we can run two time stepping methods simultaneously in order to estimate the error and adjust the step size. For example, we could pair AB3 with AB4 as practically no cost because the methods differ only in how they include known information from the recent past. The more accurate AB4 value should allow an accurate estimate of the local error in the AB3 value, and so on.

Because multistep methods rely on the solution history, though, changing the step size is more algebraically complicated than for RK methods. If $h$ is changed, then the historical values $\mathbf{u}_{i-1},\mathbf{u}_{i-2}\ldots$ and $\mathbf{f}_{i-1},\mathbf{f}_{i-2}\ldots$ are no longer given at the right moments in time to apply the iteration formula. A typical remedy is to use interpolation to re-evaluate the historical values at the appropriate times. The details are important but not especially illuminating, and we do not give them here.

## Exercises

% Must stay as #1
(problem-implicit-ab4tests)=

1. ⌨ For each IVP, solve the problem using {numref}`Function {number} <function-ab4>` with $n=100$, and plot the solution and the error $u-\hat{u}$ on separate plots.

    **(a)** $u' = -2t u, \ 0 \le t \le 2, \ u(0) = 2;\  \hat{u}(t) = 2e^{-t^2}$

    **(b)** $u' = u + t, \ 0 \le t \le 1, \ u(0) = 2;\  \hat{u}(t) = 1-t+e^t$

    **(c)** $u' = x^2/[u(1+x^3)],\ 0 \le x \le 3, \ u(0) =1;\ \hat{u}(x) =[1+(2/3)\ln (1+x^3)]^{1/2}$

    **(d)** $u''+ 9u = 9t, \: 0< t< 2\pi, \: u(0) =1,\: u'(0) = 1; \: \hat{u}(t) = t+\cos (3t)$

    **(e)** $u''+ 9u = \sin(2t), \: 0< t< 2\pi, \: u(0) =2,\: u'(0) = 1$;
    $\quad \hat{u}(t) = (1/5) \sin(3t) + 2 \cos (3t)+ (1/5) \sin (2t)$

    **(f)** $u''- 9u = 9t \: 0< t< 1, \: u(0) =2,\: u'(0) = -1; \: \hat{u}(t) = e^{3t} + e^{-3t}-t$

    **(g)** $u''+ 4u'+ 4u = t, \: 0< t< 4, \: u(0) =1,\: u'(0) = 3/4; \: \hat{u}(t) = (3t+5/4)e^{-2t} + (t-1)/4$

    **(h)** $x^2 u'' +5xu' + 4u = 0,\: 1<x<e^2, \: u(1) =1, \: u'(1) = -1; \: \hat{u}(x) = x^{-2}( 1 + \ln x)$

    **(i)** $2 x^2 u'' +3xu' - u = 0,\: 1<x<16, \: u(1) =4, \: u'(1) = -1$;
    $\quad \hat{u}(x) = 2(x^{1/2} + x^{-1})$

    **(j)** $x^2 u'' -xu' + 2u = 0,\: 1<x<e^{\pi}, \: u(1) =3, \: u'(1) = 4$;
    $\quad \hat{u}(x) = x \left[ 3 \cos \left( \ln x \right)+\sin \left( \ln x \right) \right]$

(problem-implicit-ab4converge)=
2. ⌨ For each IVP in Exercise 1, use {numref}`Function {number} <function-ab4>` for $n=10\cdot2^d$ and $d=1,\ldots,10$. Make a log-log convergence plot for the final time error $|u_n-\hat{u}(t_n)|$ versus $n$, and add a straight line indicating fourth-order convergence.

3. ⌨ Repeat Exercise 1 above  using {numref}`Function {number} <function-am2>`.

4. ⌨  Repeat Exercise 2 above using {numref}`Function {number} <function-am2>` and comparing to second-order rather than fourth-order convergence.

5. ⌨ Using {numref}`Function {number} <function-am2>` as a model, write a function `bd2` that applies the BD2 method to solve an IVP. Test the convergence of your function on one of the IVPs in [Exercise 1](#problem-implicit-ab4tests) above.

6. ⌨ For double-precision purposes, the exact solution of the IVP in {numref}`Demo %s <demo-implicit-stiff>` satisfies $\hat{u}(400)=1$.

    **(a)** Use {numref}`Function {number} <function-ab4>` with $n=600,800,1000,\ldots,2000$ and make a log-log convergence plot of the error $|u_n-1|$ as a function of $n$.

    **(b)** Repeat part (a) using {numref}`Function {number} <function-am2>`.

(problem-implicit-ivpimag)=
7. Consider the IVP
  
    ```{math}
    \mathbf{u}'(t) = \mathbf{A} \mathbf{u}(t), \quad \mathbf{A}=
    \begin{bmatrix}
      0&-4\\4&0
    \end{bmatrix}, \quad \mathbf{u}(0) =
    \begin{bmatrix}
      1\\0
    \end{bmatrix}.
    ```

    **(a)** ✍ Define $E(t) = \bigl\|\mathbf{u}(t)\bigr\|_2^2$. Show that $E(t)$ is constant. (Hint: differentiate $\mathbf{u}^T\mathbf{u}$ with respect to time and simplify it.)

    **(b)** ⌨ Use {numref}`Function {number} <function-ab4>` to solve the IVP for $t\in[0,20]$ with $n=100$ and $n=150$. Plot $|E(t)-E(0)|$ versus time for both solutions on a single log-linear graph. You should see exponential growth in time. (In this regime, AB4 is acting unstably in a sense discussed in {numref}`section-diffusion-absstab`.)

    **(c)** ⌨ Repeat part (b) with $n=400$ and $n=600$, but on a linear-linear plot. Now you should see only linear growth of $|E(t)-E(0)|$. (In this regime, AB4 is fully stable.)

    **(d)** ⌨ Repeat part (b) with AM2 instead of AB4, on a linear-linear plot. You will find that AM2 conserves energy, just like the exact solution. 
  
8. ⌨ **(a)** Modify {numref}`Function {number} <function-ab4>` to implement the AB2 method.

    **(b)** Repeat part (b) of the preceding exercise, using AB2 in place of AB4.

    **(c)** Repeat part (c) of the preceding exercise, using AB2 in place of AB4.
  
9. ⌨ **(a)** Modify {numref}`Function {number} <function-am2>` to implement the backward Euler (AM1) method.

    **(b)** Repeat part (d) of Exercise 7 above, using AM1 in place of AM2 and $n=400,800$. Does the AM1 method conserve energy?
