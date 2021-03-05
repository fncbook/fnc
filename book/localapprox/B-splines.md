%# B-splines
%\label{sec:bsplines}
%
%
%In section \ref{sec:pwlinear} we expressed a piecewise linear
%interpolant as $\sum c_k H_k(x)$, using the hat function basis. This
%basis has many attractive features, but its second order accuracy is
%mediocre, and the resulting interpolant lacks a continuous derivative.
%By changing the basis functions we can improve on both counts—at
%a cost.
%
%%  in two equivalent ways: first, from the point of view of a
%% linear function on each subinterval, and second, from the point of
%% view of a basis of globally defined functions. Each type of
%% representation has certain strengths. In this section we will
%% complement the piecewise-polynomial standpoint of
%% section~\ref{sec:splines} for cubic splines with the basis point of
%% view. The most natural basis for cubic splines are known as
%% **cubic B-splines**.
%
%We begin with the definition
%
```{math}
%  :label: bspline1
%  B(x;s,t) =
%  \begin{cases}
%    1,& \text{if $s \le x < t$}, \\
%    0, & \text{otherwise}.
%  \end{cases}
%```

%This function is also known as a characteristic or \emph{indicator
%function}, as it simply indicates whether $x$ is in the interval
%$[s,t)$. It is piecewise constant, with jumps in value at $x=s$ and
%$x=t$.
%
%Close examination (and a look back at {eq}`hatfun`) should
%convince you that we can write a piecewise linear hat function as
%
```{math}
%  H_k(x) = \frac{x-s_{k}}{s_{k+1}-s_{k}} B(x;s_{k},s_{k+1}) +
%  \frac{s_{k+2}-x}{s_{k+2}-s_{k+1}} B(x;s_{k+1},s_{k+2}),
%```

%where $s_k=t_{k-1}$.
%Specifically, this function is
%piecewise linear, nonzero only in two adjacent intervals
%$[s_{k},s_{k+1}]$ and $[s_{k+1},s_{k+2}]$, and continuous everywhere, with jumps
%in the first derivative at $s_{k}$, $s_{k+1}$, and $s_{k+2}$. We now
%rename this function $B(x;s_{k},s_{k+1},s_{k+2})$, and in this
%context we will call $s_{k},s_{k+1},s_{k+2}$ **knots**.
%
%It turns out that we can combine two of these three-knot, piecewise linear functions to
%get a four-knot, piecewise quadratic function, and so on. In general,
%we define a **B-spline of order $\mathbf{r**$} via the formula
%
```{math}
%  :label: bspline
%  B(x;s_1,s_2,\ldots,s_r) = \frac{x-s_1}{s_{r-1}-s_1} B(x;s_1,\ldots,s_{r-1}) +
%  \frac{s_{r}-x}{s_{r}-s_{2}} B(x;s_2,\ldots,s_r).
%```

%We assume that the knots satisfy $s_1<s_2<\cdots < s_r$. (However, in
%{ref}`prob-bsplines-coincidentknots` you are shown how to cope with
%knots that coincide, i.e., have multiplicity greater than one.)
%
%The formula {eq}`bspline` defines a B-spline on $r$ knots in terms
%of two B-splines on consecutive sets of $r-1$ knots. Each of these is in turn defined
%in terms of two B-splines on $r-2$ knots, but two of these shorter
%B-splines are the same function, so there are only three unique B-splines
%of $r-2$ knots:
%$B(x;s_1,\ldots,s_{r-2})$, $B(x;s_2,\ldots,s_{r-1})$, and $B(x;s_3,\ldots,s_r)$.
%The process continues recursively until one reaches $r-1$ unique B-splines on
%pairs of consecutive knots and can plug in {eq}`bspline1` for them.
%
%
````{prf:example}
%  \label{exa:Bsplinehand}
%  Let's find an explicit formula for $B(x;0,1,2,4)$ for $1 < x < 2$.
%  It's easiest to start at the bottom level, from the indicator
%  functions, and work upward. For $x$ in the given interval, the
%  indicators $B(x;0,1)$ and $B(x;2,4)$ are zero, and $B(x;1,2)\equiv
%  1$. There are two consecutive triples of knots:
%  
```{math}
\begin{split}
%    B(x;0,1,2) &= \frac{x-0}{1-0} \cdot 0 + \frac{2-x}{2-1}\cdot 1 =
%    (2-x),\\
%    B(x;1,2,4) &= \frac{x-1}{2-1} \cdot 1 + \frac{4-x}{4-2}\cdot 0 =
%    (x-1).
%  \end{split}
```

%  Finally, we arrive at
%  
```{math}
%    B(x;,0,1,2,4) = \frac{x-0}{2-0} \cdot (2-x)
%    + \frac{4-x}{4-1}\cdot(x-1) = \frac{4}{5} - \frac{(5x-8)^2}{30}.
%  ```

%````

%
%\begin{function}
%  (function-Bspline)=

````{proof:function} Bspline
****

```{code-block} julia
:lineno-start: 1

```
````
%  \caption{Evaluation of a B-spline.}
%\end{function}
%
%In {ref}`function-Bspline` we show an algorithm for evaluating a
%B-spline using {eq}`bspline1` and {eq}`bspline`
%directly. This is our first example of a recursive function—one that
%calls itself. The first line of the function is like a contract: given
%inputs that satisfy certain conditions, it returns an output
%satisfying stated conditions. In this case, given value(s) of $x$ and
%a vector of knots, {ref}`function-Bspline` returns the value of the
%B-spline. So in lines~16 and~17, where the recursive calls occur, we
%may conclude that "B1" and "B2" contain the values of
%B-splines on shorter sequences of knots. If we originally call the
%function on a vector of 4 knots, it makes two calls on vectors of 3
%knots, each of which makes two calls using vectors of 2 knots, which
%execute line~13. Only then can the calls with 3 knots be completed,
%and then finally the original call with 4 knots is completed. The
%recursive form of implementation is very simple, but (as is often the
%case) not efficient; in {ref}`prob-bsplines-implementation` you are
%asked to create a better one.
%
%One final note about line~13: it uses the fact
%that the logical operators, such as ">=" and "\&",
%return values that can be converted to numerical 1 (true) and 0
%(false) values using "double".
%
%
````{prf:example}
%  ````{admonition}
{doc}`demos/bsplinegraph`
````
%````

%
%## B-spline properties
%
%Some of the important properties of B-splines are summarized in the
%following theorem, which we give without proof.
%
%
````{proof:theorem}
%  Given a sequence of $r$ distinct knots $s_1,\ldots,s_r$ and the
%  definitions {eq}`bspline1` and {eq}`bspline`:
%  (theorem-bspline)=
%  \begin{enumerate}
%  \item $B(x;s_1,\ldots,s_r) \ge 0$ for all $x$.
%  \item On each interval $(s_k,s_{k+1})$, $B$ is a polynomial of
%    degree $r-2$ (as a function of $x$).
%  \item $B$ is nonzero only for $s_1\le x < s_{r}$. We equivalently say that $B$
%  is **supported** in the interval $[s_1,s_r)$.
%  \item For $r\ge 3$, $B$ has $r-3$ continuous derivatives in $x$. The
%     derivative of order $r-2$ has discontinuities only at the knots
%    $s_1,\ldots,s_{r}$.
%  %\item $B_{k,r}$ and $B_{j,r}$ are linearly independent if $j\neq k$.
%  \end{enumerate}
%````

%
%As a consequence of \thmref{bspline}, every linear combination of B-splines is a piecewise polynomial of degree $r-2$ and smoothness $C^{r-3}$ (i.e., $r-3$ continuous derivatives). We are next going to use them to create interpolants. But first notice that the B-splines of order greater than one are *not* cardinal functions, so we will lose one of the attractive features of hat functions.
%
%
%## Cubic spline interpolants
%
%We now return to the interpolation problem. Given interpolation nodes $t_0,\ldots,t_n$, we want to construct, say, a piecewise cubic interpolant with two continuous derivatives. Our job is to find a basis for the interpolation;, that is, a set of functions $\phi_0,\ldots,\phi_n$ such that the interpolant can be expressed as
%
```{math}
%  :label: interpbasis
%  p(x) = \sum_{k=1}^n c_k \phi_k(x),
%```

%for coefficients $c_k$ to be determined by $n+1$ interpolation conditions at the nodes.
%
%If we let the knots be the same as the nodes, then there are only $n-3$ consecutive 5-tuples to create candidate basis functions for {eq}`interpbasis`. This leaves us four short of the number needed to give us enough coefficients to match the interpolation conditions. We can borrow the trick of fictitious nodes from piecewise linear interpolation. We need four new knots in total, so a candidate knot set is $t_{-2},t_{-1},t_0,\ldots,t_n,t_{n+1},t_{n+2}$, where the first and last two nodes are located arbitrarily on their respective sides of the interval $[t_0,t_n]$.
%
%For a subtle reason to be considered in {ref}`prob-bsplines-notaknot`, however, this collection of knots will not produce what we want. Instead, we need the knot sequence
%\begin{gather*}
%  s_1=t_{-3},\,s_2=t_{-2},\,s_3=t_{-1},\,s_4=t_0,\\
%  s_5=t_2,\, s_6=t_3,\, \ldots,\, s_n=t_{n-3},\ s_{n+1}=t_{n-2},\\
%  s_{n+2}=t_n,\, s_{n+3}=t_{n+1},\, s_{n+4}=t_{n+2},\, s_{n+5}=t_{n+3}.
%\end{gather*}
%We left out nodes $t_1$ and $t_{n-1}$ and have *three* fictitious
%nodes on each side. Because $t_1$ and $t_{n-1}$ are not knots, there
%is no discontinuity in $p'''$ there. This arrangement is called a
%**not-a-knot spline**.
%
%There are $n+5$ not-a-knot knots, and thus $n+1$ consecutive
%5-tuples. These define third-order B-splines as the basis functions
%$\phi_k$ in {eq}`interpbasis`. What remains is to solve for the
%interpolation coefficients $c_k$ in that formula. By evaluation of $p$
%at each of the interpolation nodes (*not* the not-a-knot knots),
%we get a linear system of equations:
%
```{math}
%  :label: interpbasissys
%  \begin{bmatrix}
%    \phi_{1}(t_0) & \phi_{2}(t_0) & \cdots & \phi_{n+1}(t_0) \\
%    \phi_{1}(t_1) & \phi_{2}(t_1) & \cdots & \phi_{n+1}(t_1) \\
%    \vdots & \vdots & & \vdots\\
%    \phi_{1}(t_n) & \phi_{2}(t_n) & \cdots & \phi_{n+1}(t_n)
%  \end{bmatrix}
%  \begin{bmatrix}
%    c_1 \\ c_2 \\ \vdots \\ c_{n+1}
%  \end{bmatrix}
%  =
%  \begin{bmatrix}
%    y_0 \\ y_1 \\ \vdots \\ y_{n}
%  \end{bmatrix}.
%```

%The basis functions have limited support intervals. This makes many of
%the matrix entries zero. Specifically,
%
```{math}
%  B(s_j;s_i,\ldots,s_{i+4}) = 0, \qquad j<i \text{ or } j>i+4.
%```

%Consequently, the matrix in {eq}`interpbasissys` has only five
%nonzero entries per row, and upper bandwidth and lower bandwidth both
%equal to two. Because the matrix has constant bandwidth 5 for any value
%of $n$, only $O(n)$ operations are needed to find the coefficients.
%{ref}`function-spinterp` demonstrates a basic implementation.
%\begin{function}
%  (function-spinterp)=

````{proof:function} spinterp
****

```{code-block} julia
:lineno-start: 1
"""
spinterp(t,y)

Create a cubic not-a-knot spline interpolating function for data
values in `y` given at nodes in `t`.
"""
function spinterp(t,y)

    n = length(t)-1
    h = diff(t)         # differences of all adjacent pairs

    # Preliminary definitions.
    Z = zeros(n,n);
    In = I(n);  E = In[1:n-1,:];
    J = diagm(0=>ones(n),1=>-ones(n-1))
    H = diagm(0=>h)

    # Left endpoint interpolation:
    AL = [ In Z Z Z ]
    vL = y[1:n]

    # Right endpoint interpolation:
    AR = [ In H H^2 H^3 ];
    vR = y[2:n+1]

    # Continuity of first derivative:
    A1 = E*[ Z J 2*H 3*H^2 ]
    v1 = zeros(n-1)

    # Continuity of second derivative:
    A2 = E*[ Z Z J 3*H ]
    v2 = zeros(n-1)

    # Not-a-knot conditions:
    nakL = [ zeros(1,3*n) [1 -1 zeros(1,n-2)] ]
    nakR = [ zeros(1,3*n) [zeros(1,n-2) 1 -1] ]

    # Assemble and solve the full system.
    A = [ AL; AR; A1; A2; nakL; nakR ]
    v = [ vL; vR; v1; v2; 0; 0 ]
    z = A\v

    # Break the coefficients into separate vectors.
    rows = 1:n
    a = z[rows]
    b = z[n.+rows];  c = z[2*n.+rows];  d = z[3*n.+rows]
    S = [ Polynomial([a[k],b[k],c[k],d[k]]) for k = 1:n ]
    # This function evaluates the spline when called with a value
    # for x.
    function evaluate(x)
        k = findfirst(@. x<t)   # one greater than interval x belongs to
        k==1 && return NaN
        if isnothing(k)
            return x==t[end] ? y[end] : NaN
        end
        return S[k-1](x-t[k-1])
    end
    return evaluate
end
```
````
%  \caption{Cubic spline interpolation by B-splines.}
%\end{function}
%
%
````{prf:example}
%  ````{admonition}
{doc}`demos/bsplineinterp`
````
%````

%
%If we sample a function that has a continuous fourth derivative on the
%interval $[t_0,t_n]$, and $h$ is the maximum separation between
%adjacent nodes, then the convergence of the not-a-knot spline is
%fourth order, i.\ e., $\| p-f \|_\infty=O(h^4)$.  (cite deBoor?) This is
%superior to the second order convergence of piecewise linear
%interpolants. The news about condition numbers, however, is not so
%good. The condition number of not-a-knot interpolation can be at least as large
%as $2^n$, which grows very rapidly. In practice it's necessary to
%keep $n$ from being large, or introduce some constraints on the
%nodes.
%
%
%\begin{exercises}
%	\input{localfuncapprox/exercises/BSplines}
%\end{exercises}
%
%
%# Parameterization of curves
%\label{sec:parameterization}
%
%So far we have treated the data $(t_k,y_k)$ supposing that $y$ is an explicit
%function of an independent variable $x$. However, if we think of
%$(x_k,y_k)$ as points in the plane, we know that many plane curves—for example,
%circles—cannot be represented as graphs of any single-valued
%function $y=f(x)$. It's more general to introduce an independent parameter
%$t$ and treat both $x$ and $y$ as functions of $t$. In such cases, two independent interpolations with respect to $t$ can be used to represent an interpolating path.
%
%
%
````{prf:example}
%  Here we approximate an S-shaped curve using splines for both
%  coordinates.
%
%      \begin{minipage}[t]{2.5in}
%\begin{verbatim}
%% Nodal coordinate and parameter values.
%x = [63 41 27 40 59 72 62 43 27];
%y = [89 88 68 56 50 33 15 14 26];
%t = 1:9;
%
%% Interpolated values.
%ti = linspace(1,9,300);
%xi = spinterp(ti,t,x,);
%yi = spinterp(ti,t,y);
%
%plot(x,y,'o',xi,yi,'-')
%axis equal
%\end{verbatim}
%    \end{minipage}
%    \hfill \parbox[t]{1.75in}{
%      \psfrag{x}[t][t]{$x$}
%      \psfrag{y}[b][b]{$y$}
%      \raisebox{-2.1in}{\includegraphics[width=1.75in]{demoScurve}}
%    }
%````

%
%
````{prf:example}
%  \label{exa:dumbbell}
%  Here we used the mouse to click on 20 points to make a vague
%  dumbbell shape. Since we want to draw a closed curve, we treat the
%  data as periodic and use trigonometric interpolation.
%
%      \begin{minipage}[t]{2.5in}
%\begin{verbatim}
%% Nodal coordinate and parameter values.
%[x,y] = ginput(20);
%t = 0:19;
%
%% Interpolated values.
%ti = linspace(0,20,200);
%xi = triginterp(ti,t,x);
%yi = triginterp(ti,t,y);
%
%plot(x,y,'o',xi,yi,'-')
%axis equal
%\end{verbatim}
%    \end{minipage}
%    \hfill \parbox[t]{1.75in}{
%      \psfrag{x}[t][t]{$x$}
%      \psfrag{y}[b][b]{$y$}
%      \raisebox{-2.1in}{\includegraphics[width=1.75in]{dumbbell}}
%    }
%````

%
%In these examples the results were satisfactory. However, there are infinitely many ways to
%parameterize a curve.
%(Think of an ant traveling along the same curve with different nonuniform speed profiles.) The effects of a parameterization choice on the resulting curve
%are not easily predicted or expressed. Complicating matters is the fact that the
%interpolation methods in this chapter (except for the aesthetically
%challenged piecewise linear interpolants) have results that depend
%globally on each point in the data. In Example~\ref{exa:dumbbell} we
%can see some small wiggles in the curve due to this
%effect. Trigonometric interpolation is especially problematic, since a
%large, sudden change in a coordinate value will cause a Gibbs
%overshoot effect.
%
%A popular alternative allows more intuitive control
%over the shape of the interpolation curve. We begin with piecewise
%cubic interpolation. The four parameters on each cubic subinterval are
%determined by specified points *and* slopes at the endpoints of
%the subinterval. Suppose two adjacent nodes have parameter values
%$t=0$ and $t=1$. We want to find two cubic functions, $x(t)$ and
%$y(t)$, where
%
```{math}
%  x(0),\quad x(1),\quad y(0),\quad y(1),\quad \frac{d y}{d x}(0),\quad \frac{d y}{d x}(1)
%```

%are given.
%
%Note that the two cubics have eight parameters between
%them, but only six conditions for the curve are specified. The root of
%the problem is that
%
```{math}
%  \frac{d y}{d x} = \frac{dy/dt}{dx/dt}
%```

%and so if $x'$ and $y'$ are scaled by the same constant, the curve's slope
%is unchanged.  It turns out that this "problem" is actually a
%strength. The scaling factors in $x'$ and $y'$ affect how closely the
%interpolant follows the tangent line at the node. Graphically, this is
%communicated via the use of **control points**, as illustrated in
%Figure~\ref{fig:bezier}. The direction from a node to its control
%point gives the slope at the end, while the distance from the node to
%the control point sets the scaling factor in the derivatives of $x(t)$
%and $y(t)$ and influences how strongly the tangency persists.
%
%\begin{figure}[tb]
%  \centering
%  \psfrag{r0}[t][t]{$\mathbf{r}_0$}
%  \psfrag{r3}[t][t]{$\mathbf{r}_3$}
%  \psfrag{r1}[b][b]{$\mathbf{r}_1$}
%  \psfrag{r2}[b][b]{$\mathbf{r}_2$}
%  \includegraphics[width=\textwidth]{bezier}
%  \caption{Cubic B\'ezier curves. In both cases the end values and
%    slopes are the same. When the distance from node to control point is
%    increased, the relative weight of the tangency condition is
%    increased.}
%  \label{fig:bezier}
%\end{figure}
%
%At this point it becomes more convenient to switch to a vector
%notation. Let $\mathbf{r}(t)$ be the vector $[x(t),y(t)]^T$, defined for
%$0\le t \le 1$. The end nodes will be denoted $\mathbf{r}_0$ and $\mathbf{r}_3$,
%and their control points are be $\mathbf{r}_1$ and $\mathbf{r}_2$, respectively.
%We want the curve to satisfy
%
```{math}
%  \mathbf{r}(0) = \mathbf{r}_0, \qquad \mathbf{r}(1) = \mathbf{r}_3, \qquad \mathbf{r}'(0) =
%  3(\mathbf{r}_1-\mathbf{r}_0), \qquad \mathbf{r}'(1) = 3(\mathbf{r}_3-\mathbf{r}_2).
%```

%(The seemingly extraneous factors of $3$ in $\mathbf{r}'$ make the numbers
%come out cleaner in general formulas.) It is straightforward to show that
%
```{math}
%  :label: bezier
%  \mathbf{r}(t) = (1-t)^3 \mathbf{r}_0 + 3t(1-t)^2\mathbf{r}_1
%  + 3t^2(1-t) \mathbf{r}_2 + t^3 \mathbf{r}_3
%```

%is the unique cubic function satisfying these conditions.
%Equation {eq}`bezier` describes a \textbf{cubic B\'ezier
%curve}.
%
%
%In \matlab\ we might implement a simple cubic B\'ezier curve as shown
%in Function~\ref{fun:bezier}.
%\begin{function}
%  \Mfile{interpolation/Bezier.m}
%  \caption{Cubic B\'ezier curve evaluation.}
%  \label{fun:bezier}
%\end{function}
%The two columns of "node" hold the endpoints and the two
%columns of "control" are their associated control points. The
%vector "t" should have parameter values between 0 and 1. The
%last line is a little subtle, making use of vector outer products to
%create a matrix of results in which each column is a point on the
%curve for one value in "t". One convenient aspect of
%"bezier" is that it can be used exactly as written for a
%B\'ezier curve in three or more dimensions, as determined implicitly
%by the input.
%
%In practice, one could have many nodes and evaluate a different
%B\'ezier curve piecewise between each pair of neighboring nodes.
%Interior nodes would have a control point on each side. Note that each
%piece of the curve is completely determined by local data; conversely,
%moving a node affects only the segments of the curve passing through
%it.  B\'ezier interpolation is often implemented interactively, with
%the computer guessing the initial control points and the user moving
%them and the nodes to get the desired shape.
%
%
````{prf:example}
%\label{exa:bezier}
%In this example, we consider making a stylized ``V" using three
%B\'ezier curves put together.  Certainly it is easiest to
%graphically input data for nodes and control points and subsequently
%manipulate them. For this example we only need 4 nodes (and control
%points), and in the absence of graphical input, one could draw the
%desired letter on a grid to estimate the node and control point
%locations.  The following code snippet shows the nodes and control
%points that draw the letter.  Note that Function~\ref{fun:bezier} is
%called 3 times and the results are concatenated prior to plotting;
%the control points are plotted as squares. The result is shown in Figure~\ref{fig:exabezier}.
%\begin{verbatim}
%% each column is a node or control point
%nodes = [-0.3,1.3;0.1,2.5;-0.3,-1.8;0.8,2.6]';
%control=[-0.1,1.7;-0.2,1.5;-0.1,-0.8;0.5,2.4]';
%t = linspace(0,1,21);       % 21 points per curve
%letter = [];                % initialize data to plot
%curves = length(nodes)-1;   % curves are between nodes
%for k=1:(length(nodes)-1)   % one curve for each interval
%  letter = [letter,Bezier(nodes(:,k:k+1),control(:,k:k+1),t)];
%end
%% plot interpolant, control points and nodes
%plot(letter(1,:),letter(2,:),control(1,:),...
%  control(2,:),'s',nodes(1,:),nodes(2,:),'o');
%\end{verbatim}
%````

%\begin{figure}
%  \centering
%  \includegraphics[width=\textwidth]{letterexambezier}
%  \caption{B\'ezier cubic curve interpolation (see {doc}`demos/bezier`).}
%  \label{fig:exabezier}
%\end{figure}
%
%The four polynomials appearing in the definition {eq}`bezier` are
%called **Bernstein polynomials** and can be generalized to any
%degree. Furthermore, the Bernstein polynomials themselves are special
%cases of B-splines that have repeated knots (see
%Problem~\ref{pro:coincidentknots} on
%p.~\pageref{pro:coincidentknots}). In fact B-splines are usually
%preferred to B\'{e}zier curves in advanced implementations.
%
%## Problems
%\begin{exercises}
%\item Manipulate the example ``V" to see how the curves change; provide two alternative
%styles of V's.
%by moving nodes and control points.  Can you improve the appearance of this
%letter?
%\item Create a ``Z" in a similar style by modifying the example.
%\item Improve the example function to take graphical input for the nodes first
%and then the control points to generate the initial attempt to make a letter (say).
%\end{exercises}

