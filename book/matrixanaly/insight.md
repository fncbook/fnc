# From matrix to insight

Any two-dimensional array of numbers may be interpreted as a matrix. Whether or not this is the only point of view that matters to a particular application, it does lead to particular types of analysis. The related mathematical and computational tools are universally applicable and find diverse uses.

## Tables as matrices

Tables are used to represent variation of a quantity with respect to two variables. These variables may be encoded as the rows and columns of a matrix.

````{prf:example}
  Suppose we have a *corpus*, or collection of text documents. A *term-document matrix* has one column for each document and one row for each unique term appearing in the corpus. The $(i,j)$ entry of the matrix is the number of times term $i$ appears in document $j$. That is, column $j$ of the matrix is a term-frequency vector quantifying all occurrences of the indexed terms. A new document could be represented by its term-frequency vector, which is then comparable to the columns of the matrix. Or, a new term could be represented by counting its appearances in all of the documents and be compared to the rows of the matrix.

  It turns out that by finding the [singular value decomposition](svd) of the term-document matrix, the strongest patterns within the corpus can be isolated, frequently corresponding to what we interpret as textual meaning. This is known as *latent semantic analysis.*
````

````{prf:example}
  The website www.congress.gov/roll-call-votes offers data on all the votes cast in each session of the U. S. Congress. We can put members of Congress along the columns of a matrix and bills along the rows, recording a number that codes for "yea,"  "nay," "none," etc. The [singular value decomposition](svd) can reveal an objective, reproducible analysis of the partisanship and cooperation of individual members.
````



````{prf:example}
  In 2006 the online video service Netflix started an open competition for a $1 million prize. They provided a data set of 100,480,507 ratings (one to five stars) made by 480,189 users for 17,770 movies. Each rating is implicitly an entry in a 17,770-by-480,189 matrix. The object of the prize was to predict a user's ratings for movies they had not rated. This is known as a *matrix completion problem.* (It took 6 days for a contestant to improve on Netflix's private algorithm, and in 2009 the million-dollar prize was awarded to a team that had improved the performance by over 10%.)
````


## Graphs as matrices

```{index} graph nodes and edges
```
```{index} adjacency matrix
```
```{index} matrix; adjacency
```
An important concept in mathematics is that of a {term}`graph` or network. A graph consists of a set $V$ of **nodes** and a set $E$ of **edges**, each of which is an ordered pair of nodes. The natural interpretation is that the edge $(v_i,v_j)$ denotes a link from node $i$ to node $j$, in which case we say that node $i$ is adjacent to node $j$. One usually visualizes small graphs by drawing points for nodes and arrows or lines for the edges.

Graphs are useful because they are the simplest way to represent link structure—of social networks, airline routes, power grids, sports teams, and web pages, to name a few examples. They also have close ties to linear algebra. Prominent among these is the {term}`adjacency matrix` of the graph. If the graph has $n$ nodes, then its $n\times n$ adjacency matrix $\mathbf{A}$ has elements

```{math}
:label: adjmat
A_{ij} =
\begin{cases}
1,& \text{if $(v_i,v_j)\in E$ (i.e., there is an edge from node $i$ to node $j$)},\\
0,& \text{otherwise}.
\end{cases}
```

```{prf:example} Julia demo
:class: demo
:label: demos-insight-graph
{doc}`demos/insight-graph`
```

The representation of a graph by its adjacency matrix opens up the possibility for many kinds of analysis of the graph. One might ask whether the nodes admit a natural partition into clusters, for example. Or one might ask to rank the nodes in order of importance to the network as determined by some objective criteria—an application made famous by Google's PageRank algorithm, and one which is mathematically stated as an [eigenvalue problem](evd).

## Images as matrices

```{index} matrix; as image
```

Computers most often represent images as rectangular arrays of pixels, each of which is colored according to numerical values for red (R), green (G), and blue (B) components of white light. Typically these are given as integers in the range from zero (no color) to 255 (full color). Thus, an image that is $m$-by-$n$ pixels can be stored as an $m$-by-$n$-by-3 array of integer values.

We will simplify the representation by considering images represented using pixels that can take only shades of gray. We will also use floating-point numbers rather than integers, so that we can operate on them using real arithmetic, though we will stay with the convention that values should be in the range $[0,255]$. (Pixels below zero or above 255 will be colored pure black or pure white, respectively.)

```{prf:example} Julia demo
:class: demo
:label: demos-insight-image
{doc}`demos/insight-image`
```

Representation of an image as a matrix allows us to describe some common image operations in terms of linear algebra. Furthermore, the [singular value decomposition](svd) can be used to compress the information.

## Exercises

1. ✍ Consider the terms *numerical*, *analysis*, and *fun*. Write out the term-document matrix for the following statements:

    **(a)** Numerical analysis is the most fun type of analysis.

    **(b)** It's fun to produce numerical values for the digits of pi.

    **(c)** Complex analysis is a beautiful branch of mathematics.

2. ✍ Write out the adjacency matrix for the following graph on six nodes.
    
    ```{image} demos/littlegraph.png
    :alt: little graph
    :width: 500px
    :align: center
    ```

3. ✍ Here is a graph adjacency matrix.
  
    :::{math}
    \begin{bmatrix}
    0 & 1 & 0 & 1 & 0 & 1 & 0 \\
    1 & 0 & 1 & 0 & 0 & 1 & 0 \\
    0 & 1 & 0 & 0 & 1 & 0 & 1 \\
    1 & 0 & 0 & 0 & 0 & 0 & 0 \\
    0 & 0 & 1 & 0 & 0 & 1 & 0 \\
    1 & 1 & 0 & 0 & 1 & 0 & 0 \\
    0 & 0 & 1 & 0 & 0 & 0 & 0   
    \end{bmatrix}
    :::

    **(a)** How many vertices are adjacent to vertex 5?

    **(b)** How many edges are in the graph?

    **(c)** Draw the graph. 

4. ⌨ Refer to the [image manipulation demonstration](`demos/../demos/insight-image`).

    **(a)** Display the "peppers" test image upside-down.

    **(b)** Display it mirror-reversed from left to right. 

    **(c)** Display the image so that it is cropped to include only the garlic bulb in the lower right of the original picture. 
