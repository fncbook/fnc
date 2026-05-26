---
numbering:
  enumerator: 7.1.%s
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---
```{code-cell}
:tags: [remove-cell]
from numpy import *
from scipy import linalg
from scipy.linalg import norm
from matplotlib.pyplot import *
from prettytable import PrettyTable
from timeit import default_timer as timer
import sys
sys.path.append('fncbook/')
import fncbook as FNC

# This (optional) block is for improving the display of plots.
# from IPython.display import set_matplotlib_formats
# set_matplotlib_formats("svg","pdf")
# %config InlineBackend.figure_format = 'svg'
rcParams["figure.figsize"] = [7, 4]
rcParams["lines.linewidth"] = 2
rcParams["lines.markersize"] = 4
rcParams['animation.html'] = "jshtml"  # or try "html5"
```

(section-matrixanaly-insight)=

# From matrix to insight

Any two-dimensional array of numbers may be interpreted as a matrix. Whether or not this is the only point of view that matters to a particular application, it does lead to certain types of analysis. The related mathematical and computational tools are universally applicable and find diverse uses.

## Tables as matrices

Tables are used to represent variation of a quantity with respect to two variables. These variables may be encoded as the rows and columns of a matrix.

````{prf:example}
A *corpus* is a collection of text documents. A *term-document matrix* has one column for each document and one row for each unique term appearing in the corpus. The $(i,j)$ entry of the matrix is the number of times term $i$ appears in document $j$. That is, column $j$ of the matrix is a term-frequency vector quantifying all occurrences of the indexed terms. A new document could be represented by its term-frequency vector, which is then comparable to the columns of the matrix. Or, a new term could be represented by counting its appearances in all of the documents and be compared to the rows of the matrix.

It turns out that by finding the {ref}`singular value decomposition <section-matrixanaly-svd>` of the term-document matrix, the strongest patterns within the corpus can be isolated, frequently corresponding to what we interpret as textual meaning. This is known as *latent semantic analysis.*
````

::::{prf:example}
Each vote cast in the U. S. Congress is [available for download](https://www.congress.gov/roll-call-votes). We can put members of Congress along the columns of a matrix and bills along the rows, recording a number that codes for "yea,"  "nay," "none," etc. The {ref}`singular value decomposition <section-matrixanaly-svd>` can reveal an objective, reproducible analysis of the partisanship and cooperation of individual members.
::::

````{prf:example}
In 2006 the online video service Netflix started an open competition for a $1 million prize. They provided a data set of 100,480,507 ratings (one to five stars) made by 480,189 users for 17,770 movies. Each rating is implicitly an entry in a 17,770-by-480,189 matrix. The object of the prize was to predict a user's ratings for movies they had not rated. This is known as a *matrix completion problem.* (It took 6 days for a contestant to improve on Netflix's private algorithm, and in 2009 the million-dollar prize was awarded to a team that had improved the performance by over 10%.)
````

## Graphs as matrices

```{index} ! graph, ! adjacency matrix
```

::::{prf:definition} Graphs and adjacency matrices
:label: definition-adjacency-matrix
A {term}`graph` or *network* consists of a set $V$ of **nodes** and a set $E$ of **edges**, each of which is an ordered pair of nodes. If there is an edge $(v_i,v_j)$, then we say that node $i$ is **adjacent** to node $j$. The graph is **undirected** if for every edge $(v_i,v_j)$, the pair $(v_j,v_i)$ is also an edge; otherwise the graph is **directed**.

The {term}`adjacency matrix` of a graph with $n$ nodes $V$ and edge set $E$ is the $n\times n$ matrix whose elements are

```{math}
:label: adjmat
A_{ij} =
\begin{cases}
1 & \text{if $(v_i,v_j)\in E$ (i.e., node $i$ is adjacent to node $j$)},\\
0 & \text{otherwise}.
\end{cases}
```

::::

:::{note}
In an undirected graph, the edges $(v_i,v_j)$ and $(v_j,v_i)$ are equivalent and may be identified as a single edge, depending on the context.
:::

Graphs are a useful way to represent the link structure of social networks, airline routes, power grids, sports teams, and web pages, to name a few examples. The natural interpretation is that the edge $(v_i,v_j)$ denotes a link from node $i$ to node $j$, in which case we say that node $i$ is **adjacent** to node $j$. One usually visualizes small graphs by drawing points for nodes and arrows or lines for the edges.

Here are some elementary results about adjacency matrices.

::::{prf:theorem}
:label: theorem-insight-adjmat
For any {term}`graph` with {term}`adjacency matrix` $\mathbf{A}$,

1. The graph is undirected if and only if $\mathbf{A}$ is symmetric, and
2. For any positive integer $k$, the $(i,j)$ element of $\mathbf{A}^k$ is the number of ways to walk from node $i$ to node $j$ by following along exactly $k$ edges.
::::

::::{prf:proof}
:enumerated: false

Part 1 follows immediately from the definitions. Part 2 is clearly true for $k=1$. Assume inductively that it is true for $k-1$. Each walk of length $k$ from node $i$ to node $j$ must be a walk of length $k-1$ from $i$ to node $p$, then a walk of length 1 from node $p$ to node $j$. The total number of such walks is therefore

$$
\sum_{p=1}^n [\mathbf{A}^{k-1}]_{ip} \cdot A_{pj},
$$

which is the $(i,j)$ element of $\mathbf{A}^k$.
::::

::::{prf:example} Adjacency matrix
:label: demo-insight-graph

Here we create an adjacency matrix for a graph on four nodes.

```{code-cell}
A = array([[0, 1, 0, 0], [1, 0, 0, 0], [1, 1, 0, 1], [0, 1, 1, 0]])
print(A)
```

The `networkx` package has many functions for working with graphs. Here, we instruct it to create a directed graph from the adjacency matrix, then make a drawing of it.

```{index} ! Python; networkx
```

```{code-cell}
import networkx as nx
G = nx.from_numpy_array(A, create_using=nx.DiGraph)
nx.draw(G, with_labels=True, node_color="yellow")
```

Here are the counts of all walks of length 3 in the graph:

```{code-cell}
from numpy.linalg import matrix_power
print(matrix_power(A, 3))
```

If the adjacency matrix is symmetric, the result is an undirected graph: all edges connect in both directions.

```{code-cell}
A = array([[0, 1, 1, 0], [1, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]])
G = nx.from_numpy_array(A, create_using=nx.Graph)
nx.draw(G, with_labels=True, node_color="yellow")
```

::::

The representation of a graph by its adjacency matrix opens up the possibility of many kinds of analysis of the graph. One might ask whether the nodes admit a natural partition into clusters, for example. Or one might ask to rank the nodes in order of importance to the network as determined by some objective criteria—an application made famous by Google's PageRank algorithm, and one which is mathematically stated as an {ref}`eigenvalue <section-matrixanaly-evd>` problem.

## Images as matrices

```{index} image (as a matrix)
```

Computers typically represent images as rectangular arrays of pixels, each of which is colored according to numerical values for red (R), green (G), and blue (B) components of white light. Most often, these are given as integers in the range from zero (no color) to 255 (full color). Thus, an image that is $m$-by-$n$ pixels can be stored as an $m$-by-$n$-by-3 array of integer values. In Julia, we can work with an $m\times n$ matrix of 3-vectors representing entire colors.

::::{prf:example} Images as matrices
:label: demo-insight-image

```{index} ! Python; Images
```

We will use a test image from the well-known `scikit-image` package.

```{code-cell}
from skimage import data as testimages
img = getattr(testimages, "coffee")()
imshow(img)
```

The variable `img` is a matrix.

```{code-cell}
size(img)
```

However, its entries are colors, not numbers.

```{code-cell}
print(f"image has shape {img.shape}")
print(f"first pixel has value {img[0, 0]}")
```

The three values at each pixel are for intensities of red, green, and blue. We can convert each of those layers into an ordinary matrix of values between 0 and 255, which is maximum intensity.

```{code-cell}
R = img[:, :, 0]
print("upper left corner of the red plane is:")
print(R[:5, :5])
print(f"red channel values range from {R.min()} to {R.max()}")
```

It may also be convenient to convert the image to grayscale, which has just one layer of values from zero (black) to one (white).

```{code-cell}
from skimage.color import rgb2gray
A = rgb2gray(img)
A[:5, :5]
print("upper left corner of grayscale:")
print(A[:5, :5])
print(f"gray values range from {A.min()} to {A.max()}")
```

```{code-cell}
imshow(A, cmap='gray')
axis('off');
```

Some changes we make to the grayscale matrix are easy to interpret visually.

```{code-cell}
imshow(flipud(A), cmap='gray')
axis('off');
```

::::

Representation of an image as a matrix allows us to describe some common image operations in terms of linear algebra. For example, in {numref}`section-matrixanaly-svd` we will use the singular value decomposition to compress the information, and in {numref}`section-krylov-matrixfree` we will see how to apply and remove blurring effects.

## Exercises

``````{exercise}
:label: problem-insight-terms
✍ Consider the terms *numerical*, *analysis*, and *fun*. Write out the term-document matrix for the following corpus of "documents."

1. Numerical analysis is the most fun type of analysis.

2. It's fun to produce numerical values for the digits of pi.

3. Complex analysis is a beautiful branch of mathematics.
``````

`````{exercise}
:label: problem-insight-adjacency1
✍ Write out the adjacency matrix for the following graph on six nodes.

```{image} littlegraph.png
:alt: little graph
:width: 500px
:align: center
```
``````

``````{exercise}
:label: problem-insight-adjacency2
✍ Here is a graph adjacency matrix.

```{math}
:numbered: false
\begin{bmatrix}
0 & 1 & 0 & 1 & 0 & 1 & 0 \\
1 & 0 & 1 & 0 & 0 & 1 & 0 \\
0 & 1 & 0 & 0 & 1 & 0 & 1 \\
1 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 1 & 0 & 0 & 1 & 0 \\
1 & 1 & 0 & 0 & 1 & 0 & 0 \\
0 & 0 & 1 & 0 & 0 & 0 & 0   
\end{bmatrix}
```

**(a)** How many vertices are adjacent to vertex 5? (Assume node numbering starts at 1.)

**(b)** Is the graph directed or undirected?

**(c)** How many edges are in the graph?

**(d)** Draw the graph. 
``````

``````{exercise}
:label: problem-insight-image
⌨ Refer to @demo-insight-image on loading and displaying images. Choose a test image of your liking.

**(a)** Display the test image upside-down.

**(b)** Display it mirror-reversed from left to right. 

**(c)** Display the image so that it is cropped to isolate a part of the subject. 

``````

``````{exercise}
:label: problem-insight-actors
⌨ For this problem you need to download and import data:
Download [actors.mtx](actors.mtx) by clicking the link and saving (you may need to fix the file name).
``` python
import scipy.io as spio
A = spio.mmread("actors.mtx")
```

Based on data provided by the Self-Organized Networks Database at the University of Notre Dame, the matrix `A` contains information about the appearances of 392,400 actors in 127,823 movies, as given by the [Internet Movie Database](wiki:IMDb). It has $A_{ij}=1$ if actor $j$ appeared in movie $i$ and zero elements elsewhere.

**(a)** What is the maximum number of actors appearing in any one movie?

**(b)** How many actors appeared in exactly three movies?

**(c)** Define $\mathbf{C}=\mathbf{A}^T\mathbf{A}$. How many nonzero entries does $\mathbf{C}$ have? What is the interpretation of $C_{ij}$?
``````
