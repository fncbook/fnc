---
numbering: false
title: Python setup
---
(section-setup-python)=
# Setting up Python for this book

Python, and all the packages this book depends on, is free and open-source. Much of the functionality outside the core is distributed via packages that need to be installed once per system.

## Installing Python

There are tons of ways and guides. [Anaconda](https://www.anaconda.com/download) is a popular distribution that comes with many packages pre-installed. The free version is fine and, while it is overkill for this book, it is a solid choice for beginners. Command-line users are encouraged to consider [Mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html), which is a fast version of Anaconda. If you are a VS Code user, you can install Python from within the editor.

## Installing the book's functions

The book relies on functions that are distributed as a [PyPi package](https://pypi.org/project/fncbook/). You can install it by typing `pip install fncbook` at a command prompt, or by following instructions for your installation method.

## Using packages

Python uses `import` to load packages. For all the demos in the book, its package is loaded via

``` python
import fncbook as FNC
```

(`FNC` is chosen to be uniform with the Julia versions of the codes. You can choose any name you like.) Then you access its functions like `FNC.horner`, etc. You can also type `FNC.` and then press {kbd}`Tab` to see a list of available functions.

## Packages used in the book

In order to avoid repeating low-information code, the book demos are run with a few packages installed and always imported:

``` python
from numpy import *
from matplotlib.pyplot import *
from numpy.linalg import norm
from prettytable import PrettyTable
from timeit import default_timer as timer
```

Note that a lot of sources discourage using `import *` because it can lead to so-called namespace pollution. But we rely so heavily on `numpy` and `pyplot` that it is convenient here. 

Other packages are loaded as needed. These include

- `scipy`
- `scikit-image`
- `networkx`
- `rogues`

Only `scipy` is essential. The others are used for a few applications and illustrations.

## Coding environments

You *could* interact with Python only by typing in at the prompt (also called the REPL) and then pasting the results into a word processor, but you can do much, much better. The most popular ways to use Julia are:

- **[Jupyter lab](https://jupyter.org)**. This is a notebook-based interface that mixes cells having text and code, including text and graphical output. This entire book is based on the notebook architecture. You write and run code within your web browser, but the files are local.
- **[VS Code](https://code.visualstudio.com)**. This is a full-featured code editor that can be extended with lots of [Python-specific tools](https://code.visualstudio.com/docs/languages/python). It can also write and run Jupyter notebooks. This book (the version you are reading now, anyway) was written in VS Code.