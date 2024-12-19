from numpy import *
from numpy.linalg import norm
from matplotlib.pyplot import *
from prettytable import PrettyTable
from timeit import default_timer as timer
import sys
sys.path.append('pkg/')
import FNC
import importlib
importlib.reload(FNC)

# This (optional) block is for improving the display of plots.
# from IPython.display import set_matplotlib_formats
# set_matplotlib_formats("svg","pdf")
# %config InlineBackend.figure_format = 'svg'
rcParams["figure.figsize"] = [7, 4]
rcParams["lines.linewidth"] = 2
rcParams["lines.markersize"] = 4
rcParams['animation.html'] = "jshtml"  # or try "html5"