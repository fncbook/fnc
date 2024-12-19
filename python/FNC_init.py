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