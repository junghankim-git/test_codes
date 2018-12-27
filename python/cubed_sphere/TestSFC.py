#!/usr/bin/env python
from numpy import *
from scipy import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import time
from matplotlib import cm
import os
import sys
# load My Class
sys.path.append('/home/jhkim/Study/Library/Shared/Python')
from LibSpaceFillingCurve import *
from LibControlMesh import *

# nsfc
# n=2:hi4bert, n=3:peano, n=5:cinco
nsfc   = 5
ncols  = nsfc
nrows  = nsfc



isfc = SFC(nsfc)


# Faces (all canvas domain)
sfc    = array([[0 for i in arange(nrows)] for j in arange(ncols)])
csfc   = array([[0 for i in arange(nrows)] for j in arange(ncols)])
cidx   = array([[[0,0] for i in arange(nrows)] for j in arange(ncols)])

isfc.GenSFC()
a, b = isfc.GetMajorJoiner(2, 0, 1)
isfc.GenPrimitiveCurve(nsfc, 0, 0, 0, 0, 1, sfc)


mymesh = Mesh()
mymesh.PrintMesh2D(2, b, '## in Main ##')
mymesh.PrintMesh2D(nsfc, sfc, '## in Main ##')
