#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib import cm
import os
import sys
import random
SHARE_DIR=os.environ['PYTHON_SHARE_DIR']
sys.path.append(SHARE_DIR)
from class_beta_dist import *



beta = beta_dist(50,2,2)
beta.draw_pdf()
