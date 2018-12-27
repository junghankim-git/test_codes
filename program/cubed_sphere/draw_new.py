#!/usr/bin/env python

# load My Class
import sys
SHARE_DIR='../../share/python'
sys.path.append(SHARE_DIR)
from class_cubed_sphere import *


infilename = './cs.nc'

cs = cubed_sphere(infilename)

# ne=15
cs.draw_paper_figure_1_1('figure_01_1')

# ne=4
#cs.draw_paper_figure_1_2('figure_01_2')

# ne=2, nprocs=5
cs.draw_paper_figure_2_1('figure_02_1')
#cs.draw_paper_figure_2_2('figure_02_2')
cs.draw_paper_figure_3_1('figure_03_1')
cs.draw_paper_figure_3_2('figure_03_2')

# ne=15, nprocs=16
cs.draw_paper_figure_3_2('figure_03_2')
#cs.draw_paper_figure_t1('figure_t1')

cs.show()
