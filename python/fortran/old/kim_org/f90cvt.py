#!/app/compilers/gcc/4.7.1/applib1/Python/2.7.3/bin/python
###!/opt/local/bin/python
import os
import sys
from optparse import OptionParser
SHARE_DIR=os.environ['PYTHON_SHARE_DIR']
sys.path.append(SHARE_DIR)
from f90convert import *

opt = OptionParser()

# action: 'store', 'store_const', 'append', 'count', 'callback'
opt.add_option('-s', '--source',      dest='src',  default='./src', action='store',  help='source')
opt.add_option('-d', '--destination', dest='dst',  default='./dst', action='store',  help='destination')

(options, args) = opt.parse_args()

converter = F90Converter(options.src, options.dst)

converter.driver()


