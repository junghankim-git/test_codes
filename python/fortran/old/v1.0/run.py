#!/app/compilers/gcc/4.7.1/applib1/Python/2.7.3/bin/python
###!/opt/local/bin/python
import os
import sys
from mymodule import *

infilepath = './codes'
oufilepath = './result'

converter = F90Converter(infilepath, oufilepath)

converter.Driver()


