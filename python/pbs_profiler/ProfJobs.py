#!/usr/bin/env python
#-------------------------------------------------------------------------------
# MODULE KiapsGrid
#
#> @brief
#>  - show status of batch jobs
#>
#> @date 13JAN2015
#>  - JH KIM(jh.kim@kiaps.org) : First written 
#>    main module : /home/jhkim/Study/Library/Shared/Python/ProfJob/LibProfJob.py
#> @date 13FEB2015
#>  - JH KIM(jh.kim@kiaps.org) : Modified the main module
#> @date 17FEB2015
#>  - JH KIM(jh.kim@kiaps.org) : Added exclusive nodes and jobs (option: -j) informations
#-------------------------------------------------------------------------------

import os
import sys
from optparse import OptionParser
#SHARE_DIR=os.environ['PYTHON_SHARE_DIR']
SHARE_DIR='/home/jhkim/work/share/python'
sys.path.append(SHARE_DIR)
from class_pbs_jobs import *

opt = OptionParser()

# action: 'store', 'store_const', 'append', 'count', 'callback'
opt.add_option('-s', '--status',  dest='status',  default=True,  action='store_false', help='Not print status (default)')
opt.add_option('-q', '--qstat',   dest='qstat',   default=False, action='store_true',  help='Print qstat')
opt.add_option('-u', '--user',    dest='user',    default=False, action='store_true',  help='Print user infomation')
opt.add_option('-j', '--job',     dest='job',     default=False, action='store_true',  help='Print job infomation')

(options, args) = opt.parse_args()


pbs = pbs_jobs('Gaon2')
pbs.set_job_info()

if options.qstat:
   pbs.print_qstat()

if options.status:
   pbs.print_status(options.user)

if options.job:
   pbs.print_jobs()

if options.user:
   pbs.set_user_info()
   pbs.print_users()

if not options.user:
   pbs.print_jobs_cpus()

pbs.print_nodes_cpus()
pbs.print_available()

print '------------------------------------------------------------------------------'
print '# Contact: Junghan Kim (jh.kim@kiaps.org)'
print '------------------------------------------------------------------------------'
