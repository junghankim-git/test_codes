# -*- coding: utf-8 -*-
# for '║'
import os
import sys
import time
import datetime
import subprocess
import numpy as np
sys.path.append('/home/jhkim/work/share/python')
from class_parser import *

#######################################################################################
# Description: 현재 배치 잡들을 상태를 확인하고 가용한 노드와 cpu 정보를 출력한다.
#
# Date: 15JAN2015
#   - JH KIM: 처음 작성
#
# Date: 29JAN2015
#   - JH KIM: exclusive 노드 정보 추가
#
# Date: 05FEB2015
#   - JH KIM: serial job, offline node 정보 추가
#
# Date: 06FEB2015
#   - JH KIM: Gaon1에 대한 호환성 추가
#
# Date: 12FEB2015
#   - JH KIM: ending job에 대한 정보 추가
#
# Date: 17FEB2015
#   - JH KIM: Job 들에 대한 walltime 정보 추가
#
# Date: 18MAR2015
#   - JH KIM: 확장된 Gaon2에 적용
#
# Date: 27APR2015
#   - JH KIM: normal, normal20의 구분
#
# Date: 02JUN2015
#   - JH KIM: job에 queue 정보 추가
#
#######################################################################################

class Jobs:
   def __init__(self):
      self.id       = -1
      self.name     = ''
      self.queue    = ''
      self.user     = ''
      self.ncpus    = -1
      self.ncpus_n  = []
      self.nnodes   = -1
      self.nodes    = []
      self.walltime = [-1.0, -1.0, -1.0, -1.0, -1.0] # total wall time [sec], days, hours, minute, seconds
      self.status   = 'N'
      self.exclu    = False

   def CopyJob(self, tjob):
      self.id       = tjob.id
      self.name     = tjob.name
      self.queue    = tjob.queue
      self.user     = tjob.user
      self.ncpus    = tjob.ncpus
      self.ncpus_n  = tjob.ncpus_n
      self.nnodes   = tjob.nnodes
      self.nodes    = tjob.nodes
      self.walltime = tjob.walltime
      self.status   = tjob.status
      self.exclu    = tjob.exclu



class Users:
   def __init__(self):
      self.id       = -1
      self.name     = ''
      self.nTotJobs = 0
      self.nRunJobs = 0
      self.nHolJobs = 0
      self.nWaiJobs = 0
      self.nEndJobs = 0
      self.TotCPUs  = 0
      self.RunCPUs  = 0
      self.HolCPUs  = 0
      self.WaiCPUs  = 0
      self.EndCPUs  = 0

class pbs_jobs:
   def __init__(self,system='Gaon2'):
      # system info
      self.system  = system
      if system=='Gaon1':
         self.nNodes   = 5
         self.maxCPUs  = [64 for i in range(self.nNodes)]
      else:
         # 1-92: 16, 93-134:20
         self.mCPUs    = [16, 20]
         self.nNodes   = 135
         self.maxCPUs  = [0 for i in range(self.nNodes)]
         for i in range(self.nNodes): 
            if i < 92:
              self.maxCPUs[i] = self.mCPUs[0]
            else:
              self.maxCPUs[i] = self.mCPUs[1]

      # Choose the exceted nodes
      self.excetedNodes = []
      self.statusNodes  = []
      self.GetExceptedNodes()

      # Parser Object
      self.parser     = parser()

      # Process Object
      self.proc       = os.popen('qstat -f')

      # Number of jobs
      self.nTotJobs   = 0      # number of total   jobs
      self.nRunJobs   = 0      # number of running jobs
      self.nHolJobs   = 0      # number of holding jobs
      self.nWaiJobs   = 0      # number of Waiting jobs
      self.nEndJobs   = 0      # number of Ending  jobs

      # Number of CPUs
      self.nTotCPUs   = 0      # number of Total     CPUs
      self.nRunCPUs   = 0      # number of Running   CPUs
      self.nHolCPUs   = 0      # number of Holding   CPUs
      self.nWaiCPUs   = 0      # number of Waiting   CPUs
      self.nEndCPUs   = 0      # number of Waiting   CPUs

      self.nExcCPUs   = 0      # number of Exclusive CPUs
      self.nTotLines  = 0      # number of lines of 'qstat -f'

      # Job Information: ID, user, status <= dummy informations
      self.TotJobs = []
      self.RunJobs = []
      self.HolJobs = []
      self.WaiJobs = []
      self.EndJobs = []

      # Nodes info
      self.nCPUs_node = [0 for i in range(self.nNodes)]        # used cpus in each nodes
      self.Users_node = ['' for i in range(self.nNodes)]       # users in each nodes
      self.Exclu_node = [False for i in range(self.nNodes)]    # used cpus in each nodes

      # Available Nodes info
      self.nAvaNodes    = 0
      self.nCPUsAvaNode = []
      self.mCPUsAvaNode = []
      self.nameAvaNode  = []
      # Exclusive Nodes info
      self.nExcNodes    = 0
      self.nCPUsExcNode = []
      self.nameExcNode  = []

      ctime = time.localtime(time.time())
      c_year  = ctime.tm_year
      c_month = ctime.tm_mon
      c_day   = ctime.tm_mday
      c_hour  = ctime.tm_hour
      c_min   = ctime.tm_min
      c_sec   = ctime.tm_sec
      self.ctime = datetime.datetime(c_year, c_month, c_day, c_hour, c_min, c_sec)

   def GetExceptedNodes(self):
      par     = parser()
      pro     = os.popen('pbsnodes -a')
      par.initialize_proc(pro)
      par.make_blocked_string('Mom = ')

      for iblck in range(par.nblcks):
         stat = par.get_value_block(iblck, 'state = ', '\n').strip()
         serq = par.get_value_block(iblck, 'queue = ', '\n').strip()
         if stat.find('offline') >= 0:
            inode = par.get_value_block(iblck, 'Mom = node', '\n').strip()
            if int(inode) <= self.nNodes:
              self.excetedNodes.append(int(inode))
              self.statusNodes.append('offline')
         if stat.find('down') >= 0:
            inode = par.get_value_block(iblck, 'Mom = node', '\n').strip()
            #print stat, serq, inode, self.nNodes
            if inode=='-1: none': continue
            if int(inode) <= self.nNodes:
              self.excetedNodes.append(int(inode))
              self.statusNodes.append('down')
         if serq.find('serialq') >= 0:
            inode = par.get_value_block(iblck, 'Mom = node', '\n').strip()
            if int(inode) <= self.nNodes:
              self.excetedNodes.append(int(inode))
              self.statusNodes.append('serial queue')
         if serq.find('quickq') >= 0:
            inode = par.get_value_block(iblck, 'Mom = node', '\n').strip()
            if int(inode) <= self.nNodes:
              self.excetedNodes.append(int(inode))
              self.statusNodes.append('quick queue')
         if serq.find('special') >= 0:
            inode = par.get_value_block(iblck, 'Mom = node', '\n').strip()
            if int(inode) <= self.nNodes:
              self.excetedNodes.append(int(inode))
              self.statusNodes.append('special queue')



   def GetExceptedNodes_old(self):
      par     = parser()
      pro     = os.popen('qhosts')
      par.initialize_proc(pro)

      lineString  = par.string.splitlines()
      for iline in range(par.nlines):
         nodestr = par.get_value_line_pattern(lineString[iline], 'offline', 'node', 'offline')
         if nodestr != '':
            print int(nodestr)
            if int(nodestr) <= self.nNodes: self.excetedNodes.append(int(nodestr))

         nodestr = par.get_value_line_pattern(lineString[iline], 'serialq', 'node', '    ')
         if nodestr != '':
            if int(nodestr) <= self.nNodes: self.excetedNodes.append(int(nodestr))


   def set_job_info(self):
      #proc   = os.popen('qstat -f')
      parser = self.parser
      parser.initialize_proc(self.proc)
      parser.make_blocked_string('Job Id:')

      # number of ...
      self.nTotLines = parser.nlines
      self.nTotJobs  = parser.nblcks
      self.nRunJobs  = parser.get_counts('job_state = R')
      self.nHolJobs  = parser.get_counts('job_state = H')
      self.nWaiJobs  = parser.get_counts('job_state = Q')
      self.nEndJobs  = parser.get_counts('job_state = E')

      if (self.nTotJobs != self.nRunJobs+self.nHolJobs+self.nWaiJobs+self.nEndJobs):
         #print 'Warning: Occur "job-end" event.'
         print 'Need a synchronizing....  try again'
         #parser.print_string()
         #quit()

      # Job Information: ID, user, status, ...
      self.TotJobs   = [Jobs() for i in range(self.nTotJobs)]
      self.RunJobs   = [Jobs() for i in range(self.nRunJobs)]
      self.HolJobs   = [Jobs() for i in range(self.nHolJobs)]
      self.WaiJobs   = [Jobs() for i in range(self.nWaiJobs)]
      self.EndJobs   = [Jobs() for i in range(self.nEndJobs)]

      njobs   = self.nTotJobs
      irjob   = 0
      ihjob   = 0
      iqjob   = 0
      iejob   = 0
      for ijob in range(njobs):
         if self.system == 'Gaon2':
            self.TotJobs[ijob].id       = int(parser.get_value_block(ijob, 'Job Id:', '.gaon2_').strip())
         elif  self.system == 'Gaon1':
            self.TotJobs[ijob].id       = int(parser.get_value_block(ijob, 'Job Id:', '.master').strip())
         self.TotJobs[ijob].name     = parser.get_value_block(ijob, 'Job_Name =', '\n').strip()
         self.TotJobs[ijob].queue    = parser.get_value_block(ijob, 'queue =', '\n').strip()
         self.TotJobs[ijob].user     = parser.get_value_block(ijob, 'Job_Owner =', '@').strip()
         self.TotJobs[ijob].ncpus    = int(parser.get_value_block(ijob, 'Resource_List.ncpus =', '\n').strip())
         self.TotJobs[ijob].nnodes   = int(parser.get_value_block(ijob, 'Resource_List.nodect =', '\n').strip())
         exclu                       = parser.get_value_block(ijob, 'Resource_List.place =', '\n').strip()
         if exclu.find('excl') >= 0:
            self.TotJobs[ijob].exclu = True
         self.TotJobs[ijob].status   = parser.get_value_block(ijob, 'job_state =', '\n').strip()
         self.nTotCPUs = self.nTotCPUs + self.TotJobs[ijob].ncpus

         if self.TotJobs[ijob].status == 'R':
            # node information
            node_txt = parser.get_value_block(ijob, 'exec_vnode =', 'Hold_Types =').strip()
            # DEBUG
            #print node_txt
            self.SaveNodeInfo(self.TotJobs[ijob], node_txt)  # save 'node', 'ncpus_n' attributes
            # walltime information
            start_txt = parser.get_value_block(ijob, 'stime = ', '\n').strip()
            self.SaveWallInfo(self.TotJobs[ijob], start_txt)  # save 'walltime'


            self.RunJobs[irjob].CopyJob(self.TotJobs[ijob])
            irjob = irjob + 1

            self.nRunCPUs = self.nRunCPUs + self.TotJobs[ijob].ncpus
         elif self.TotJobs[ijob].status == 'H':
            self.HolJobs[ihjob].CopyJob(self.TotJobs[ijob])
            ihjob = ihjob + 1

            self.nHolCPUs = self.nHolCPUs + self.TotJobs[ijob].ncpus
         elif self.TotJobs[ijob].status == 'Q':
            self.WaiJobs[iqjob].CopyJob(self.TotJobs[ijob])
            iqjob = iqjob + 1

            self.nWaiCPUs = self.nWaiCPUs + self.TotJobs[ijob].ncpus
         elif self.TotJobs[ijob].status == 'E':
            # node information
            node_txt = parser.get_value_block(ijob, 'exec_vnode =', 'Hold_Types =').strip()
            self.SaveNodeInfo(self.TotJobs[ijob], node_txt)  # save 'node', 'ncpus_n' attributes
            # walltime information
            start_txt = parser.get_value_block(ijob, 'stime = ', '\n').strip()
            self.SaveWallInfo(self.TotJobs[ijob], start_txt)  # save 'walltime'


            self.EndJobs[iejob].CopyJob(self.TotJobs[ijob])
            iejob = iejob + 1

            self.nEndCPUs = self.nEndCPUs + self.TotJobs[ijob].ncpus
         else:
            print 'could not determined job status...'
#            quit()

      self.SaveAvaNodeInfo()


   def SaveNodeInfo(self, job, nodetxt):
      parser  = self.parser
      if self.system == 'Gaon1':
         s_node  = '(node'
         e_node  = '.cm.cluster:ncpus'
         s_cpu   = 'ncpus='
         e_cpu   = ')'
      elif self.system == 'Gaon2':
         s_node  = '(node'
         e_node  = ':ncpus'
         s_cpu   = 'ncpus='
         e_cpu   = ')'
      nodetxt = nodetxt.replace('\n','')
      nodetxt = nodetxt.replace('\t','')

      nnodes  = job.nnodes
      i = 0
      for i in range(nnodes):
         inode = int(parser.get_value(nodetxt, s_node, e_node, i))
         ncpus = int(parser.get_value(nodetxt, s_cpu, e_cpu, i))
         #print inode
         job.nodes.append(inode)
         job.ncpus_n.append(ncpus)
         # DEBUG
         #print ncpus
         self.nCPUs_node[inode-1] = self.nCPUs_node[inode-1] + ncpus
         self.Users_node[inode-1] = self.Users_node[inode-1] + ', ' + job.user
         self.Exclu_node[inode-1] = job.exclu


   def SaveWallInfo(self, job, starttxt):
      months = {'Jan':1, 'Feb':2, 'Mar':3, 'Apr':4, 'May':5, 'Jun':6, 'Jul':7, 'Aug':8, 'Sep':9, 'Oct':10, 'Nov':11, 'Dec':12,}
      splittxt = starttxt.split()

      s_year  = int(splittxt[4])
      s_day   = int(splittxt[2])
      s_month = months[splittxt[1]]
      timestring = splittxt[3].split(':')
      s_hour  = int(timestring[0])
      s_min   = int(timestring[1])
      s_sec   = int(timestring[2])
      stime = datetime.datetime(s_year, s_month, s_day, s_hour, s_min, s_sec)
      dtime = self.ctime - stime
      job.walltime[0] = dtime.total_seconds()

      min  = 0
      hour = 0
      day  = 0

      min  = int(job.walltime[0])/60
      sec  = int(job.walltime[0])%60
      if min > 0:
         hour = min/60
         min  = min%60
         if hour > 0:
            day  = hour/24
            hour = hour%24
      job.walltime[1] = day
      job.walltime[2] = hour
      job.walltime[3] = min
      job.walltime[4] = sec
    


   def SaveAvaNodeInfo(self):
      maxcpus = self.maxCPUs
      nnodes  = self.nNodes
      for inode in range(nnodes):
         if self.nCPUs_node[inode] != maxcpus[inode]:
            if self.excetedNodes.count(inode+1) == 0:
               if self.Exclu_node[inode] == True:
                 self.nameExcNode.append(str(inode+1))
                 self.nCPUsExcNode.append(maxcpus[inode] - self.nCPUs_node[inode])
                 self.nExcNodes = self.nExcNodes + 1
               else:
                 self.nameAvaNode.append(str(inode+1))
                 self.nCPUsAvaNode.append(maxcpus[inode] - self.nCPUs_node[inode])
                 self.mCPUsAvaNode.append(maxcpus[inode])
                 self.nAvaNodes = self.nAvaNodes + 1
            else:
               self.nameExcNode.append(str(inode+1))
               self.nCPUsExcNode.append(maxcpus[inode] - self.nCPUs_node[inode])
               self.nExcNodes = self.nExcNodes + 1


   def set_user_info(self):
      parser     = self.parser
      self.Users = []
      njobs      = self.nTotJobs

      users = []
      ncpus = []
      state = []
      for ijob in range(njobs):
         users.append(self.TotJobs[ijob].user)
         ncpus.append(self.TotJobs[ijob].ncpus)
         state.append(self.TotJobs[ijob].status)

      nusers = 0
      for ijob in range(njobs):
         isPre = False
         for iuser in range(nusers):
            if self.Users[iuser].name == users[ijob]:
               self.Users[iuser].nTotJobs = self.Users[iuser].nTotJobs + 1
               self.Users[iuser].TotCPUs  = self.Users[iuser].TotCPUs + ncpus[ijob]
               if state[ijob] == 'R':
                  self.Users[iuser].nRunJobs = self.Users[iuser].nRunJobs + 1
                  self.Users[iuser].RunCPUs  = self.Users[iuser].RunCPUs + ncpus[ijob]
               elif state[ijob] == 'H':
                  self.Users[iuser].nHolJobs = self.Users[iuser].nHolJobs + 1
                  self.Users[iuser].HolCPUs  = self.Users[iuser].HolCPUs + ncpus[ijob]
               elif state[ijob] == 'Q':
                  self.Users[iuser].nWaiJobs = self.Users[iuser].nWaiJobs + 1
                  self.Users[iuser].WaiCPUs  = self.Users[iuser].WaiCPUs + ncpus[ijob]
               elif state[ijob] == 'E':
                  self.Users[iuser].nEndJobs = self.Users[iuser].nEndJobs + 1
                  self.Users[iuser].EndCPUs  = self.Users[iuser].EndCPUs + ncpus[ijob]
               isPre = True
         if not isPre:
            userobj = Users()
            userobj.id = nusers
            userobj.name     = users[ijob]
            userobj.nTotJobs = 1
            userobj.TotCPUs  = ncpus[ijob]
            if state[ijob] == 'R':
               userobj.nRunJobs = 1
               userobj.RunCPUs  = ncpus[ijob]
            elif state[ijob] == 'H':
               userobj.nHolJobs = 1
               userobj.HolCPUs  = ncpus[ijob]
            elif state[ijob] == 'Q':
               userobj.nWaiJobs = 1
               userobj.WaiCPUs  = ncpus[ijob]
            elif state[ijob] == 'E':
               userobj.nEndJobs = 1
               userobj.EndCPUs  = ncpus[ijob]
            self.Users.append(userobj)
            nusers = nusers + 1

      self.nUsers = nusers



   def print_qstat(self):
      p = os.popen('qstat | wc -l')
      output = p.readline()
      nRunJobs = int(output)
      
      p = os.popen('qstat')
      print ' '
      for i in range(nRunJobs):
         output = p.readline()
         jobinfo = output.replace('\n','')
         print jobinfo
      print ' '



   def print_status(self, isUser=False):
      nnodes  = self.nNodes
      maxcpus = self.maxCPUs
      print ' '
      print '# Nodes information'
      print '#  - * : exclusive node'
      print '#  - % : excepted  node (serial job, checking, repairing, ...)'
      print '============================================='
      print '|   Node   |    max   |   used   |   avail  |'
      print '|-------------------------------------------|'
      for inode in range(nnodes):
         self.Users_node[inode] = self.Users_node[inode].replace(', ', '', 1)
         if self.Exclu_node[inode] == False:
           tag_str = ' '
         else:
           tag_str = '*'

         if self.excetedNodes.count(inode+1) == 0:
            tag_str = tag_str+' '
         else:
            tag_str = tag_str+'%'
         if isUser:
            print '|  node%003d |      %2d  | %2s   %2d  |      %2d  |' % \
             (inode+1, maxcpus[inode], tag_str, self.nCPUs_node[inode], maxcpus[inode]-self.nCPUs_node[inode]), self.Users_node[inode]
         else:
            print '|  node%003d |      %2d  | %2s   %2d  |      %2d  |' % \
             (inode+1, maxcpus[inode], tag_str, self.nCPUs_node[inode], maxcpus[inode]-self.nCPUs_node[inode])
      print '|-------------------------------------------|'
      print '|    Total |    %4d  |    %4d  |    %4d  |' % \
        (sum(maxcpus), sum(self.nCPUs_node[:]), sum(maxcpus)-sum(self.nCPUs_node[:]))
      print '============================================='
      print ' '


   def print_users(self):
      nusers = self.nUsers
      print ' '
      print '# Users information'
      print '======================================================================================================'
      print '|    User  | nTJobs | nRJobs | nEJobs | nHJobs | nWJobs |  TCPUs |  RCPUs |  ECPUs |  HCPUs |  WCPUs |'
      print '|----------------------------------------------------------------------------------------------------|'
      for iuser in range(nusers):
         print '| %8s |    %3d |    %3d |    %3d |    %3d |    %3d |   %4d |   %4d |   %4d |   %4d |   %4d |' % \
           (self.Users[iuser].name, self.Users[iuser].nTotJobs, self.Users[iuser].nRunJobs, self.Users[iuser].nEndJobs, self.Users[iuser].nHolJobs, self.Users[iuser].nWaiJobs, self.Users[iuser].TotCPUs, self.Users[iuser].RunCPUs, self.Users[iuser].EndCPUs, self.Users[iuser].HolCPUs, self.Users[iuser].WaiCPUs)
      print '|----------------------------------------------------------------------------------------------------|'
      print '|   Total  |    %3d |    %3d |    %3d |    %3d |    %3d |   %4d |   %4d |   %4d |   %4d |   %4d |' % \
          (self.nTotJobs, self.nRunJobs, self.nEndJobs, self.nHolJobs, self.nWaiJobs, self.nTotCPUs, self.nRunCPUs, self.nEndCPUs, self.nHolCPUs, self.nWaiCPUs)
      print '======================================================================================================'
      print ' '
      print '    - nTotXXXs : # of   total jobs or cpus in jobs.'
      print '    - nRunXXXs : # of  runnig jobs or cpus in jobs. (job state: R)'
      print '    - nHolXXXs : # of holding jobs or cpus in jobs. (job state: H)'
      print '    - nWaiXXXs : # of waiting jobs or cpus in jobs. (job state: W)'
      print '    - nEndXXXs : # of  ending jobs or cpus in jobs. (job state: E)'
      print ' '


   def print_jobs_cpus(self):
      print ' '
      print '# Jobs & CPUs information'
      print '======================================================================================================'
      print '|    User  | nTJobs | nRJobs | nEJobs | nHJobs | nWJobs |  TCPUs |  RCPUs |  ECPUs |  HCPUs |  WCPUs |'
      print '|----------------------------------------------------------------------------------------------------|'
      print '|   Total  |    %3d |    %3d |    %3d |    %3d |    %3d |   %4d |   %4d |   %4d |   %4d |   %4d |' % \
            (self.nTotJobs, self.nRunJobs, self.nEndJobs, self.nHolJobs, self.nWaiJobs, self.nTotCPUs, self.nRunCPUs, self.nEndCPUs, self.nHolCPUs, self.nWaiCPUs)
      print '======================================================================================================'
      print ' '
      print '    - nTotXXXs : # of   total jobs or cpus in jobs.'
      print '    - nRunXXXs : # of  runnig jobs or cpus in jobs. (job state: R)'
      print '    - nHolXXXs : # of holding jobs or cpus in jobs. (job state: H)'
      print '    - nWaiXXXs : # of waiting jobs or cpus in jobs. (job state: W)'
      print '    - nEndXXXs : # of  ending jobs or cpus in jobs. (job state: E)'
      print ' '


   def print_jobs(self):
      nrjobs = self.nRunJobs
      nejobs = self.nEndJobs
      print ' '
      print '# Job information (for the running and ending jobs)'
      print '=================================================================================='
      print '|    ID   |      User    |      Name    |     Queue    |  nCores  |   Wall-time  |'
      print '|---------------------------------------------------------------------------------'
      for ijob in range(nrjobs):
         print '| %7d |  %10s  |  %10s  |  %10s  |    %4d  | %2dd %02d:%02d:%02d |' %\
            (self.RunJobs[ijob].id, self.RunJobs[ijob].user, self.RunJobs[ijob].name[0:10], self.RunJobs[ijob].queue, self.RunJobs[ijob].ncpus, \
             self.RunJobs[ijob].walltime[1], self.RunJobs[ijob].walltime[2], self.RunJobs[ijob].walltime[3], self.RunJobs[ijob].walltime[4])
      for ijob in range(nejobs):
         print '| %7d |  %10s  |  %10s  |  %10s  |    %4d  | %2dd %02d:%02d:%02d |' %\
            (self.EndJobs[ijob].id, self.EndJobs[ijob].user, self.EndJobs[ijob].name[0:10], self.EndJobs[ijob].queue, self.EndJobs[ijob].ncpus, \
             self.EndJobs[ijob].walltime[1], self.EndJobs[ijob].walltime[2], self.EndJobs[ijob].walltime[3], self.EndJobs[ijob].walltime[4])
      print '=================================================================================='


   def print_nodes_cpus(self):
      nnodes  = self.nNodes
      maxcpus = self.maxCPUs

      enodes  = self.nExcNodes
      anodes  = self.nAvaNodes
      rnodes  = nnodes - enodes - anodes

      ncpus   = sum(maxcpus)
      rcpus   = self.nRunCPUs + self.nEndCPUs
      ecpus   = sum(self.nCPUsExcNode)
      acpus   = sum(self.nCPUsAvaNode)
      print ' '
      print '# Total Statistic'
      print '===================================================================================================='
      print '|          | TotNodes | RunNodes | ExcNodes | AvaNodes |  TotCPUs |  RunCPUs |  ExcCPUs |  AvaCPUs |'
      print '|--------------------------------------------------------------------------------------------------|'
      print '|   Total  |    %4d  |    %4d  |    %4d  |    %4d  |    %4d  |    %4d  |    %4d  |    %4d  |' %    \
            (nnodes, rnodes, enodes, anodes, ncpus, rcpus, ecpus, acpus)
      print '===================================================================================================='
      print ' '
      print '    - nTotXXXXs : # of     total nodes or cpus in the machine.'
      print '    - nRunXXXXs : # of    runnig nodes or cpus in the machine. (include ending event)'
      print '    - nExcXXXXs : # of exclusive nodes or cpus in the machine. (include excepted nodes)'
      print '    - nAvaXXXXs : # of available nodes or cpus in the machine.'
      print ' '



   def print_available(self):
      navanode = self.nAvaNodes
      #print self.nAvaNodes, self.nCPUsAvaNode, self.nameAvaNode
      nsort, ncpus_n, nnodes = self.SortSumList(self.nAvaNodes, self.nCPUsAvaNode, self.mCPUsAvaNode)
      #print nsort, ncpus_n, nnodes
      print ' '
      print '# Available maximum nodes & cpus combinations (will be updated)'
      print '  % excepted nodes (serialq or repairing, ...):', self.excetedNodes, self.statusNodes
      print '============================================='
      print '|          |  nNodes  | CPUs/node|   nCPUs  |'
      print '|          | (select=)| (ncpus=) |  (-n XX) |'
      print '|-------------------------------------------|'
      for i in range(nsort):
        print '|   %3d    |    %4d  |    %4d  |    %4d  |' % \
              (i+1, nnodes[i], ncpus_n[i], nnodes[i]*ncpus_n[i])
      print '============================================='
      print ' '
      if nsort > 0:
         print '  % example'
         string = '    #PBS -l select='+str(nnodes[0])+':ncpus='+str(ncpus_n[0])+':mpiprocs='+str(ncpus_n[0])+':ompthreads=1'
         print string
         string = '    ...'
         print string
         string = '    mpirun -n '+str(nnodes[0]*ncpus_n[0])+' -f $PBS_NODEFILE ./my_program_name' 
         print string
         print ' '


   def SortSumList(self, nAvaNodes, nCPUsAvaNode, mCPUsAvaNode):
      ncounts     = 0
      res_nnodes  = []
      res_ncpus_n = []
      res_mcpus_n = []


      # 같은 개수의 cpu끼리 묶음 [2,2,2,4,4] => 2, [3,2], [2,4]
      for ilist in range(nAvaNodes):
         isPre = False
         for idiff in range(ncounts):
            if res_ncpus_n[idiff] == nCPUsAvaNode[ilist]:
               res_nnodes[idiff] = res_nnodes[idiff] + 1
               isPre = True
         if not isPre:
            res_nnodes.append(1)
            res_ncpus_n.append(nCPUsAvaNode[ilist])
            res_mcpus_n.append(mCPUsAvaNode[ilist])
            ncounts = ncounts + 1


      # max cpu에 대한 오름차순 정리
      tmp = []
      for ii in range(ncounts):
         tmp.append([res_mcpus_n[ii], res_ncpus_n[ii], res_nnodes[ii]])
      tmp = sorted(tmp, key=self.getKey)
      for ii in range(ncounts):
         res_mcpus_n[ii] = tmp[ii][0]
         res_ncpus_n[ii] = tmp[ii][1]
         res_nnodes[ii] = tmp[ii][2]


      # Queue 구분
      nDivs  = len(self.mCPUs)
      staDiv = []
      couDiv = []
      for idiv in range(nDivs):
         if res_mcpus_n.count(self.mCPUs[idiv]) > 0:
            staDiv.append(res_mcpus_n.index(self.mCPUs[idiv]))
            couDiv.append(res_mcpus_n.count(self.mCPUs[idiv]))
         else:
            staDiv.append(0)
            couDiv.append(0)
  

      # 오름차순 정리 2, [3,2], [2,4] => 2, [5,2], [2,4]
      for idiv in range(nDivs):
         tmp = []
         sta = staDiv[idiv]
         end = staDiv[idiv]+couDiv[idiv]
         for ii in range(sta,end):
            tmp.append([res_ncpus_n[ii], res_nnodes[ii]])
         tmp = sorted(tmp, key=self.getKey)
         for ii in range(sta,end):
            res_ncpus_n[ii] = tmp[ii-sta][0]
            res_nnodes[ii] = tmp[ii-sta][1]
            for i in range(sta+ii+1,end):
               res_nnodes[ii] = res_nnodes[ii] + tmp[i-sta][1]

      
      return ncounts, res_ncpus_n, res_nnodes


   def getKey(self, item):
      return item[0]



   def print_status_eASCII(self, isUser=False):
      print ' '
      print '# Nodes information'
      print '╔══════════╦══════════╦══════════╦══════════╗'
      print '║   Node   │    max   │   used   │   avail  ║'
      print '╠──────────┼──────────┼──────────┼──────────╣'
      nnodes  = self.nNodes
      maxcpus = self.maxCPUs
      for inode in range(nnodes):
         self.Users_node[inode] = self.Users_node[inode].replace(', ', '', 1)
         if isUser:
            print '║  node%003d │      %2d  │      %2d  │      %2d  ║' % \
            (inode+1, maxcpus[inode], self.nCPUs_node[inode], maxcpus[inode]-self.nCPUs_node[inode]), self.Users_node[inode]
         else:
            print '║  node%003d │      %2d  │      %2d  │      %2d  ║' % \
             (inode+1, maxcpus[inode], self.nCPUs_node[inode], maxcpus[inode]-self.nCPUs_node[inode])
      print '╠──────────┼──────────┼──────────┼──────────╣'
      print '║    Total │    %4d  │    %4d  │    %4d  ║' % \
        (sum(maxcpus), sum(self.nCPUs_node[:]), sum(maxcpus)-sum(self.nCPUs_node[:]))
      print '╚══════════╩══════════╩══════════╩══════════╝'
      print ' '


   def print_users_eASCII(self):
      nusers = self.nUsers
      print ' '
      print '# Users information'
      print '╔══════════╦══════════╦══════════╦══════════╦══════════╦══════════╦══════════╦══════════╦══════════╗'
      print '║    User  ║ nTotJobs │ nRunJobs │ nHolJobs │ nWaiJobs ║  TotCPUs │  RunCPUs │  HolCPUs │  WaiCPUs ║'
      print '╠──────────╬──────────┼──────────┼──────────┼──────────╬──────────┼──────────┼──────────┼──────────╣'
      for iuser in range(nusers):
         print '║ %8s ║      %2d  │      %2d  │      %2d  │      %2d  ║    %4d  │    %4d  │    %4d  │    %4d  ║' % \
           (self.Users[iuser].name, self.Users[iuser].nTotJobs, self.Users[iuser].nRunJobs, self.Users[iuser].nHolJobs, self.Users[iuser].nWaiJobs, self.Users[iuser].TotCPUs, self.Users[iuser].RunCPUs, self.Users[iuser].HolCPUs, self.Users[iuser].WaiCPUs)
      print '╠──────────╬──────────┼──────────┼──────────┼──────────╬──────────┼──────────┼──────────┼──────────╣'
      print '║   Total  ║      %2d  │      %2d  │      %2d  │      %2d  ║    %4d  │    %4d  │    %4d  │    %4d  ║' % \
          (self.nTotJobs, self.nRunJobs, self.nHolJobs, self.nWaiJobs, self.nTotCPUs, self.nRunCPUs, self.nHolCPUs, self.nWaiCPUs)
      #print '╚==================================================================================================╝'
      print '╚══════════╩══════════╩══════════╩══════════╩══════════╩══════════╩══════════╩══════════╩══════════╝'
      print ' '


   def print_jobs_cpus_eASCII(self):
      print ' '
      print '# Jobs & CPUs information'
      print '╔══════════╦══════════╦══════════╦══════════╦══════════╦══════════╦══════════╦══════════╦══════════╗'
      print '║          ║ nTotJobs │ nRunJobs │ nHolJobs │ nWaiJobs ║  TotCPUs │  RunCPUs │  HolCPUs │  WaiCPUs ║'
      print '╠──────────╬──────────┼──────────┼──────────┼──────────╬──────────┼──────────┼──────────┼──────────╣'
      print '║   Total  ║      %2d  │      %2d  │      %2d  │      %2d  ║    %4d  │    %4d  │    %4d  │    %4d  ║' % \
            (self.nTotJobs, self.nRunJobs, self.nHolJobs, self.nWaiJobs, self.nTotCPUs, self.nRunCPUs, self.nHolCPUs, self.nWaiCPUs)
      print '╚══════════╩══════════╩══════════╩══════════╩══════════╩══════════╩══════════╩══════════╩══════════╝'
      print ' '


   def print_nodes_cpus_eASCII(self):
      nnodes  = self.nNodes
      maxcpus = self.maxCPUs
      print ' '
      print '# Total Statistic'
      print '╔══════════╦══════════╦══════════╦══════════╦══════════╦══════════╦══════════╗'
      print '║          ║ TotNodes │ RunNodes │ AvaNodes ║  TotCPUs │  RunCPUs │  AvaCPUs ║'
      print '╠──────────╬──────────┼──────────┼──────────╬──────────┼──────────┼──────────╣'
      print '║   Total  ║    %4d  │    %4d  │    %4d  ║    %4d  │    %4d  │    %4d  ║' % \
            (self.nNodes, self.nNodes-self.nAvaNodes, self.nAvaNodes, nnodes*maxcpus[inode], self.nRunCPUs, nnodes*maxcpus[inode]-self.nRunCPUs)
      print '╚══════════╩══════════╩══════════╩══════════╩══════════╩══════════╩══════════╝'
      print ' '



   def print_available_eASCII(self):
      navanode = self.nAvaNodes
      #print self.nAvaNodes, self.nCPUsAvaNode, self.nameAvaNode
      nsort, values, ncounts = self.SortSumList(self.nAvaNodes, self.nCPUsAvaNode)
      #print nsort, values, ncounts
      print ' '
      print '# Available maximum nodes & cpus combinations'
      print '  % excepted nodes:', self.excetedNodes
      print '╔══════════╦══════════╦══════════╦══════════╗'
      print '║          ║  nNodes  │ CPUs/node│   nCPUs  ║'
      print '║          ║ (select=)│ (ncpus=) │  (-n XX) ║'
      print '╠──────────╬──────────┼──────────┼──────────╣'
      for i in range(nsort):
        print '║   %3d    ║    %4d  │    %4d  │    %4d  ║' % \
              (i+1, ncounts[i], values[i], ncounts[i]*values[i])
      print '╚══════════╩══════════╩══════════╩══════════╝'
      print ' '



# old version


   def print_nodes_cpus_old(self):
      nnodes  = self.nNodes
      maxcpus = self.maxCPUs
      print ' '
      print '# Total Statistic'
      print '╔══════════╦══════════╦══════════╦══════════╦══════════╦══════════╦══════════╦══════════╦══════════╦══════════╦══════════╗'
      print '║          ║ nTotJobs │ nRunJobs │ nHolJobs │ nWaiJobs ║  TotCPUs │  RunCPUs │  HolCPUs │  WaiCPUs ║  SysCPUs │  AvaCPUs ║'
      print '╠──────────╬──────────┼──────────┼──────────┼──────────╬──────────┼──────────┼──────────┼──────────╬──────────┼──────────╣'
      print '║   Total  ║      %2d  │      %2d  │      %2d  │      %2d  ║    %4d  │    %4d  │    %4d  │    %4d  ║    %4d  │    %4d  ║' % \
            (self.nTotJobs, self.nRunJobs, self.nHolJobs, self.nWaiJobs, self.nTotCPUs, self.nRunCPUs, self.nHolCPUs, self.nWaiCPUs, nnodes*maxcpus[inode], nnodes*maxcpus[inode]-self.nRunCPUs)
      print '╚══════════╩══════════╩══════════╩══════════╩══════════╩══════════╩══════════╩══════════╩══════════╩══════════╩══════════╝'
      print ' '



