# -*- coding: utf-8 -*- import sys
import sys
sys.path.append('/home/jhkim/Study/Library/Shared/Python')
from LibParser import *
from LibPlot import *
#######################################################################################
# Description: Profiling of KIM
#
# Date: 13MAR2015
#   - JH KIM: First, for KIM version 0.21.01
#
# Date: 18MAR2015
#   - JH KIM: KIM version 0.21.02
#             Add Subroutine: WrtieFile, ReadFile
#
# Date: 24MAR2015
#   - JH KIM: KIM version 0.21.02
#             Modified for Difference and formats
#
# Date 03APR2015
#   - JH KIM: Initailize with input path (Initialize_DIR, Dir2CoresDirs)
#
#######################################################################################

class Component:

   def __init__(self, compname_in, names_in):
      self.compName  = compname_in

      self.nNames    = len(names_in)+1
      self.Names     = names_in
      self.Names.append('Total')

      self.nKeys     = [0 for i in range(self.nNames)]
      self.Keys      = None
      self.Keys_Sign = None

      self.isIni     = False


   def CopyComponent(self, comp):
      self.compName  = comp.compName
      self.nNames    = comp.nNames
      self.Names     = comp.Names
      self.nKeys     = comp.nKeys
      self.Keys      = comp.Keys
      self.Keys_Sign = comp.Keys_Sign
      self.nCores    = comp.nCores
      self.Cores     = comp.Cores
      self.Walltime  = comp.Walltime
      self.Mean      = comp.Mean
      self.Dev       = comp.Dev
      self.SpeedBase = comp.SpeedBase
      self.Speed     = comp.Speed
      self.DevR      = comp.DevR
      self.Ratio     = comp.Ratio


   def CopyMeta(self, comp):
      #self.compName  = comp.compName
      #self.nNames    = comp.nNames
      #self.Names     = comp.Names
      self.nKeys     = comp.nKeys
      self.Keys      = comp.Keys
      self.Keys_Sign = comp.Keys_Sign


   def SetKeys(self, nKeys_in, Keys_in, Sign_in):
      if len(nKeys_in) != self.nNames-1:
         print 'ERROR :: check the dimension of sub components... :', self.compName
         quit()
      self.nKeys = nKeys_in
      self.nKeys.append(0)
      self.Keys  = Keys_in
      self.Keys.append([])
      self.Keys_Sign = Sign_in
      self.Keys_Sign.append([])
      print '## Loading %13s ...  has %2d components.'%('[ '+self.compName+' ]', self.nNames)
      #self.PrintKeysInfo()


   def PrintKeysInfo(self):

      print 'Key Info:'
      for i in range(self.nNames):
         print '%2i : %s'%(i, self.Names[i])
         for j in range(self.nKeys[i]):
            print '  Key%2i : %s (sign: %2d)'%(j, self.Keys[i][j], self.Keys_Sign[i][j])
      #print ' '


   def IniCores(self, cores_in, isDiff=False):
      if not isDiff:
         self.isIni     = True
         self.nCores    = len(cores_in)
         self.Cores     = cores_in     # x-axis
   
         self.Walltime  = [[[0.0 for k in range(self.Cores[j])] for j in range(self.nCores)] for i in range(self.nNames)]
         self.Mean      = [[0.0 for j in range(self.nCores)] for i in range(self.nNames)]
         self.Dev       = [[0.0 for j in range(self.nCores)] for i in range(self.nNames)]
         self.SpeedBase = 0.0
         self.Speed     = [[0.0 for j in range(self.nCores)] for i in range(self.nNames)]
         self.DevR      = [[0.0 for j in range(self.nCores)] for i in range(self.nNames)]
         self.Ratio     = [[0.0 for j in range(self.nCores)] for i in range(self.nNames)]
      else:
         if not self.isIni:
            print 'component must be initialized...'
            quit()
         self.nCoresD    = len(cores_in)
         self.CoresD     = cores_in     # x-axis
   
         self.WalltimeD  = [[[0.0 for k in range(self.CoresD[j])] for j in range(self.nCoresD)] for i in range(self.nNames)]
         self.MeanD      = [[0.0 for j in range(self.nCoresD)] for i in range(self.nNames)]
         self.DevD       = [[0.0 for j in range(self.nCoresD)] for i in range(self.nNames)]
         self.SpeedBaseD = 0.0
         self.SpeedD     = [[0.0 for j in range(self.nCoresD)] for i in range(self.nNames)]
         self.DevRD      = [[0.0 for j in range(self.nCoresD)] for i in range(self.nNames)]
         self.RatioD     = [[0.0 for j in range(self.nCoresD)] for i in range(self.nNames)]

         # mapping
         self.mapD       = [-1 for i in range(self.nCoresD)]
         for i in range(self.nCoresD):
            for j in range(self.nCores):
               if self.CoresD[i] == self.Cores[j]:
                  self.mapD[i] = j
                  continue


   def ReadData(self, filedirs, filehead, filetail, isDiff=False, is1File=False):
      if not isDiff:
         nCores = self.nCores
         Cores  = self.Cores
      else:
         nCores = self.nCoresD
         Cores  = self.CoresD

      if len(filedirs) != nCores:
         print 'check the directories of input files...'
         print filedirs
         quit()

      print '## Read Profile data.... :', self.compName
      for ncores in range(nCores):
         filename = filedirs[ncores]+'/'+filehead+'Summary'+filetail
         if is1File: filename = filedirs[ncores]+'/'+filehead+filetail
         print '  - START Parsing (nprocs: %6d) Directory : %s'%(Cores[ncores], filename)
         infile = open(filename, 'r')
         self.Parsing(infile, ncores, isDiff=isDiff, is1File=is1File)
         infile.close()
         print '  - END   Parsing (nprocs: %6d) Directory : %s'%(Cores[ncores], filename)
         print ' '
      self.Statistics(isSum=True,isDiff=isDiff)


   def ReadData_files(self, filedirs, filehead, filetail):

      if len(filedirs) != self.nCores:
         print 'check the directories of input files...'
         quit()

      print '## Read Profile data.... :', self.compName
      for ncores in range(self.nCores):
         print '  - START Parsing (nprocs: %6d) Directory : %s'%(self.Cores[ncores], filename)
         for icore in range(self.Cores[ncores]):
            filename = filedirs[ncores]+'/'+filehead+'_%06i'%(icore)+filetail
            infile = open(filename, 'r')
            self.Parsing(infile, ncores, icore)
            infile.close()
         print '  - END   Parsing (nprocs: %6d) Directory : %s'%(self.Cores[ncores], filename)
         print ' '
      #self.PrintIniInfo()
      self.Statistics(isSum=False,isDiff=False)


   def Parsing(self, infile_in, ncores_in, icore_in=-1, isDiff=False, is1File=False):
      if not is1File:
         if icore_in < 0:
            string = infile_in.read()
            idx = 2
         else:
            string = infile_in.read(17600)   # 17000 - skip the bellow message
            idx = 5
      else:
         #string = infile_in.read(17600)   # 17000 - skip the bellow message
         string = infile_in.read()   # 17000 - skip the bellow message
         string = string.split('************ PROCESS      0 (     0) ************')[1]
         idx = 3
         #print string
      #if string.find('Run HybDA') >= 0: idx = 3
      parser = Parser(string)
      for i in range(self.nNames-1):
         for j in range(self.nKeys[i]):
            n = parser.GetCounts(self.Keys[i][j])
            if n < 1:
               print 'check key: ', n, self.Keys[i][j]
               print 'in ', ncores_in, icore_in
               quit()
            linestring = parser.GetValueString(self.Keys[i][j], '\n')
            linestring = linestring.split()
            if icore_in < 0:
               if not isDiff:
                  self.Mean[i][ncores_in] = self.Mean[i][ncores_in] + self.Keys_Sign[i][j]*float(linestring[idx])
               else:
                  self.MeanD[i][ncores_in] = self.MeanD[i][ncores_in] + self.Keys_Sign[i][j]*float(linestring[idx])
            else:
               if not isDiff:
                  self.Walltime[i][ncores_in][icore_in] = self.Walltime[i][ncores_in][icore_in] + self.Keys_Sign[i][j]*float(linestring[idx])
               else:
                  self.WalltimeD[i][ncores_in][icore_in] = self.WalltimeD[i][ncores_in][icore_in] + self.Keys_Sign[i][j]*float(linestring[idx])


   def PrintIniInfo(self):
      print 'Initailize %13s ...  has %2d components.'%('[ '+self.compName+' ]', self.nNames)
      for i in range(self.nNames):
         #print '    %2i : %05d'%(i, self.nCores)
         for j in range(self.nCores):
            print self.Names[i], self.Walltime[i][j]
      print ' '


   def Statistics(self,isSum=True,isDiff=False):
      import numpy as np
      if not isDiff:
         nCores = self.nCores
      else:
         nCores = self.nCoresD

      # Save Average, Standard deviation
      if not isSum:
         for iname in range(self.nNames-1):
            for ncore in range(nCores):
               self.Mean[iname][ncore] = np.mean(self.Walltime[iname][ncore][:]).tolist()
               self.Dev[iname][ncore]  = np.std(self.Walltime[iname][ncore][:]).tolist()

      if not isDiff:
         # Save 'Total'
         for ncore in range(nCores):
            for iname in range(self.nNames-1):
               self.Mean[self.nNames-1][ncore] = self.Mean[self.nNames-1][ncore] + self.Mean[iname][ncore]
   
         # Save Speedup, Ratio
         for iname in range(self.nNames):
            base = self.Mean[iname][0]*self.Cores[0]
            for ncore in range(nCores):
               if self.Mean[iname][ncore] != 0:
                  self.Speed[iname][ncore] = base/self.Mean[iname][ncore]
                  self.Ratio[iname][ncore] = self.Mean[iname][ncore]/self.Mean[self.nNames-1][ncore]*100.0
               else:
                  self.Speed[iname][ncore] = 0.0
                  self.Ratio[iname][ncore] = 0.0
      else:
         # Save 'Total'
         for ncore in range(nCores):
            for iname in range(self.nNames-1):
               self.MeanD[self.nNames-1][ncore] = self.MeanD[self.nNames-1][ncore] + self.MeanD[iname][ncore]
   
         # Save Speedup, Ratio
         for iname in range(self.nNames):
            base = self.MeanD[iname][0]*self.Cores[0]
            for ncore in range(nCores):
               if self.MeanD[iname][ncore] != 0:
                  self.SpeedD[iname][ncore] = base/self.MeanD[iname][ncore]
                  self.RatioD[iname][ncore] = self.MeanD[iname][ncore]/self.MeanD[self.nNames-1][ncore]*100.0
               else:
                  self.SpeedD[iname][ncore] = 0.0
                  self.RatioD[iname][ncore] = 0.0


   def PrintScreen(self):
      # self.Mean, self.Ratio, self. Speed
      print ' '

      print '[[  %12s  ]]'%{self.compName}
      string = ' nprocs  '
      for iname in range(self.nNames):
         string = string + self.Names[iname] + '   '
      print string

      for ncore in range(self.nCores):
         string = ' %04d |'%(self.Cores[ncore])
         for iname in range(self.nNames):
            #string = string + str(self.Mean[iname][ncore]) + ' ' + str(self.Ratio[iname][ncore]) + ' ' + str(self.Speed[iname][ncore])
            string = string + ' %5.3e %5.3e %5.3e |'%(self.Mean[iname][ncore], self.Ratio[iname][ncore], self.Speed[iname][ncore])
         print string
      print ' '
      


   def PrintWiki(self):

      nCores = self.nCores
      Cores  = self.Cores

      print ' '
      print '====  %12s  ===='%(self.compName)
      if self.compName == 'Run':
         print '- Dynamics, Physcis: 통신, IO, 동기화 과정이 제외'
         print '- Communication: 역학의 모든 통신 시간의 합'
         print '- !Input/Output: 역학과 물리 과정의 모든 I/O 시간의 합'
         print '- Syncronize: 역학의 각 과정에서 발생하는 동기화 시간의 합'
      elif self.compName == 'Dynamics':
         print '- Syncronize: 역학-물리 접합에서 나타나는 동기화 시간'
         print '- 각 과정에서 발생하는 통신과 I/O 시간은 제외'
      elif self.compName == 'Physics':
         print '- 각 물리 컴포넌트에서 발생하는 I/O 시간은 제외'
      # line 1
      string = '||'*(3*self.nNames+2)
      print string
      # line 2
      string = '||=  =||'
      for iname in range(self.nNames):
         string = string + '||'*2 + '= ' + self.Names[iname] + ' =||'
      print string
      # line 3
      string='||= ncpus =||'+'= walltime [s] =||= Ratio [%] =||= Speed-up =||'*self.nNames
      print string
      # line 4
      string = '||'*(3*self.nNames+2)
      print string

      width_core = self.GetWidthNumber(Cores[-1])
      width_wall = self.GetWidthNumber(int(self.Mean[-1][0]))

      for ncore in range(nCores):
         #string = '|| %04d ||'%(Cores[ncore])
         string = '|| {0:0{width}d} ||'.format(Cores[ncore],width=width_core)
         for iname in range(self.nNames):
            #string = string + '  %5.3e  ||  %05.2f  ||  %07.1f  ||'%(self.Mean[iname][ncore], self.Ratio[iname][ncore], self.Speed[iname][ncore])
            string = string + '  {0:0{width_wall}.2f}  ||  {1:05.2f}  ||  {2:0{width_core}.1f}  ||'.format(self.Mean[iname][ncore], self.Ratio[iname][ncore], self.Speed[iname][ncore], width_wall=width_wall+3, width_core=width_core+3)
         print string

      # line final
      string = '||'*(3*self.nNames+2)
      print string
      print ' '
      if self.compName == 'Run': print '[wiki:KIM_performance_gaon2_model Figs]'
      if self.compName == 'Dynamics': print '[wiki:KIM_performance_gaon2_Dynamics Figs]'
      print '[[br]]'
      print ' '
      


   def PrintWiki_Diff(self):

      nCores = self.nCoresD
      Cores  = self.CoresD

      width  = self.GetWidthNumber(Cores[-1])

      print ' '
      print '====  %12s  ===='%(self.compName)
      # line 1
      string = '||'*(3*self.nNames+2)
      print string
      # line 2
      string = '||=  =||'
      for iname in range(self.nNames):
         string = string + '||'*2 + '= ' + self.Names[iname] + ' =||'
      print string
      # line 3
      string='||= ncpus =||'+'= walltime [s] =||= Ratio [%] =||= Speed-up =||'*self.nNames
      print string
      # line 4
      string = '||'*(3*self.nNames+2)
      print string

      width_core = self.GetWidthNumber(Cores[-1])
      width_wall = self.GetWidthNumber(int(self.Mean[-1][0]))

      stringd = ''
      for ncore in range(nCores):
         icore = self.mapD[ncore]
         string  = '||  {0:0{width}d}  ||'.format(self.Cores[icore],width=width_core)
         stringd = '||  {0:0{width}d}  ||'.format(Cores[ncore],width=width_core)
         #string  = '|| %04d ||'%(self.Cores[icore])
         #stringd = '|| %04d ||'%(Cores[ncore])
         stringr = '||  Speed-up  ||'
         for iname in range(self.nNames):
#            string  = string + '  %5.3e  ||  %05.2f  ||  %07.1f  ||'%(self.Mean[iname][icore], self.Ratio[iname][icore], self.Speed[iname][icore])
#            stringd = stringd+ '  %5.3e  ||  %05.2f  ||  %07.1f  ||'%(self.MeanD[iname][ncore], self.RatioD[iname][ncore], self.SpeedD[iname][ncore])
            string  = string + '  {0:0{width_wall}.2f}  ||  {1:05.2f}  ||  {2:0{width_core}.1f}  ||'.format(self.Mean[iname][icore], self.Ratio[iname][icore], self.Speed[iname][icore], width_wall=width_wall+3, width_core=width_core+3)
            stringd = stringd + '  {0:0{width_wall}.2f}  ||  {1:05.2f}  ||  {2:0{width_core}.1f}  ||'.format(self.MeanD[iname][ncore], self.RatioD[iname][ncore], self.SpeedD[iname][ncore], width_wall=width_wall+3, width_core=width_core+3)
            #if self.Mean[iname][icore] != 0.0:
            if self.MeanD[iname][ncore] != 0.0:
               #result = (self.Mean[iname][icore]-self.MeanD[iname][ncore])/self.Mean[iname][icore]*100.0
               result = self.Mean[iname][icore]/self.MeanD[iname][ncore]
            else:
               result = 0.0
            stringr = stringr+ '  %5.2f %1s  ||  -  ||  -  ||'%(result, ' ')
         print string
         print stringd
         print stringr
 

      # line final
      string = '||'*(3*self.nNames+2)
      print string
      print '[[br]]'
      print ' '


   def GetWidthNumber(self, number):
      inum = 0
      nmax = 7
      for i in range(nmax):
         num = 10**(i+1)
         if number/num == 0:
            inum = i+1
            break
      return inum



   def WriteFile(self, filedir=''):
      import os
      if filedir == '':
         filedir = './data'

      if not os.path.isdir(filedir):
         os.system('mkdir '+filedir)

      print '## Write statistics ... :', self.compName
      # self.Mean, self.Ratio, self. Speed
      outfilename = filedir+'/'+self.compName+'.ascii'
      outfile = open(outfilename, 'w')
      string = ''
      for ncore in range(self.nCores):
         string = string + ' %04d |'%(self.Cores[ncore])
         for iname in range(self.nNames):
            string = string + ' %f %f %f '%(self.Mean[iname][ncore], self.Ratio[iname][ncore], self.Speed[iname][ncore])
         string = string+'\n'
      outfile.write(string)
      outfile.close()



   def ReadFile(self, filedir=None):
      import os
      if filedir == None:
         filedir = '/home/jhkim/Study/Library/Python/ProfileKIM/data'

      infilename = filedir+'/'+self.compName+'.ascii'
      if not os.path.isfile(infilename):
         print 'ERROR :: file not found... :', infilename
         quit()

      print '## Read  statistics ... :', self.compName
      infile = open(infilename, 'r')
      string = infile.read()
      s_string = string.split('\n')

      nCores = len(s_string)-1
      Cores = []
      for ncore in range(nCores):
         ss_string = s_string[ncore].split()
         Cores.append(int(ss_string[0]))
      self.IniCores(Cores)

      for ncore in range(self.nCores):
         ss_string = s_string[ncore].split()
         for iname in range(self.nNames):
            self.Mean[iname][ncore]  = float(ss_string[iname*3+2])
            self.Ratio[iname][ncore] = float(ss_string[iname*3+3])
            self.Speed[iname][ncore] = float(ss_string[iname*3+4])
  
      infile.close()
