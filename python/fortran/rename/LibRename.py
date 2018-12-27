import os
import sys
import time
#SHARE_DIR=os.environ['PYTHON_SHARE_DIR']
SHARE_DIR='/home/jhkim/work/share/python'
sys.path.append(SHARE_DIR)
from LibErrorHandler import INFO, WARN, ERROR, Logger

#######################################################################################
# Description: 
#
# Date: -----2015
#   - JH KIM: First written
#
# Date: 09APR2015
#   - JH KIM: Functions for viewing difference between resources
#
#######################################################################################


class AnalysisFiles:
   def __init__(self, targetDir='.', extension=[], exclude=[], onQuiet=True, onPrtPattern=True, onRpcPattern=False, onMoveFile=False):
      if targetDir[-1] == '/': targetDir = targetDir[0:-1]
      if not isinstance(extension, list):
         Logger('extension must be list...', ERROR, '__init__', self)

      self.targetDir    = targetDir
      self.extension    = extension
      self.exclude    = exclude
      self.findWord     = ''
      self.onQuiet      = onQuiet      # Print file name
      self.onPrtPattern = onPrtPattern # Print pattern word
      self.onRpcPattern = onRpcPattern # Replace pattern word
      self.onMoveFile   = onMoveFile   # Move file name

      self.source = ''
      self.sedin = ''


   def RunProcess(self, sedinfile):

      if not os.path.isfile(sedinfile):
        print 'check sed in file...', sedinfile
        quit()
      self.sedin = sedinfile.strip()
      
      self.ListDirectory(self.targetDir)


   def CheckExtension(self, filename):
      next = len(self.extension)
      if next == 0: return True
      for iext in self.extension:
         ne = len(iext)
         if filename[-ne:] == iext:
            return True
      return False


   def CheckException(self, filename):
      next = len(self.exclude)
      if next == 0: return False
      for iext in self.exclude:
         ne = len(iext)
         if filename[-ne:] == iext:
            return True
      return False


   def IsASCII(self, filename):
      command  = 'file '+filename
      process  = os.popen(command)
      att_file = process.read()
      how1  = att_file.find('ASCII')
      how2  = att_file.find('script')
      how3  = att_file.find('text')
      if how1 < 0 and how2 < 0 and how3 < 0:
         return False
      else:
         return True



   def PrintPattern(self, filename):
      if not self.onPrtPattern: return 0
      file = open(filename, 'r')

      if self.onRpcPattern:
        outfilename = filename+'.tmp'
        outfile = open(outfilename, 'w')
      found = False
      iline = 0
      for line in file:
         iline = iline + 1
         if line.find(self.source) >= 0:
            found = True
            print ' => LINE: '+str(iline)+' :: '+line.replace('\n','')
            if self.onRpcPattern: outfile.write(line.replace(self.source, self.target))
         else:
            if self.onRpcPattern: outfile.write(line)
      if found:
         print '================  Matched in File: '+filename+' ================='
         print ' '
      file.close()
      if self.onRpcPattern:
         if found:
            outfile.close()
            os.system('mv '+outfilename+' '+filename)
         else:
            os.system('rm -f '+outfilename)


   def DoRename(self, filename):
      command = 'sed -f '+self.sedin+' '+filename+' > tmp.txt'
      print command
      os.system(command)
      command = 'mv tmp.txt '+filename
      print command
      os.system(command)



   def MoveFile(self, dir, file):
      if file == self.source:
         infilename  = dir+'/'+file
         outfilename = dir+'/'+self.target
         command = 'mv '+infilename+' '+outfilename
         os.system(command)


   def ListDirectory(self, dir='./'):
      if dir[-1] == '/': dir = dir[0:-1]
      self.PrintDir(dir, True)
      ndirs = 0
      dirs  = []
      for ifile in os.listdir(dir):
         filename = dir+'/'+ifile
         if ifile == '.' or ifile =='..': continue
         if self.CheckException(ifile): continue
         if os.path.isdir(filename):
            dirs.append(filename)
            ndirs = ndirs + 1
         else:
            self.PrintFile(filename,True)
            if self.CheckExtension(filename):
               self.MoveFile(dir, ifile)
               if self.IsASCII(filename):
                  #self.PrintPattern(filename)
                  self.DoRename(filename)
            self.PrintFile(filename,False)
      self.PrintDir(dir, False)

      for i in range(ndirs):
         self.ListDirectory(dirs[i])


   def PrintDir(self, dirname, isStart=True):
      if not self.onQuiet:
         if isStart:
           print '='*80
           print '+ START - Directory :: {0:80s}'.format(dirname)
         else:
           print '+ END   - Directory :: {0:80s}'.format(dirname)
           print '='*80
           print '\n\n'


   def PrintFile(self, filename, isStart=True):
      if not self.onQuiet:
         if isStart:
           print '-'*80
           print '+ START - File :: {0:80s}'.format(filename)
         else:
           print '+ END   - File :: {0:80s}'.format(filename)
           print '-'*80



