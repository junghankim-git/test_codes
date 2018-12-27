import os
import sys
import time
import filecmp

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


class diff_files:
   def __init__(self, sourceDir='.', targetDir='.', extension=[], exclude=[], onQuiet=False, onPrtPattern=True, isMoveFile=False):
      if sourceDir[-1] == '/': sourceDir = sourceDir[0:-1]
      if targetDir[-1] == '/': targetDir = targetDir[0:-1]
      if not isinstance(extension, list):
         print 'extension must be list...'
         quit()

      self.width        = 77
      self.sourceDir    = sourceDir
      self.targetDir    = targetDir
      self.extension    = extension
      self.exclude      = exclude
      self.findWord     = ''
      self.onQuiet      = onQuiet      # Print file name
      self.onPrtPattern = onPrtPattern # Print pattern word
      self.isMoveFile   = isMoveFile   #

      self.source = ''
      self.target = ''


   def run_process(self, args):
      print '='*self.width
      self.PrintDirectory(self.sourceDir, self.targetDir)
      self.ListDirectory(self.sourceDir, self.targetDir)
      print '='*self.width


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


# START here
   def ListDirectory(self, src_dir='./', trg_dir='./'):
      if src_dir[-1] == '/': src_dir = src_dir[0:-1]
      if trg_dir[-1] == '/': trg_dir = trg_dir[0:-1]
      self.PrintDetailDir(src_dir, trg_dir, True)
      n_src_dirs = 0
      src_dirs  = []
      n_trg_dirs = 0
      trg_dirs  = []
      for ifile in os.listdir(src_dir):
         if ifile == '.' or ifile =='..': continue
         if self.CheckException(ifile): continue
         src_fname = src_dir+'/'+ifile
         trg_fname = trg_dir+'/'+ifile
         if os.path.isdir(src_fname):
            src_dirs.append(src_fname)
            n_src_dirs = n_src_dirs + 1
            trg_dirs.append(trg_fname)
            n_trg_dirs = n_trg_dirs + 1
            if os.path.isdir(trg_fname):
               self.PrintState(ifile, ifile, 'same')
            else:
               self.PrintState(ifile, ifile, 'source')
         else:
            self.PrintDetailFile(src_fname,True)
            if self.CheckExtension(src_fname):
               if not os.path.exists(trg_fname):
                  self.PrintState(ifile, ifile, 'source')
               else:
                  if filecmp.cmp(src_fname, trg_fname):
                     self.PrintState(ifile, ifile, 'same')
                  else:
                     self.PrintState(ifile, ifile, 'diff')
               #if self.IsASCII(src_fname):
               #else:
               #   self.PrintState(src_fname, trg_fname, 'not readable')
            self.PrintDetailFile(src_fname,False)
      self.PrintDetailDir(src_dir, trg_dir, False)

      for i in range(n_src_dirs):
         print '='*self.width
         self.PrintDirectory(src_dirs[i], trg_dirs[i])
         self.ListDirectory(src_dirs[i], trg_dirs[i])


   def IsSame(self, src, trg):
      return True


   def PrintDetailDir(self, dirname, dirname2, isStart=True):
      if not self.onQuiet:
         if isStart:
           print '='*80
           print '+ START - Directory :: {0:80s}'.format(dirname)
           print '+                   :: {0:80s}'.format(dirname2)
         else:
           print '+ END   - Directory :: {0:80s}'.format(dirname)
           print '                    :: {0:80s}'.format(dirname2)
           print '='*80
           print '\n\n'


   def PrintDetailFile(self, filename, isStart=True):
      if not self.onQuiet:
         if isStart:
           print '-'*80
           print '+ START - File :: {0:80s}'.format(filename)
         else:
           print '+ END   - File :: {0:80s}'.format(filename)
           print '-'*80


   def PrintDirectory(self, src_dir, trg_dir):
      #print '| source : {0:107} |'.format(src_dir)
      #print '| target : {0:107} |'.format(trg_dir)
      #print '='*120
      print '| source : {0:84} |'.format(src_dir)
      print '| target : {0:84} |'.format(trg_dir)
      print '='*self.width


   def PrintState(self, src, trg, value):
      #print '| {0:50s} | {1:50s} | {2:10s} |'.format(src, trg, value)
      #print '-'*120
      print '| {0:30s} | {1:30s} | {2:7s} |'.format(src, trg, value)
      print '-'*self.width
