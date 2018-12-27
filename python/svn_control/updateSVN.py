#!/usr/bin/env python
#-------------------------------------------------------------------------------
#
#> @brief
#>  - 
#>
#> @date 11MAR2016
#>  - JH KIM(jh.kim@kiaps.org) : First written
import os
import sys
from optparse import OptionParser
#sys.path.append('/home/jhkim/Study/Library/Shared/Python')


########################################################################################################################
# Classes and Functions
########################################################################################################################

class KIM_svn:

  def __init__(self, updateFile):
    self.checked    = False
    self.KIMDIR     = ''
    self.SRCDIR     = ''
    self.updateFile = updateFile

    if not os.path.isfile(updateFile):
      print ' [ err] updateFile is not file...', updateFile
      quit()



  def CheckUpdateFile(self):
    # Load update file
    ufile = open(self.updateFile, 'r')
    self.KIMDIR = self.GetDirName(ufile, 'KIMDIR')
    self.SRCDIR = self.GetDirName(ufile, 'SRCDIR')
    print '[ log] KIMDIR = ', self.KIMDIR
    print '[ log] SRCDIR = ', self.SRCDIR
    self.mlines, self.filelist = self.GetFilelist(ufile)
    print '[ log] total line = ', self.mlines
    #print '[ log] filelist   = ', self.filelist
    ufile.close()

    # Check update contents
    isready = self.CheckFiles(self.filelist)
    if isready:
      print '[ log] ready for update...'
    else:
      print '[ err] check update list...'
      quit()



  def GetDirName(self, infile, name):
    infile.seek(0, 0)
    for line in infile:
      uline = line.strip()
      start = uline.find(name)
      if start >= 0 and uline[start-1] != '$':
        eqpos = uline.find('=')
        if eqpos > start:
          dirname = uline[eqpos+1:]
          if os.path.isdir(dirname):
            return dirname
            break
          else:
            print '[ err] '+dirname+' is not directory...'
            quit()
        else:
          print '[ err] check dirname= ...'
          quit()



  def GetFilelist(self, infile):
    infile.seek(0, 0)
    mlines = [0 for i in range(4)]
    filelist = [['' for i in range(2)] for j in range(200)]
    for line in infile:
      uline = line.strip()
      startM = uline[0:4].find('(M)')
      startD = uline[0:4].find('(D)')
      startA = uline[0:4].find('(A)')
      if startM >= 0 or startD >= 0 or startA >= 0:
        ikim = uline.find('$KIMDIR/')
        if ikim >= 0:
          if startM >= 0:
            mlines[0] = mlines[0] + 1 
            filelist[mlines[3]][0] = '(M)'
          if startD >= 0:
            mlines[1] = mlines[1] + 1 
            filelist[mlines[3]][0] = '(D)'
          if startA >= 0:
            mlines[2] = mlines[2] + 1 
            filelist[mlines[3]][0] = '(A)'
          filelist[mlines[3]][1] = uline[ikim+8:].strip()
          mlines[3] = mlines[3] + 1
        else:
          print '[ err] (M) or (D) or (A) must has $KIMDIR...'
          quit()
    return mlines, filelist[:][0:mlines[3]]



  def CheckFiles(self, filelist):
    isready = True
    nfiles = len(filelist)

    print ' '
    print '=============================================================================================='
    print '| KIMDIR = %-81s |'%(self.KIMDIR)
    print '| SRCDIR = %-81s |'%(self.SRCDIR)
    print '=============================================================================================='
    print '|                                     filename                                         | tag |'
    print '=============================================================================================='
    for ifile in range(nfiles):
      kimfile   = self.KIMDIR+'/'+filelist[ifile][1]
      srcfile   = self.SRCDIR+'/'+filelist[ifile][1]
      iskimfile = os.path.isfile(kimfile)
      issrcfile = os.path.isfile(srcfile)
      iskimdir  = os.path.isdir(kimfile)
      issrcdir  = os.path.isdir(srcfile)
      print '| %-84s | %3s |'%(filelist[ifile][1], filelist[ifile][0])
      if filelist[ifile][0] == '(M)':
        if not iskimfile or not issrcfile:
          isready = False
          print '[ err] (M) check file..', kimfile, srcfile
      elif filelist[ifile][0] == '(D)':
        if not iskimfile and not iskimdir:
          isready = False
          print '[ err] (D) check file..', kimfile
      elif filelist[ifile][0] == '(A)':
        if iskimfile or iskimdir or (not issrcfile and not issrcdir):
          isready = False
          print '[ err] (A) check file..', kimfile, srcfile
      else:
        isready = False
        print '[ err] check tag (M, D, A)....', filelist[ifile][0]
    print '=============================================================================================='
    print ' '



  def UpdateSVN(self):
    nfiles = len(self.filelist)

    print ' '
    print '=============================================================================================='
    print '| KIMDIR = %-81s |'%(self.KIMDIR)
    print '| SRCDIR = %-81s |'%(self.SRCDIR)
    print '=============================================================================================='
    print '|                                     filename                                         | tag |'
    print '=============================================================================================='
    for ifile in range(nfiles):
      kimfile = self.KIMDIR+'/'+self.filelist[ifile][1]
      srcfile = self.SRCDIR+'/'+self.filelist[ifile][1]
      print '| %-84s | %3s |'%(self.filelist[ifile][1], self.filelist[ifile][0])
      if self.filelist[ifile][0] == '(M)':
        command = 'cp '+srcfile+' '+kimfile
        os.system(command)
      elif self.filelist[ifile][0] == '(D)':
        command = 'svn rm '+kimfile
        os.system(command)
      elif self.filelist[ifile][0] == '(A)':
        command = 'cp -rf '+srcfile+' '+kimfile
        os.system(command)
        command = 'svn add '+kimfile
        os.system(command)
    print '=============================================================================================='
    print ' '
       



########################################################################################################################
# Main
########################################################################################################################

opt = OptionParser()
# action: 'store', 'store_const', 'append', 'count', 'callback'
opt.add_option('-s', '--file',    dest='ufile',     default='./kim.txt',  action='store',       help='update list file')
#opt.add_option('-k', '--kimdir',  dest='kimdir',    default='none',  action='store',       help='KIMDIR path')
opt.add_option('-c', '--check',   dest='check',     default=True,  action='store_false', help='Print check')
opt.add_option('-u', '--update',  dest='update',    default=False, action='store_true',  help='svn update')

(options, args) = opt.parse_args()

#print 'Options = ', options
#print 'Args    = ', args

uKIM = KIM_svn(options.ufile)

if options.check:
  uKIM.CheckUpdateFile()

if options.update:
  uKIM.UpdateSVN()




