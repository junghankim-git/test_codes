import os
import sys
import time
sys.path.append('/home/jhkim/work/share/python')
from logger import INFO, WARN, FATAL, logger

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


class class_analysis_files:
    def __init__(self, target_path='.', extension=[], exclude=[], on_quiet=True, \
                 on_print_ptn=True, on_replace_ptn=False, on_move_file=False):
        if target_path[-1]=='/': target_path = target_path[0:-1]
        if not isinstance(extension,list):
           logger('extension must be list...',FATAL,'__init__',self)
  
        self.target_path    = target_path
        self.extension      = extension
        self.exclude        = exclude
        self.on_quiet       = on_quiet        # print file name
        self.on_print_ptn   = on_print_ptn    # print pattern word
        self.on_replace_ptn = on_replace_ptn  # replace pattern word
        self.on_move_file   = on_move_file    # move file name
        self.pattern        = ''
        self.replace        = ''
 
 
    def run_process(self, args):
        if not isinstance(args,list):
            logger('arguments must be list...',FATAL,'run_process',self)
        if self.on_move_file:
            if len(args)!=2: logger('# of arguments must be 2 ...',FATAL,'run_process - 3',self)
            self.on_print_ptn = False
            self.on_replace_ptn = False
        if self.on_replace_ptn:
            if len(args)!=2: logger('# of arguments must be 2 ...',FATAL,'run_process - 2',self)
            self.on_print_ptn
        if self.on_print_ptn:
            if len(args)<1: logger('# of arguments must greater than 0 ...',FATAL,'run_process - 1',self)
  
        if len(args)>0: self.pattern = args[0]
        if len(args)==2: self.replace = args[1]
        
        print(' ')
        print(' # Processing: START')
        print(' ')
        self.process_directory(self.target_path)
        print(' ')
        print(' # Processing: END')
 
 
    def check_extension(self, filename):
        next = len(self.extension)
        if next==0: return True

        isext = False
        for iext in self.extension:
            ne = len(iext)
            if filename[-ne:]==iext:
                isext = True
                exit
        return isext
 
 
    def check_exception(self, filename):
        next = len(self.exclude)
        if next==0: return False
        for iext in self.exclude:
            ne = len(iext)
            if filename[-ne:]==iext:
                return True
            else:
                return False
 
 
    def check_ascii(self, filename):
        command  = 'file '+filename
        process  = os.popen(command)
        att_file = process.read()
        how1 = att_file.find('ASCII')
        how2 = att_file.find('script')
        how3 = att_file.find('text')
        how4 = att_file.find('FORTRAN')
        if how1<0 and how2<0 and how3<0 and how4<0:
            return False
        else:
            return True
 
 
 
    def process_pattern(self, filename):
        file = open(filename,'r')
        if self.on_replace_ptn:
            outfilename = filename+'.tmp'
            outfile = open(outfilename,'w')
        found = False
        iline = 0
        for line in file:
            iline = iline+1
            if line.find(self.pattern)>=0:
                found = True
                print(' => line(%5d) :: %s'%(iline,line.replace('\n','')))
                if self.on_replace_ptn: outfile.write(line.replace(self.pattern,self.replace))
            else:
                if self.on_replace_ptn: outfile.write(line)
        if found: self.print_found(filename)
        file.close()
        if self.on_replace_ptn:
            if found:
                outfile.close()
                os.system('mv '+outfilename+' '+filename)
            else:
                os.system('rm -f '+outfilename)




    def print_found(self, filename):
            print('%-15s matched in file: %-31s %15s '%('='*15,filename,'='*15))
            print(' ')
 

 
    def print_path(self, dirname, is_start=True):
        if is_start:
            print('='*80)
            print('+ start - path: {0:80s}'.format(dirname))
        else:
            print('+ end   - path: {0:80s}'.format(dirname))
            print('='*80)
            print('\n\n')
 
 

    def print_file(self, filename, is_start=True):
        if is_start:
            print('-'*80)
            print('+ start - file: {0:80s}'.format(filename))
        else:
            print('+ end   - file: {0:80s}'.format(filename))
            print('-'*80)
            print(' ')
 
 
 

    def move_file(self, dir, file):
        if file==self.pattern:
            infilename  = dir+'/'+file
            outfilename = dir+'/'+self.replace
            command = 'mv '+infilename+' '+outfilename
            os.system(command)


 
 
    def process_directory(self, dir='./'):
        if dir[-1]=='/': dir = dir[0:-1]
        if not self.on_quiet: self.print_path(dir,True)
        ndirs = 0
        dirs  = []
        for ifile in os.listdir(dir):
            filename = dir+'/'+ifile
            if ifile=='.' or ifile=='..': continue
            if self.check_exception(ifile): continue
            if os.path.isdir(filename):
                dirs.append(filename)
                ndirs = ndirs+1
            else:
                if not self.on_quiet: self.print_file(filename,True)
                #
                if self.check_extension(filename) and self.check_ascii(filename) and self.on_print_ptn:
                    self.process_pattern(filename)
                #
                if not self.on_quiet: self.print_file(filename,False)
        if not self.on_quiet: self.print_path(dir,False)
  
        for i in range(ndirs):
             self.process_directory(dirs[i])


