#!/usr/bin/env python

import os
import sys
import re
from optparse import OptionParser
from f90convert import *

class F90Grammer:

    def __init__(self, filename):
        self.ntypes       = 7
        self.types_info   = [[[False,-1,-1] for j in range(self.ntypes)] for i in range(self.ntypes)]
        self.name         = []
        self.type         = []
        self.need         = []
        self.goto         = []
        self.end          = []
        self.rename       = []
        self.read_grammer_file(filename)
        self.n            = len(self.name)
        self.reset()
        self.debug        = False
        self.typ_idx      = 6
  
  
    def reset(self):
        self.types     = [0]
        self.intents   = [0]
        self.ctype     = 0
        self.cintent_h = 0
        self.is_cont   = False
  
  
    def read_grammer_file(self, filename):
  
        # intents (new)
        infile = open(filename, 'r')

        for line in infile:
            sline = line.split()

            if len(sline) == 6 and sline[0] == '@DEF':
                otype = int(sline[1])
                ntype = int(sline[2])

                if int(sline[3]) == 0:
                    depend = False
                else:
                    depend = True

                intenth = int(sline[4])
                intentv = int(sline[5])
                self.types_info[otype][ntype][0] = depend
                self.types_info[otype][ntype][1] = intenth
                self.types_info[otype][ntype][2] = intentv

        infile.close()

        infile = open(filename, 'r')

        for line in infile:

            if line.strip() != '' and line.strip()[0] != '@':
                sline = (line.strip()).split()

                if len(sline) != 6:
                    print 'check grammer file...'
                    print line.strip()
                    quit()
                else:
                    self.name.append(sline[0])
                    self.type.append(int(sline[1]))
                    self.need.append(sline[2])
                    self.goto.append(int(sline[3]))
                    self.end.append(int(sline[4]))
                    self.rename.append(sline[5])
    
        infile.close()
  
  
  
  
    def find_index(self, string):
        ustring = string.upper()
        idx = -1

        for i in range(self.n):
            if ustring == self.name[i]:
                idx = i
                break
            elif self.name[i] == '1234' and ustring.isdigit():
                idx = i
                break
    
        return idx
  
  
  
  
    def replace_line_before(self, line):
        # replace
        applied_line    = line
        applied_line    = applied_line.replace('(', ' ( ')
        applied_line    = applied_line.replace(')', ' ) ')
        applied_line    = applied_line.replace('*', ' * ')
        applied_line    = applied_line.replace('+', ' + ')
        applied_line    = applied_line.replace('-', ' - ')
        applied_line    = applied_line.replace('/', ' / ')
        applied_line    = applied_line.replace('**', ' ** ')
        applied_line    = applied_line.replace('=', ' = ')
        applied_line    = applied_line.replace(':', ' : ')
        applied_line    = applied_line.replace(',', ' , ')
#    applied_line    = applied_line.replace('.EQ.', ' .EQ. ')
#    applied_line    = applied_line.replace('.NE.', ' .NE. ')
#    applied_line    = applied_line.replace('.GT.', ' .GT. ')
#    applied_line    = applied_line.replace('.GE.', ' .GE. ')
#    applied_line    = applied_line.replace('.LT.', ' .LT. ')
#    applied_line    = applied_line.replace('.LE.', ' .LE. ')
        applied_line    = applied_line.replace('==', ' == ')
        applied_line    = applied_line.replace('!=', ' != ')
        applied_line    = applied_line.replace('>', ' > ')
        applied_line    = applied_line.replace('>=', ' >= ')
        applied_line    = applied_line.replace('<', ' < ')
        applied_line    = applied_line.replace('<=', ' <= ')
#    applied_line    = applied_line.replace('.NOT.', ' .NOT. ')
#    applied_line    = applied_line.replace('.NOT.', ' .NOT. ')
        return applied_line



    def replace_line_after(self, line):
        applied_line = line
        applied_line = applied_line.replace('= =', '==')
        applied_line = applied_line.replace('= >', '=>')
        applied_line = applied_line.replace('< =', '<=')
        applied_line = applied_line.replace(': :', '::')
        applied_line = applied_line.replace(' (', '(')
        applied_line = applied_line.replace('( ', '(')
        applied_line = applied_line.replace(' )', ')')
        applied_line = applied_line.replace(') %', ')%')
        applied_line = applied_line.replace(' ,', ',')
        # add
        applied_line = applied_line.replace(' => ', '=>')
        applied_line = applied_line.replace('/ /', '//')
        applied_line = applied_line.replace(' // ', '//')
        applied_line = applied_line.replace('len = ', 'len=')
        applied_line = applied_line.replace('kind = ', '')
        applied_line = applied_line.replace('kind=', '')
        applied_line = applied_line.replace(' .and. ', '.and.')
        applied_line = applied_line.replace(' .and.', '.and.')
        applied_line = applied_line.replace('.and. ', '.and.')
        applied_line = applied_line.replace(' .or. ', '.or.')
        applied_line = applied_line.replace(' .or.', '.or.')
        applied_line = applied_line.replace('.or. ', '.or.')
        applied_line = applied_line.replace(' .not. ', '.not.')
        applied_line = applied_line.replace(' .not.', '.not.')
        applied_line = applied_line.replace('.not. ', '.not.')
    
        applied_line = applied_line.replace(' != ', '!=')
        applied_line = applied_line.replace(' !=', '!=')
        applied_line = applied_line.replace('!= ', '!=')
        applied_line = applied_line.replace(' .ne. ', '.ne.')
        applied_line = applied_line.replace(' .ne.', '.ne.')
        applied_line = applied_line.replace('.ne. ', '.ne.')
#
        applied_line = applied_line.replace(' == ', '==')
        applied_line = applied_line.replace(' ==', '==')
        applied_line = applied_line.replace('== ', '==')
        applied_line = applied_line.replace(' .eq. ', '.eq.')
        applied_line = applied_line.replace(' .eq.', '.eq.')
        applied_line = applied_line.replace('.eq. ', '.eq.')
        applied_line = applied_line.replace(' > ', '>')
        applied_line = applied_line.replace(' >', '>')
        applied_line = applied_line.replace('> ', '>')
        applied_line = applied_line.replace(' .gt. ', '.gt.')
        applied_line = applied_line.replace(' .gt.', '.gt.')
        applied_line = applied_line.replace('.gt. ', '.gt.')
        applied_line = applied_line.replace(' >= ', '>')
        applied_line = applied_line.replace(' .ge. ', '.ge.')
        applied_line = applied_line.replace(' .ge.', '.ge.')
        applied_line = applied_line.replace('.ge. ', '.ge.')
        applied_line = applied_line.replace(' < ', '<')
        applied_line = applied_line.replace(' <', '<')
        applied_line = applied_line.replace('< ', '<')
        applied_line = applied_line.replace(' .lt. ', '.lt.')
        applied_line = applied_line.replace(' .lt.', '.lt.')
        applied_line = applied_line.replace('.lt. ', '.lt.')
        applied_line = applied_line.replace(' <= ', '<=')
        applied_line = applied_line.replace(' <=', '<=')
        applied_line = applied_line.replace('<= ', '<=')
        applied_line = applied_line.replace(' .le. ', '.le.')
        applied_line = applied_line.replace(' .le.', '.le.')
        applied_line = applied_line.replace('.le. ', '.le.')
        applied_line = applied_line.replace('un = ', 'un=')
        applied_line = applied_line.replace('len = ', 'len=')
        applied_line = applied_line.replace('intent(in)', 'intent(in   )')
        applied_line = applied_line.replace('intent(out)', 'intent(  out)')
        applied_line = applied_line.replace('dim = ', 'dim=')
        applied_line = applied_line.replace(':, :', ':,:')
        applied_line = applied_line.replace('end do', 'enddo')
        applied_line = applied_line.replace('else if', 'elseif')
        applied_line = applied_line.replace('end if', 'endif')
        applied_line = applied_line.replace('if(', 'if (')
        applied_line = applied_line.replace('elseif(', 'elseif (')
        applied_line = applied_line.replace(' + ', '+')
        applied_line = applied_line.replace(' +', '+')
        applied_line = applied_line.replace('+ ', '+')
        applied_line = applied_line.replace(' - ', '-')
        applied_line = applied_line.replace(' -', '-')
        applied_line = applied_line.replace('- ', '-')
        applied_line = applied_line.replace(' * ', '*')
        applied_line = applied_line.replace(' *', '*')
        applied_line = applied_line.replace('* ', '*')
        applied_line = applied_line.replace(' / ', '/')
        applied_line = applied_line.replace(' /', '/')
        applied_line = applied_line.replace('/ ', '/')
        applied_line = applied_line.replace('intent', '                      intent')
#
        applied_line = re.sub('^%s'%'   end subroutine', '!\n   return\n   end subroutine', applied_line)
        applied_line = re.sub('^%s'%'   end function', '!\n   return\n   end function', applied_line)
        applied_line = re.sub('^%s'%'   implicit none','!-------------------------------------------------------------------------------\n   implicit none\n!',applied_line)
        applied_line = re.sub('^%s'%'   subroutine','!-------------------------------------------------------------------------------\n!\n!\n!-------------------------------------------------------------------------------\n   subroutine',applied_line)
        applied_line = re.sub('^%s'%'   function','!-------------------------------------------------------------------------------\n!\n!\n!-------------------------------------------------------------------------------\n   function',applied_line)
#        applied_line = re.sub('^%s'%'   module', '!-------------------------------------------------------------------------------\n   module', applied_line)
        applied_line = re.sub('^%s'%'   end module', '!-------------------------------------------------------------------------------\n!\n!\n!-------------------------------------------------------------------------------\n   end module', applied_line)
        #applied_line = re.sub('^%s'%' ', ' ', applied_line)
    
        return applied_line



    def get_type_indexes(self, upper_line, splict_upper_line):
        indexes  = []
        types    = []
        needs    = []
        gotos    = []
        ends     = []
        renames  = []
        nn       = len(splict_upper_line)

        for i in range(nn):
            ig = self.find_index(splict_upper_line[i])
            #if ig >= 0:
            indexes.append(i)
            types.append(self.type[ig])
            needs.append(self.need[ig])
            gotos.append(self.goto[ig])
            ends.append(self.end[ig])
            renames.append(self.rename[ig])

        # type, isend
        nmatched = len(types)

        if self.is_cont:
            type  = -2
            isend = False
            self.is_cont = False
        elif nmatched == 0 or indexes[0] != 0:
            type  = 5
            isend = False
        else:
            if needs[0] == 'none' or upper_line.count(needs[0]) > 0:
                type = types[0]
            else:
                type = gotos[0]

            if ends[0] != 0:
                isend = True
            else:
                isend = False
    
        return type, types, isend, indexes
  
  
  
    def end_type(self, linetype, line):
        # Pop
        if linetype == -1:

            if self.types.count(self.typ_idx) > 0:

                if self.ctype == self.typ_idx:
                    before_intent = self.intents[-1]
                    tmp = self.types.pop()
                    tmp = self.intents.pop()
                else:
                    before_intent = self.intents[-2]
                    tmp = self.types.pop()
                    tmp = self.types.pop()
                    tmp = self.intents.pop()
                    tmp = self.intents.pop()

                if self.types[-1] == self.typ_idx:
                    self.intents[-1] = before_intent

            else:

                while self.types[-1] != 0 and self.types[-1] != 1 and self.types[-1] != 2:
                    tmp = self.types.pop()
                    tmp = self.intents.pop()

                if line.strip() != 'CONTAINS' and line.strip() != 'contains':
                  #if line.strip() != 'CONTAINS':
                    before_intent = self.intents[-1]
                    tmp = self.types.pop()
                    tmp = self.intents.pop()
                    if self.types[-1] != 0: self.intents[-1] = before_intent

        else:   #if linetype != -1:
            tmp = self.types.pop()
            tmp = self.intents.pop()
            self.types[-1] = linetype

        self.ctype     = self.types[-1]
        self.cintent_h = self.intents[-1]
        applied_line   = ' '*self.cintent_h + line
    
        return applied_line
  

  
    def preprocess(self, line):
    
        if line.lstrip()[0] == '!' or line.lstrip()[0] == '#':
            type         = -3
            isend        = False
            applied_line = line.strip() + '\n'
        else:
            nspc    = 3 # leading space
            line    = self.replace_line_before(line)

            if self.is_cont: # continuity
                line = ' '+line.lstrip()

            s_line  = line.split()
            uline   = line.upper()
            s_uline = uline.split()
            lline   = line.lower()
            s_lline = lline.split()
            nstrs   = len(s_uline)
            type, types, isend, indexes = self.get_type_indexes(uline, s_uline)
            applied_line = ''

            for i in range(nstrs):
                #if types[i] != 1 and types[i] != 2 and types[i] != 3 and types[i] != 4 and types[i] != 6:
                #  applied_line = applied_line + s_lline[i].lstrip() + ' '
                #else:
                #  applied_line = applied_line + s_lline[i].lstrip() + ' '
                if s_line[i][0] == "'":
                  applied_line = applied_line + s_line[i].lstrip() + ' '
                else:
                  applied_line = applied_line + s_lline[i].lstrip() + ' '

            # remove continuation
            if applied_line.count('&') > 0:
                print 'is cont!!!(0)'
                applied_line = applied_line.replace('&','').rstrip()

                if self.is_cont:
                    applied_line = self.replace_line_after(applied_line)
                    print 'is cont!!!(1)'
                else:
                    applied_line = self.replace_line_after(' '*nspc+applied_line)

                self.is_cont = True
            else:
                applied_line = applied_line.rstrip()+'\n'
                applied_line = self.replace_line_after(' '*nspc+applied_line)
    
        return type, isend, applied_line



    def set_horizontal_intent(self, ntype, isend, line):
        otype = self.otype

        if   ntype == -3:
            applied_line = line
        elif ntype == -2:
            applied_line = line
        elif ntype == -1 or isend:
            applied_line = self.end_type(ntype, line)
        else:
            if self.types_info[otype][ntype][1] == -1:
                print 'could not change types... ', otype, ntype
                quit()

            self.ctype = ntype

            if self.types_info[otype][ntype][0]:
                self.cintent_h = self.cintent_h + self.types_info[otype][ntype][1]
            else:
                self.cintent_h = self.types_info[otype][ntype][1]

            if ntype != otype or ntype == self.typ_idx:
                self.types.append(self.ctype)
                self.intents.append(self.cintent_h)

            applied_line = ' '*self.cintent_h + line
    
        # for debug
        if self.debug:
            print 'line    = ', line.strip()
            print 'type    = ', ntype
            print 'types   = ', self.types
            print 'intents = ', self.intents
            print 'ctype   = ', self.ctype
            print 'cintent = ', self.cintent_h
            print ' '
            applied_line = ':%02d:%02d:'%(ntype, self.ctype) + applied_line
    
        return applied_line


 
    def set_vertical_intent(self, line):
  
      otype    = self.otype
      ntype    = self.ctype
      intent_v = 0

      if   ntype == -3:
        pass
      elif ntype == -2:
        pass
      else:
        intent_v = self.types_info[otype][ntype][2]
  
      applied_line = '\n'*intent_v + line
  
      return applied_line



    def generate_line(self, line):

        if line.strip() == '':
            applied_line = None
        else:
            self.otype = self.ctype
            type, isend, applied_line = self.preprocess(line)
            #applied_line              = self.set_horizontal_intent(type, isend, applied_line)
            #applied_line              = self.set_vertical_intent(applied_line)

        return applied_line
  



class F90Converter:



    def __init__(self, infilepath, oufilepath, extension=[], exclude=[]):
        self.infilepath = infilepath
        self.oufilepath = oufilepath
        self.extension  = extension
        self.exclude    = exclude
        self.grammer    = F90Grammer('/home/jhkim/work/python/fortran/fortran.dat')



    def driver(self):
        infiles = os.listdir(self.infilepath)

        for ifile in infiles:
            infilename = self.infilepath+'/'+ifile
            oufilename = self.oufilepath+'/'+ifile

            if not os.path.isdir(infilename):
                print 'Processing file : ', infilename
                self.convert(infilename, oufilename)
            else:
                print infilename, ' is directory'



    def convert(self, infilename, oufilename):
        infile = open(infilename, 'r')
        oufile = open(oufilename, 'w')
        oldline = ''
        iscont  = False

        for line in infile:
            applied_line = self.grammer.generate_line(line)
            if applied_line != None:
                oufile.write(applied_line)
            else:
                oufile.write('\n')

        infile.close()
        oufile.close()
        self.grammer.reset()




def main():

    opt = OptionParser()
    # action: 'store', 'store_const', 'append', 'count', 'callback'
    opt.add_option('-s', '--source',      dest='src',  default='./src', action='store',  help='source')
    opt.add_option('-d', '--destination', dest='dst',  default='./dst', action='store',  help='destination')
    (options, args) = opt.parse_args()
    
    converter = F90Converter(options.src, options.dst)
    converter.driver()



if __name__ == '__main__':
    main()

    print(" 1) '-------'")
    print(" 2) ': '")
    print(" 3) 'intent'")
    print(" 4) ', '")




