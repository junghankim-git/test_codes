#!/usr/bin/env python

import os
import sys
import re
from optparse import OptionParser
from f90convert import *




class F90Grammer:

    def __init__(self, filename):
        self.use_capital  = False
        self.merge_cont   = False
        self.base_hintent = 0
        self.base_vintent = 0
        self.vintent_mark = ''
        self.ntypes       = 7
        self.types_info   = [[[False,-1,-1] for j in range(self.ntypes)] for i in range(self.ntypes)]
        self.name         = []
        self.type         = []
        self.need         = []
        self.goto         = []
        self.end          = []
        self.rename       = []
        self.read_grammer_file(filename)
        self.n_chars      = len(self.name)
        self.reset()
        self.debug        = False
        # constants
        self.typ_idx      = 5
        self.max_width    = 80




    def reset(self):
        self.types     = [0]     # all   type stateus
        self.intents   = [0]     # all intent stateus
        self.ctype     = 0       # current type
        self.cintent_h = 0       # current horizontal intent
        self.is_cont   = False   # continuous state




    def read_grammer_file(self, filename):
  
        # intents (new)
        infile = open(filename, 'r')

        for line in infile:
            sline = line.split()

            if len(sline) == 2 and sline[0] == '@CAPITAL':
                self.use_capital = True

            if len(sline) == 2 and sline[0] == '@BASE_HINTENT':
                self.base_hintent = int(sline[1])

            if len(sline) == 2 and sline[0] == '@BASE_VINTENT':
                self.base_vintent = int(sline[1])

            if len(sline) == 2 and sline[0] == '@VINTENT_MARK':
                self.vintent_mark = sline[1]

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

        print(' ')
        print(' # default setting ')
        print('  - capitailize            : {}'.format(self.use_capital))
        print('  - base horizontal intent : {}'.format(self.base_hintent))
        print('  - base vertical   intent : {}'.format(self.base_vintent))
        print('  - vertical intent mark   : {}'.format(self.vintent_mark))
        print(' ')

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




    def find_char_index(self, string):
        ustring = string.upper()
        idx = 6

        for i in range(self.n_chars):
            if ustring == self.name[i]:
                idx = i
                break
            elif self.name[i] == '1234' and ustring.isdigit():
                idx = i
                break

        return idx




    def replace_line_before(self, line):
        # replace
        new_line = line
        new_line = new_line.replace('(', ' ( ')
        new_line = new_line.replace(')', ' ) ')
        new_line = new_line.replace('*', ' * ')
        new_line = new_line.replace('+', ' + ')
        new_line = new_line.replace('-', ' - ')
        new_line = new_line.replace('/', ' / ')
        new_line = new_line.replace('**', ' ** ')
        new_line = new_line.replace('=', ' = ')
        new_line = new_line.replace(':', ' : ')
        new_line = new_line.replace(',', ' , ')
        new_line = new_line.replace('==', ' == ')
        new_line = new_line.replace('!=', ' != ')
        new_line = new_line.replace('>', ' > ')
        new_line = new_line.replace('>=', ' >= ')
        new_line = new_line.replace('<', ' < ')
        new_line = new_line.replace('<=', ' <= ')
#        new_line = new_line.replace('.EQ.', ' .EQ. ')
#        new_line = new_line.replace('.NE.', ' .NE. ')
#        new_line = new_line.replace('.GT.', ' .GT. ')
#        new_line = new_line.replace('.GE.', ' .GE. ')
#        new_line = new_line.replace('.LT.', ' .LT. ')
#        new_line = new_line.replace('.LE.', ' .LE. ')
#        new_line = new_line.replace('.NOT.', ' .NOT. ')
#        new_line = new_line.replace('.NOT.', ' .NOT. ')
        return new_line



    def replace_line_after(self, line):
        new_line = line
        new_line = new_line.replace('= =', '==')
        new_line = new_line.replace('= >', '=>')
        new_line = new_line.replace('< =', '<=')
        new_line = new_line.replace(': :', '::')
        new_line = new_line.replace(' (', '(')
        new_line = new_line.replace('( ', '(')
        new_line = new_line.replace(' )', ')')
        new_line = new_line.replace(') %', ')%')
        new_line = new_line.replace(' ,', ',')
        # add
        new_line = new_line.replace(' => ', '=>')
        new_line = new_line.replace('/ /', '//')
        new_line = new_line.replace(' // ', '//')
        new_line = new_line.replace('len = ', 'len=')
        new_line = new_line.replace('kind = ', '')
        new_line = new_line.replace(' .and. ', '.and.')
        new_line = new_line.replace(' .and.', '.and.')
        new_line = new_line.replace('.and. ', '.and.')
        new_line = new_line.replace(' .or. ', '.or.')
        new_line = new_line.replace(' .or.', '.or.')
        new_line = new_line.replace('.or. ', '.or.')
        new_line = new_line.replace(' .not. ', '.not.')
        new_line = new_line.replace(' .not.', '.not.')
        new_line = new_line.replace('.not. ', '.not.')
        new_line = new_line.replace(' != ', '!=')
        new_line = new_line.replace(' !=', '!=')
        new_line = new_line.replace('!= ', '!=')
        new_line = new_line.replace(' .ne. ', '.ne.')
        new_line = new_line.replace(' .ne.', '.ne.')
        new_line = new_line.replace('.ne. ', '.ne.')
        new_line = new_line.replace(' == ', '==')
        new_line = new_line.replace(' ==', '==')
        new_line = new_line.replace('== ', '==')
        new_line = new_line.replace(' .eq. ', '.eq.')
        new_line = new_line.replace(' .eq.', '.eq.')
        new_line = new_line.replace('.eq. ', '.eq.')
        new_line = new_line.replace(' > ', '>')
        new_line = new_line.replace(' >', '>')
        new_line = new_line.replace('> ', '>')
        new_line = new_line.replace(' .gt. ', '.gt.')
        new_line = new_line.replace(' .gt.', '.gt.')
        new_line = new_line.replace('.gt. ', '.gt.')
        new_line = new_line.replace(' >= ', '>')
        new_line = new_line.replace(' .ge. ', '.ge.')
        new_line = new_line.replace(' .ge.', '.ge.')
        new_line = new_line.replace('.ge. ', '.ge.')
        new_line = new_line.replace(' < ', '<')
        new_line = new_line.replace(' <', '<')
        new_line = new_line.replace('< ', '<')
        new_line = new_line.replace(' .lt. ', '.lt.')
        new_line = new_line.replace(' .lt.', '.lt.')
        new_line = new_line.replace('.lt. ', '.lt.')
        new_line = new_line.replace(' <= ', '<=')
        new_line = new_line.replace(' <=', '<=')
        new_line = new_line.replace('<= ', '<=')
        new_line = new_line.replace(' .le. ', '.le.')
        new_line = new_line.replace(' .le.', '.le.')
        new_line = new_line.replace('.le. ', '.le.')
        new_line = new_line.replace('un = ', 'un=')
        new_line = new_line.replace('len = ', 'len=')
        new_line = new_line.replace('intent(in)', 'intent(in   )')
        new_line = new_line.replace('intent(out)', 'intent(  out)')
        new_line = new_line.replace('dim = ', 'dim=')
        new_line = new_line.replace(':, :', ':,:')
        new_line = new_line.replace(':, :', ':,:')
        new_line = new_line.replace('end do', 'enddo')
        new_line = new_line.replace('else if', 'elseif')
        new_line = new_line.replace('end if', 'endif')
        new_line = new_line.replace('if(', 'if (')
        new_line = new_line.replace('elseif(', 'elseif (')
        new_line = new_line.replace('only : ', 'only: ')
        new_line = new_line.replace(' : ', ':')
        new_line = new_line.replace(' + ', '+')
        new_line = new_line.replace(' +', '+')
        new_line = new_line.replace('+ ', '+')
        new_line = new_line.replace(' - ', '-')
        new_line = new_line.replace(' -', '-')
        new_line = new_line.replace('- ', '-')
        new_line = new_line.replace(' * ', '*')
        new_line = new_line.replace(' *', '*')
        new_line = new_line.replace('* ', '*')
        new_line = new_line.replace(' / ', '/')
        new_line = new_line.replace(' /', '/')
        new_line = new_line.replace('/ ', '/')
#        new_line = new_line.replace('intent', '                      intent')
#
        '''
        new_line = re.sub('^%s'%'   end subroutine', '!\n   return\n   end subroutine', new_line)
        new_line = re.sub('^%s'%'   end function', '!\n   return\n   end function', new_line)
        new_line = re.sub('^%s'%'   implicit none','!-------------------------------------------------------------------------------\n   implicit none\n!',new_line)
        new_line = re.sub('^%s'%'   subroutine','!-------------------------------------------------------------------------------\n!\n!\n!-------------------------------------------------------------------------------\n   subroutine',new_line)
        new_line = re.sub('^%s'%'   function','!-------------------------------------------------------------------------------\n!\n!\n!-------------------------------------------------------------------------------\n   function',new_line)
#        new_line = re.sub('^%s'%'   module', '!-------------------------------------------------------------------------------\n   module', new_line)
        new_line = re.sub('^%s'%'   end module', '!-------------------------------------------------------------------------------\n!\n!\n!-------------------------------------------------------------------------------\n   end module', new_line)
        #new_line = re.sub('^%s'%' ', ' ', new_line)
        '''

        return new_line



    def get_leading_type(self, upper_line, split_upper_line):
        indexes  = []
        types    = []
        needs    = []
        gotos    = []
        ends     = []
        renames  = []
        nn       = len(split_upper_line)

        # make types array
        for i in range(nn):
            ig = self.find_char_index(split_upper_line[i])
            #if ig >= 0:
            indexes.append(i)
            types.append(self.type[ig])
            needs.append(self.need[ig])
            gotos.append(self.goto[ig])
            ends.append(self.end[ig])
            renames.append(self.rename[ig])

        # determine the type, is_end (end character)
        nmatched = len(types)

        if self.is_cont:
            itype  = -2
            is_end = False

        elif nmatched == 0 or indexes[0] != 0: # otherwise type = 6: normal
            itype  = 6
            is_end = False

        else:
            # check need statement
            if needs[0] == 'none' or upper_line.count(needs[0]) > 0:
                itype = types[0]
            else:
                itype = gotos[0]

            if ends[0] == 0:
                is_end = False
            else:  # else, case 
                is_end = True

        return itype, is_end




    def end_type(self, linetype, line):

        # enddo, endif, end, ...
        if linetype == -1:

            # in find 'if', 'do', 'interface', ... in line types
            if self.types.count(self.typ_idx) > 0:

                # else, endif
                if self.ctype == self.typ_idx:
                    before_intent = self.intents[-1]
                    tmp = self.types.pop()
                    tmp = self.intents.pop()
                # normal case
                else:
                    before_intent = self.intents[-2]
                    tmp = self.types.pop()
                    tmp = self.types.pop()
                    tmp = self.intents.pop()
                    tmp = self.intents.pop()

                if self.types[-1] == self.typ_idx:
                    self.intents[-1] = before_intent

            # normal case
            else:

                while self.types[-1] != 0 and self.types[-1] != 1 and self.types[-1] != 2:
                    tmp = self.types.pop()
                    tmp = self.intents.pop()

                if line.strip() != 'CONTAINS' and line.strip() != 'contains':
                    before_intent = self.intents[-1]
                    tmp = self.types.pop()
                    tmp = self.intents.pop()
                    if self.types[-1] != 0: self.intents[-1] = before_intent

        # else
        else:   #if linetype != -1:
            if self.types[-1] != self.typ_idx:
                tmp = self.types.pop()
                tmp = self.intents.pop()
                self.types[-1] = linetype

        self.ctype     = self.types[-1]
        self.cintent_h = self.intents[-1]
        new_line       = ' '*self.cintent_h + line

        return new_line




    def preprocess(self, line):
        new_line = ''

        if line.lstrip()[0] == '!' or line.lstrip()[0] == '#':

            itype     = -3
            is_end    = False
            new_line  = line.strip() + '\n'

        else:

            line     = self.replace_line_before(line)
            s_line   = line.split()
            u_line   = line.upper()
            l_line   = line.lower()
            u_s_line = u_line.split()
            l_s_line = l_line.split()
            nstrs    = len(u_s_line)
            itype, is_end = self.get_leading_type(u_line, u_s_line)

            # code or comment
            for i in range(nstrs):
                if s_line[i][0] == "'":
                  new_line = new_line + s_line[i].lstrip() + ' '
                else:
                  new_line = new_line + l_s_line[i].lstrip() + ' '

            # base horizontal intent
            if not self.is_cont:
                new_line = ' '*self.base_hintent + new_line

            # continuation
            if self.merge_cont:
                if new_line.count('&') > 0:
                    idx          = new_line.index('&')
                    new_line     = new_line[:idx]
                    self.is_cont = True
                else:
                    new_line     = new_line.rstrip()+'\n'
                    self.is_cont = False
            else:
                if new_line.count('&') > 0:
                    new_line     = new_line.rstrip()+'\n'
                    if self.is_cont:
                       i_cont       = new_line.index('&')
                       if i_cont < self.max_width:
                           new_line = ' '*(self.max_width-i_cont) + new_line
                    self.is_cont = True
                else:
                    new_line     = new_line.rstrip()+'\n'
                    if self.is_cont:
                        i_cont = len(new_line)
                        if i_cont < self.max_width:
                            new_line = ' '*(self.max_width-i_cont) + new_line
                    self.is_cont = False

            '''
            if new_line.count('&') > 0:

                if self.merge_cont: # method 1 : remove after '&'

                    idx          = new_line.index('&')
                    new_line     = new_line[:idx]

                else:               # method 2 : not remove

                    new_line     = new_line.rstrip()+'\n'
                    i_cont       = new_line.index('&')
                    if i_cont < 80:
                        new_line = ' '*(80-i_cont) + new_line

                self.is_cont = True

            else:

                new_line     = new_line.rstrip()+'\n'
                self.is_cont = False
            '''

            new_line = self.replace_line_after(new_line)
    
        return itype, is_end, new_line



    def set_horizontal_intent(self, itype, is_end, line):
        otype = self.otype

        if   itype == -3:
            new_line = line
        elif itype == -2:
            new_line = line
        elif itype == -1 or is_end:
            new_line = self.end_type(itype, line)
        else:
            if self.types_info[otype][itype][1] == -1:
                print 'could not change types... ', otype, itype
                quit()

            self.ctype = itype

            if self.types_info[otype][itype][0]:
                self.cintent_h = self.cintent_h + self.types_info[otype][itype][1]
            else:
                self.cintent_h = self.types_info[otype][itype][1]

            if itype != otype or itype == self.typ_idx:
                self.types.append(self.ctype)
                self.intents.append(self.cintent_h)

            new_line = ' '*self.cintent_h + line

        # for debug
        if self.debug:
            print 'line    = ', line.strip()
            print 'type    = ', itype
            print 'types   = ', self.types
            print 'intents = ', self.intents
            print 'ctype   = ', self.ctype
            print 'cintent = ', self.cintent_h
            print ' '
            new_line = ':%02d:%02d:'%(itype, self.ctype) + new_line

        return new_line


 
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

        new_line = (self.vintent_mark+'\n')*intent_v + line

        return new_line


 
    def apply_kim_convension(self, itype, is_end, line):
  
        otype    = self.otype
        ctype    = self.ctype
        new_line = line
        #print('otype = {}, ctype = {}'.format(self.otype,self.ctype))

        if otype==0 and ctype==1:
            new_line = self.add_one_line(new_line)
        if otype==1 and ctype==0:
            new_line = self.add_four_line(new_line)
        if otype==1 and ctype==2:
            new_line = self.add_four_line(new_line)
        if otype==1 and ctype==3:
            #new_line = self.add_one_line(new_line)
            new_line = self.add_history(new_line)
        if otype==2 and ctype==3:
            new_line = self.add_one_line(new_line)
        if line.count('end module')>0:
            new_line = new_line+self.draw_one_line()
        if ctype==6 or (otype==6 and ctype==5):
            new_line = new_line.replace(', ',',')

        return new_line




    def draw_one_line(self):
        return '!'+'-'*(self.max_width-1)+'\n'

    def add_one_line(self, line):
        return self.draw_one_line()+line

    def add_four_line(self, line):
        #return '!\n'+2*('!'+'-'*(self.max_width-1)+'\n')+'!\n'+line
        return self.draw_one_line()+'!\n!\n!\n'+self.draw_one_line()+line

    def add_history(self, line):
        head = \
        '''!-------------------------------------------------------------------------------
!
!  abstract: 
!
!  history log :
!    201x-0x-xx  your   name    code comment
!
!  structure:
!
!-------------------------------------------------------------------------------
'''
        return head + line



    def generate_line(self, line):

        if line.strip() == '':
            new_line = None
        else:
            self.otype = self.ctype
            itype, is_end, new_line = self.preprocess(line)
            new_line                = self.set_horizontal_intent(itype, is_end, new_line)
            new_line                = self.set_vertical_intent(new_line)
            new_line                = self.apply_kim_convension(itype, is_end, new_line)

        return new_line
  



class F90Converter:



    def __init__(self, infilepath, oufilepath, extension=[], exclude=[]):
        self.infilepath = infilepath
        self.oufilepath = oufilepath
        self.extension  = extension
        self.exclude    = exclude
        self.grammer    = F90Grammer('/home/jhkim/work/python/fortran/v5.0/fortran.dat')



    def driver(self):
        infiles = os.listdir(self.infilepath)

        print(' ')
        print(' # start f90 converter')
        for ifile in infiles:
            infilename = self.infilepath+'/'+ifile
            oufilename = self.oufilepath+'/'+ifile

            if not os.path.isdir(infilename):
                print('  - processing file : {}'.format(infilename))
                self.convert(infilename, oufilename)
            else:
                print infilename, ' is directory'
        print(' ')



    def convert(self, infilename, oufilename):
        infile = open(infilename, 'r')
        oufile = open(oufilename, 'w')
        oldline = ''
        iscont  = False

        for line in infile:
            new_line = self.grammer.generate_line(line)
            if new_line != None:
                oufile.write(new_line)
            else:
                pass
                #oufile.write('\n')

        infile.close()
        oufile.close()
        self.grammer.reset()




def main():

    opt = OptionParser()
    # action: 'store', 'store_const', 'append', 'count', 'callback'
    opt.add_option('-s','--source',     dest='src', default='./src',action='store', help='source')
    opt.add_option('-d','--destination',dest='dst', default='./dst',action='store', help='destination')
    (options, args) = opt.parse_args()
    
    converter = F90Converter(options.src, options.dst)
    converter.driver()



if __name__ == '__main__':
    main()

    print(" ")
    print(' # your jobs')
    print("   0) write history with new format")
    print("   1) align statements: only, dimension, intent...")
    print("   2) align width.. (max: 80)")
    print("   3) clean comments")
    print("   4) add newlines")
    print(" ")




