#!/usr/bin/env python

import os
import sys
import re
from optparse import OptionParser
from f90convert import *




class F90Grammer:

    def __init__(self, filename, continuation, var_cont):
        self.use_capital  = False
        self.continuation = continuation
        self.var_cont     = var_cont
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
        self.shell_idx    = 5
        self.max_width    = 80




    def reset(self):
        self.types     = [0]     # all   type stateus
        self.intents   = [0]     # all intent stateus
        self.isendtype = False   # lately end type?
        self.ctype     = 0       # current type
        self.cintent_h = 0       # current horizontal intent




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
        new_line = new_line.replace(' /= ', '/=')
        new_line = new_line.replace(' /=', '/=')
        new_line = new_line.replace('/= ', '/=')
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

        if nmatched == 0 or indexes[0] != 0: # otherwise type = 6: normal
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




    def is_shell_lately_type(self):
        if self.types[-1] == self.shell_idx or self.types[-1] == 0 or \
           self.types[-1] == 1              or self.types[-1] == 2:
            return True
        else:
            return False




    def find_idx_lately_shell(self):

        it = -1; i0 = -1; i1 = -1; i2 = -1
        if self.types.count(self.shell_idx) > 0:
            it = len(self.types)-self.types[::-1].index(self.shell_idx)-1
        if self.types.count(0) > 0:  i0 = len(self.types)-self.types[::-1].index(0)-1
        if self.types.count(1) > 0:  i1 = len(self.types)-self.types[::-1].index(1)-1
        if self.types.count(2) > 0:  i2 = len(self.types)-self.types[::-1].index(2)-1
        ii = max(it,i0,i1,i2)

        return ii




    def end_type(self, linetype, line):

        # enddo, endif, endfunction, contains, ...
        if linetype == -1:

            # for multiple shell
            ii = self.find_idx_lately_shell()
            before_intent = self.intents[ii]

            # remove contents util front of the shell
            while not self.is_shell_lately_type():
                tmp = self.types.pop()
                tmp = self.intents.pop()

            if line.strip() != 'CONTAINS' and line.strip() != 'contains':
                # remove a shell
                tmp = self.types.pop()
                tmp = self.intents.pop()

            self.isendtype = True


        # else, case
        else:   #if linetype != -1:

            while not self.is_shell_lately_type():
                tmp = self.types.pop()
                tmp = self.intents.pop()
                self.types[-1] = linetype
            before_intent = self.intents[-1]

            self.isendtype = False


        # determine the current type and horizontal intent
        self.ctype     = self.types[-1]
        if self.types[-1] == self.shell_idx: # if before == 5
            self.cintent_h = before_intent
        else:
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
            new_line = ' '*self.base_hintent + new_line.rstrip() + '\n'

            new_line = self.replace_line_after(new_line)
    
        return itype, is_end, new_line



    def apply_horizontal_intent(self, itype, is_end, line):
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

            if not self.isendtype: # not after endif, enddo, ...
              if self.types_info[otype][itype][0]:
                  self.cintent_h = self.cintent_h + self.types_info[otype][itype][1]
              else:
                  self.cintent_h = self.types_info[otype][itype][1]

            if itype != otype or itype == self.shell_idx:
                self.types.append(self.ctype)
                self.intents.append(self.cintent_h)

            self.isendtype = False
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


 
    def apply_vertical_intent(self, line):

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
            new_line = self.add_history(new_line)
        if otype==2 and ctype==3:
            new_line = self.add_one_line(new_line)
        if line.count('end module')>0 or line.count('end program')>0:
            new_line = new_line+self.draw_one_line()
        if ctype==6 or (otype==6 and ctype==5):
            new_line = new_line.replace(', ',',')

        return new_line




    def draw_one_line(self):
        return '!'+'-'*(self.max_width-1)+'\n'

    def add_one_line(self, line):
        return self.draw_one_line()+line

    def add_four_line(self, line):
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




    def apply_continuation(self, itype, line):

        if itype == -1 or itype == -2 or itype == -3:
            return line

        if self.var_cont and (itype == 3 or itype == 4):
            # devide head and lines
            if itype == 3:
                div_str = 'only:'
                ndivs = line.count(div_str)
            elif itype == 4 and (len(line) >= 9 and line[3:9] != 'public'):
                div_str = '::'
                ndivs = line.count(div_str)
            else:
                ndivs = 0

            if ndivs == 1:
                lines = line.split(div_str)
                head  = lines[0]+div_str+' &'
                nls   = 1
                lines = lines[1]
            else:
                head  = line
                nls   = 0
                lines = ''

            # processing head
            comment = ''
            if head[-1] == '&':

                comment, head = self.split_comment(head)

                if len(head) <= self.max_width:
                    head = head[:-1]+' '*(self.max_width-len(head))+'&'

            # processing lines
            if nls > 0:
                nn      = lines.count(',')+1
                mlines  = lines.split(',')
                lines   = ''
                pos_var = 60
                for i in range(nn):
                    if i == nn-1:
                        mlines[i] = ' '*pos_var+mlines[i]
                        lines = lines + mlines[i]
                    else:
                        mlines[i] = ' '*pos_var+mlines[i]+','+' '*(self.max_width-len(mlines[i])-pos_var-2)
                        lines = lines + mlines[i]+'&'+'\n'

            # new_lines = head + lines
            new_line = comment + head
            for i in range(nls):
                new_line = new_line+'\n'+lines


        # normal continuation
        else:

            line = line[:-1]
            new_line = self.generate_newline(line)
            new_line = new_line + '\n'

        return new_line




    def generate_newline(self, line):

        # split comment and nomal line(remove \n)
        comment = ''
        if line.count('\n') > 0:
            comment, line = self.split_comment(line)

        new_line = comment + self.add_cont_line(line)

        return new_line




    def add_cont_line(self, line):
        '''
        " abc, def, gec, .............., abc, def, gec .............., abc, def, gec "
        to
        " abc, def, gec, .............., &
          abc, def, gec, .............., &
          abc, def, gec "
        '''

        comment = ''
        if   line.count('!') == 0:
            pass
        elif line.count('!') == 1:
            ic = line.index('!')
            comment = line[ic:]
            line    = line[:ic]
        else:
            print line
            print 'add_cont_line: too many ! mark'

        new_line = ''
        if len(line) > self.max_width:
            idx      = self.find_idx_near_max(line)
            before   = line[:idx+1]+' '*(self.max_width-len(line[:idx+1])-1)+'&'
            after    = 15*' '+line[idx+1:]
            if len(after) < self.max_width+3: # final line
               after    = (self.max_width-len(after)-1)*' ' + after
            new_line = before + '\n' + self.add_cont_line(after) + comment
        else:
            new_line = line + comment

        return new_line




    def find_idx_near_max(self, line):

        divisor = ['//', ',', '+', '-', ';', '=']
        ndivs   = len(divisor)
        has_div = True
        idxs    = [0 for i in range(ndivs)]
        idxs_n  = [0 for i in range(ndivs)]
        for i in range(ndivs):
            n = line.count(divisor[i])
            base = 0
            for j in range(n):
                pos  = line.index(divisor[i],idxs[i]+1)
                if pos > self.max_width-2:
                    idxs_n[i] = pos
                    break
                idxs[i] = pos
        idx = max(idxs)
        if idx == 0: idx = max(idxs_n)
        return idx




    def split_comment(self, line):
        '''
        '!\n!\n!------------\n   my code line'
         ->
                 comment         +     new_line
        '!\n!\n!------------\n'  + '   my code line'
        '''
        nhoh  = line.count('\n')

        comment  = ''
        new_line = line
        if nhoh > 0:
            split_line = line.split('\n')
            for i in range(nhoh):
                comment = comment + split_line[i] + '\n'
            new_line = split_line[-1]

        return comment, new_line




    def remove_continuation(self, line, new_line):

        if   new_line.count('&') == 0:
            is_cont = False
        elif new_line.count('&') == 1:

            is_comment = False
            if new_line.count("'") == 2:
                i = new_line.index("'")
                j = new_line.index("'",i+1)
                m = new_line.index("&")
                if m > i and m < j:
                    is_cont  = False
                    is_comment = True
            elif new_line.count('"') == 2:
                i = new_line.index('"')
                j = new_line.index('"',i+1)
                m = new_line.index('&')
                if m > i and m < j:
                    is_cont  = False
                    is_comment = True
           
            if not is_comment:
                if line.strip()[0] == '&' and not line.strip()[-1] == '&':
                    new_line = new_line.replace('&','')
                    is_cont  = False
                else:
                    idx      = new_line.index('&')
                    new_line = new_line[:idx]
                    is_cont  = True
        elif new_line.count('&') == 2:
            new_line     = new_line.replace('&','',1)
            idx          = new_line.index('&')
            new_line     = new_line[:idx]
            is_cont = True
        else:
            print("too many '&' mark...")
            quit()

        return is_cont, new_line



    def convert(self, infilename, oufilename):
        infile = open(infilename, 'r')
        oufile = open(oufilename, 'w')

        is_cont = False
        for line in infile:
            if line.strip() == '':
                continue

            # remove the continuation
            if is_cont:
                if line.strip()[0] == '!': line = '    &\n'
                new_line = new_line + line.strip()
            else:
                new_line = line.strip()

            if new_line.strip()[0] != '#' and new_line.strip()[0] != '!':
                is_cont, new_line = self.remove_continuation(line, new_line)
            if is_cont: continue

            self.otype = self.ctype
            itype, is_end, new_line = self.preprocess(new_line)
            new_line                = self.apply_horizontal_intent(itype, is_end, new_line)
            new_line                = self.apply_vertical_intent(new_line)
            new_line                = self.apply_kim_convension(itype, is_end, new_line)
            if self.continuation:
                new_line            = self.apply_continuation(itype, new_line)
            oufile.write(new_line)

        infile.close()
        oufile.close()
        self.reset()



class F90Converter:



    def __init__(self, infilepath, oufilepath, continuation=False, var_cont=False, extension=[], exclude=[]):
        self.infilepath = infilepath
        self.oufilepath = oufilepath
        self.extension  = extension
        self.exclude    = exclude
        self.grammer    = F90Grammer('/home/jhkim/work/python/fortran/v8.0/fortran.dat', continuation, var_cont)
        if not os.path.isdir(oufilepath):
            os.system('mkdir -p '+oufilepath)



    def driver(self):
        infiles = os.listdir(self.infilepath)

        print(' ')
        print(' # start f90 converter')
        for ifile in infiles:
            infilename = self.infilepath+'/'+ifile
            oufilename = self.oufilepath+'/'+ifile

            if not os.path.isdir(infilename):
                print('  - processing file : {}'.format(infilename))
                self.grammer.convert(infilename, oufilename)
            else:
                print infilename, ' is directory'
        print(' ')




def main():

    opt = OptionParser()
    # action: 'store', 'store_const', 'append', 'count', 'callback'
    opt.add_option('-s','--source',     dest='src', default='./src',action='store',      help='source')
    opt.add_option('-d','--destination',dest='dst', default='./dst',action='store',      help='destination')
    opt.add_option('-c','--continue',   dest='cont',default=False,  action='store_true', help='continuation')
    opt.add_option('-v','--variable',   dest='var', default=False,  action='store_true', help='var_cont')
    (options, args) = opt.parse_args()
    
    converter = F90Converter(options.src, options.dst, options.cont, options.var)
    converter.driver()

    print(" ")
    print(' # your jobs')
    print("   0) write history with new format")
    print("   1) align statements: only, dimension, intent...")
    print("   2) clean comments")
    print("   3) add newlines (!)")
    print("   4) align width.. (max: 80)")
    print(" ")



if __name__ == '__main__':
    main()




