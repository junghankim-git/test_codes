import os
import sys
import subprocess

class parser:
    def __init__(self, string=''):
        if string == '':
          self.string  = ''
          self.nlines  = 0
          self.is_init = False
          self.bstring = []
          self.is_blck = False
        else:
          self.string  = string
          self.nlines  = string.count('\n')
          self.is_init = True
          self.bstring = []
          self.is_blck = False




    def initialize(self, string):
        self.string  = string
        self.nlines  = string.count('\n')
        self.is_init = True

 
 
    def initialize_proc(self, process):
        stringarray = process.readlines()
        nlines      = len(stringarray)
  
        string = ''
        for i in range(nlines):
            string = string + stringarray[i]
  
        self.string  = string
        self.nlines  = nlines
        self.is_init = True

 
    def initialized(self):
        if not self.is_init:
            print 'need to initialize the parser...'
            quit()
 
 
    def get_counts(self, fword):
        self.initialized()
        ncounts = self.string.count(fword)
        return ncounts
 
 
    def make_blocked_string(self, div_word):
        self.initialized()
        nlines      = self.string.count('\n')
        nblcks      = self.string.count(div_word)
        lineString  = self.string.splitlines()
        blockString = ['' for i in range(nblcks)]
        iblck       = 0
        isStart     = False
        for iline in range(nlines):
            #if lineString[iline][0:len(div_word)] == div_word:
            if lineString[iline].find(div_word) >= 0:
                if isStart:
                    iblck = iblck + 1
                else:
                    isStart = True
            blockString[iblck] = blockString[iblck]+lineString[iline]+'\n'
        self.nblcks   = nblcks
        self.bstring  = blockString
        self.is_blck  = True
 
 
    # string = 'xxxxxxxx sword resulttxt eword': return =>  resulttxt 
    def get_value_string(self, sword, eword, num=0):
        self.initialized()
        return self.get_value(self.string, sword, eword, num)
 
 
    def get_value_block(self, iblck, sword, eword, num=0):
        self.initialized()
        return self.get_value(self.bstring[iblck], sword, eword, num)
 
 
    def get_value(self, string, sword, eword, num=0):
        s = string.find(sword)
        e = string.find(eword, s)
        for n in range(num):
            s = string.find(sword, s+1)
            e = string.find(eword, e+1)
  
        if s != -1 and e != -1:
            res_string = string[s+len(sword):e]
        else:
            res_string = '-1: none'
        return res_string
 
 
    def get_value_line_pattern(self, line, pword, sword, eword, num=0):
        found = line.find(pword)
        if found >= 0:
            s = line.find(sword)
            e = line.find(eword, s)
            for n in range(num):
                s = line.find(sword, s+1)
                e = line.find(eword, e+1)
            if s != -1 and e != -1:
                res_string = line[s+len(sword):e]
            else:
                res_string = ''
        else:
            res_string = ''
        return res_string
 
 
    def print_string(self):
        print self.string
 
 
    def print_block_string(self, iblck):
        print self.bstring[iblck]
 
 
 
