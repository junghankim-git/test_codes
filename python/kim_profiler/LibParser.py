import os
import sys
import subprocess
sys.path.append('/home/jhkim/Study/Library/Shared/Python')

class Parser:
   def __init__(self, string=''):
      if string == '':
        self.string  = ''
        self.nLines  = 0
        self.isIni   = False

        self.bstring = []
        self.isBlcks = False
      else:
        self.string  = string
        self.nLines  = string.count('\n')
        self.isIni   = True

        self.bstring = []
        self.isBlcks = False


   def Ini(self, string):
      self.string = string
      self.nLines = string.count('\n')
      self.isIni  = True


   def IniFromProcessor(self, process):
      stringarray = process.readlines()
      nlines      = len(stringarray)

      string = ''
      for i in range(nlines):
         string = string + stringarray[i]

      self.string = string
      self.nLines = nlines
      self.isIni  = True

   def IsIni(self):
      if not self.isIni:
         print 'need to initialize the parser...'
         quit()


   def GetCounts(self, fword):
      self.IsIni()
      ncounts    = self.string.count(fword)
      return ncounts


   def MakeBlockedString(self, div_word):
      self.IsIni()
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
      self.nBlcks   = nblcks
      self.bstring  = blockString
      self.isBlcks  = True


   # string = 'xxxxxxxx sword resulttxt eword': return =>  resulttxt 
   def GetValueString(self, sword, eword, num=0):
      self.IsIni()
      return self.GetValue(self.string, sword, eword, num)


   def GetValueBlock(self, iblck, sword, eword, num=0):
      self.IsIni()
      return self.GetValue(self.bstring[iblck], sword, eword, num)


   def GetValue(self, string, sword, eword, num=0):
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


   def GetValueLine_Pattern(self, line, pword, sword, eword, num=0):
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


   def PrintString(self):
      print self.string


   def PrintBlock(self, iblck):
      print self.bstring[iblck]





#### old version ####


   # string = 'xxxxxxxx sword resulttxt eword': return =>  resulttxt 
   def GetValue_old(self, sword, eword):
      nlines     = self.string.count('\n')
      splitline  = self.string.splitlines()
      res_string = ''
      for iline in range(nlines):
         splitline[iline] = splitline[iline]+'\n'
         s          = splitline[iline].find(sword)
         e          = splitline[iline].find(eword)
         if s != -1 and e != -1:
            res_string = splitline[iline][s+len(sword):e]
      return res_string


   def GetValueFromBlock_old(self, iblck, sword, eword):
      nlines     = self.bstring[iblck].count('\n')
      splitline  = self.bstring[iblck].splitlines()
      res_string = ''
      for iline in range(nlines):
         splitline[iline] = splitline[iline]+'\n'
         s          = splitline[iline].find(sword)
         e          = splitline[iline].find(eword)
         if s != -1 and e != -1:
            res_string = splitline[iline][s+len(sword):e]
      return res_string

   def ProcessToString_old(self, nlines, process):
      string = ''
      for i in range(nlines):
         line   = process.readline()
         string = string + line
      return string
   
   
   def StringToBlock_old(self, string, div_word):
      nlines      = string.count('\n')
      nblcks      = string.count(div_word)
      lineString  = string.splitlines()
      blockString = ['' for i in range(nblcks)]
      iblck       = 0
      isStart     = False
      for iline in range(nlines):
         if lineString[iline][0:len(div_word)] == div_word:
            if isStart:
               iblck = iblck + 1
            else:
               isStart = True
         blockString[iblck] = blockString[iblck]+lineString[iline]+'\n'
      return blockString


   def IsOneWord_old(self, string, word):
      s = string.find(sword)
      if s == -1:
         return True
      else:
         return False

   # string = 'xxxxxxxx sword resulttxt eword': return =>  resulttxt 
   def GetValue_old(self, string, sword, eword):
      nlines     = string.count('\n')
      splitline  = string.splitlines()
      res_string = ''
      for iline in range(nlines):
         splitline[iline] = splitline[iline]+'\n'
         s          = splitline[iline].find(sword)
         e          = splitline[iline].find(eword)
         if s != -1 and e != -1:
            res_string = splitline[iline][s+len(sword):e]
      return res_string
