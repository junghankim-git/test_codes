import os
import sys
from inspect import currentframe, getframeinfo

INFO  = 0
WARN  = -1
FATAL = -2

def logger(message='', info=0, method=None, class_in=None):
    if info == WARN:
        info_string = '(  Warn)'
    elif info == FATAL:
        info_string = '( Error)'
    else:
        info_string = '(  Info)'
    classstring   = ''
    methodstring  = ''
    if not class_in==None:
        classstring = class_in.__class__.__name__
    else:
        classstring = 'None Class'
    if not method==None:
        methodstring = method
    else:
        methodstring  = 'None Method'
    twidth = 28
    width  = twidth - len(classstring) - len(methodstring)
    print('{0:7s} [ {1:s}:{2:s}{3:{width:d}s} ]:: {4:s}'.format(info_string,classstring,methodstring,' '*width,message,width=width))
    if info == -2:
       quit()




class ErrorHandler:
    def __init__(self):
        self.infos = [0, -1, -2] # 0:info,  -1:warning,  -2:fatal
 
    def Logger(self, class_in, method_name_in, message_in, info):
        if info == WARN:
          info_string = '@-1@'
        elif info == FATAL:
          info_string = '@-2@'
        else:
          info_string = '@--@'
  
        print info_string+' ['+class_in.__class__.__name__+'] ' + method_name_in + ':: ' + message_in
  
        if info == -2:
           quit()
 
    def GetMethodName(self):
        caller = currentframe().f_back
        func_name = getframeinfo(caller)[2]
        return func_name
