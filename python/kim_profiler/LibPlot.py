import os
import sys
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import cm
sys.path.append('/home/jhkim/Study/Library/Shared/Python')

class MatPlotLib:

   def __init__(self,colsize,rowsize,title,tsize=18,nrows=1,ncols=1):
      self.nrows  = nrows
      self.ncols  = ncols
      self.plt    = plt
      self.pltfig = plt.figure(figsize=(colsize,rowsize))
      self.pltfig.suptitle(title,fontsize=tsize)
      self.axis = [[self.pltfig.add_subplot(nrows,ncols,ncols*i+j+1) for j in range(ncols)] for i in range(nrows)]
      self.plot = [[None for j in range(ncols)] for i in range(nrows)]

      # format
      self.ncolors  = 7
      self.colors   = ['k', 'b', 'g', 'r', 'y', 'm', 'c']
      self.nmarkers = 8
      self.markers  = ['o', 'v', 's', '8', 'p', 'h', 'x', '+']
      self.nlines   = 2
      self.lines    = ['-', '--']

   def IniAxis(self,irow,icol,xmin,xmax,ymin,ymax,title,xlabel,ylabel,lgrid=False,lframe=True):
      if irow >= self.nrows:
         eh.Logger(self, self.IniAxis.__name__, 'Check \'irow\'...', -2)
      if icol >= self.ncols:
         eh.Logger(self, self.IniAxis.__name__, 'Check \'icol\'...', -2)
      self.axis[irow][icol].clear()
      self.axis[irow][icol].set_xlim(xmin,xmax)
      self.axis[irow][icol].set_ylim(ymin,ymax)
      self.axis[irow][icol].set_title(title)
      self.axis[irow][icol].set_xlabel(xlabel)
      self.axis[irow][icol].set_ylabel(ylabel)
      if lgrid: self.axis[irow][icol].grid(True)
      if lframe == False: self.axis[irow][icol].set_frame_on(False)
      return self.axis[irow][icol]

   def IniAxisSpec(self,axis,xmin,xmax,ymin,ymax,title,xlabel,ylabel,lframe=True):
      axis.clear()
      axis.set_xlim(xmin,xmax)
      axis.set_ylim(ymin,ymax)
      axis.set_title(title)
      axis.set_xlabel(xlabel)
      axis.set_ylabel(ylabel)
      if lframe == False: axis.set_frame_on(False)
   
   
   def Plot(self,irow,icol,xval,yval,fmt='k-',lw=1.0,label=None): 
      #self.axis[irow][icol].plot(xval,yval,fmt,label=label)
      self.plot[irow][icol], = self.axis[irow][icol].plot(xval,yval,fmt,lw=lw,label=label)
      return self.plot[irow][icol]
   
   
   def Bar(self,irow,icol,xval,yval,bottom=0.0,color=None,edgecolor=None,hatch=None,label=None): 
      self.plot[irow][icol] = self.axis[irow][icol].bar(xval,yval,bottom=bottom,color=color,edgecolor=edgecolor,hatch=hatch,label=label)
      return self.plot[irow][icol]
   
   
   def PlotSetData(self,irow,icol,xval,yval): 
      self.plot[irow][icol].set_data(xval,yval)
   
   
   def PlotError(self,irow,icol,xval,yval,xerr=None,yerr=None,fmt='k-',label=None): 
      self.axis[irow][icol].errorbar(xval,yval,xerr=xerr,yerr=yerr,fmt=fmt,label=label)
   
   
   def PlotVLine(self,irow,icol,x,ymin,ymax,c='k',ls='-', label=None):
      ret = self.axis[irow][icol].axvline(x=x, ymin=ymin, ymax=ymax, c=c, ls=ls, label=label)
      return ret
   
   def PlotHLine(self,irow,icol,y,xmin,xmax,c='k',ls='-', label=None):
      ret = self.axis[irow][icol].axhline(y=y, xmin=xmin, xmax=xmax, c=c, ls=ls, label=label)
      return ret


   def PlotLegend(self,irow,icol,loc,fontsize=None,shadow=True):
      self.axis[irow][icol].legend(loc=loc,fontsize=fontsize,shadow=shadow)


   def GetPlotFormats(self,ndims,nfmts,ncolors=None,nmarkers=None):

      ifc = [0 for idim in range(ndims)]
      ifm = [0 for idim in range(ndims)]
      ifl = [0 for idim in range(ndims)]
      if ndims == 1:
         fmt_out = ['' for ifmt in range(nfmts)]
         for ifmt in range(nfmts):
            if ifc[0] == self.ncolors:  ifc[0] = 0
            if ifm[0] == self.nmarkers: ifm[0] = 0
            if ifl[0] == self.nlines:   ifl[0] = 0
            fmt_out[ifmt] = self.colors[ifc[0]]+self.markers[ifm[0]]+self.lines[ifl[0]]
            if (ifmt+1)/self.ncolors != ifm[0]: ifm[0] = ifm[0] + 1
            ifc[0] = ifc[0] + 1
      elif ndims == 2:
         fmt_out = [['' for ifmt in range(nfmts)] for idim in range(ndims)]
         for idim in range(ndims):
            for ifmt in range(nfmts):
               if ifc[idim] == self.ncolors:  ifc[idim] = 0
               if ifm[idim] == self.nmarkers: ifm[idim] = 0
               if idim == 0:
                  fmt_out[idim][ifmt] = self.colors[ifc[idim]]+self.markers[ifm[idim]]+self.lines[0]
               else:
                  fmt_out[idim][ifmt] = self.colors[ifc[idim]]+self.markers[ifm[idim]]+self.lines[1]
               if (ifmt+1)/self.ncolors != ifm[idim]: ifm[idim] = ifm[idim] + 1
               ifc[idim] = ifc[idim] + 1
      else:
         eh.Logger(self, self.GetPlotFormats.__name__, 'Check \'ndims\'...', -2)

      return fmt_out


   def SaveFigure(self,filename):
      self.pltfig.savefig(filename)


   def Show(self, isShow=None):
      plt.show(isShow)


   def Draw(self):
      plt.draw()
