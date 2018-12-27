# -*- coding: utf-8 -*- import sys
import sys
sys.path.append('/home/jhkim/Study/Library/Shared/Python')
from LibParser import *
from LibPlot import *
sys.path.append('/home/jhkim/Study/Library/Shared/Python/ProfileKIM')
from LibProfComp import *

class ProfileVAR:

   def __init__(self):
      self.isKeys   = False
      self.isCores  = False
      self.isData   = False
      self.nCores   = -1
      self.FileDir  = ''
      self.FileHead = 'ProfileInfo'
      self.FileTail = '.prf'


      # Component's names and Parsing information (key, dependent on KIM version)
      # Model
      names_model   = ['Initialize', 'Observation', 'Outer Loop', 'Other']
      nkeys_mod     = [1, 1, 1, 4]
      keys_mod      = [['Initialize HybDA'], ['Observation'], ['Run HybDA'], ['Total','Initialize HybDA','Observation','Run HybDA']]
      keys_sign_mod = [[1], [1], [1], [1,-1,-1,-1]]
 
      # Outer Loopt
      names_otl     = ['Make RHS b', 'Check Adj', 'Make LHS Ax0', 'InnerLoop', 'TlmToModel', 'Other']
      nkeys_otl     = [1, 1, 1, 1, 3, 8]
      keys_otl      = [['Make_RHS_b'],['Check_Adj'],['Make_LHS_Ax0'],['Innerloop'],['TlmModelToEig_outer 3','TlmEigToModel_outer', 'Check_RMSE'],['Run HybDA', 'Make_RHS_b', 'Check_Adj', 'Make_LHS_Ax0', 'Innerloop','TlmModelToEig_outer 3','TlmEigToModel_outer','Check_RMSE']]
      keys_sign_otl = [[1], [1], [1], [1], [1,1,1], [1,-1,-1,-1,-1,-1,-1,-1]]

      # Inner Loop
      names_inl     = ['Make LHS Ax', 'Cal CostFun', 'Other']
      nkeys_inl     = [1, 1, 3]
      keys_inl      = [['Make_LHS_Ax  '],['Cal_CostFun'],['Innerloop','Make_LHS_Ax  ','Cal_CostFun']]
      keys_sign_inl = [[1], [1], [1,-1,-1]]

      # Inner Loop
      names_lhs     = ['TlmEigToObs', 'AdjEigToObs', 'Other']
      nkeys_lhs     = [1, 1, 3]
      keys_lhs      = [['TlmEigToObs_Make_Ax 4'],['AdjEigToObs_Make_Ax 4'],['Make_LHS_Ax  ','TlmEigToObs_Make_Ax 4','AdjEigToObs_Make_Ax 4']]
      keys_sign_lhs = [[1], [1], [1,-1,-1]]

      # TlmEigToObs (Inner Loop)
      names_TEO     = ['Bcast', 'SpecToCtrl', 'TlmCtrlToModel', 'ModelToObs', 'Other']
      nkeys_TEO     = [1, 1, 1, 1, 5]
      keys_TEO      = [['Run Bcast4SpecToCtrl 4'],['Run SpecToCtrl 4'],['Run TlmCtrlToModel 4'],['Run TlmModelToObs 4'],['TlmEigToObs_Make_Ax 4','Run Bcast4SpecToCtrl 4','Run SpecToCtrl 4','Run TlmCtrlToModel 4','Run TlmModelToObs 4']]
      keys_sign_TEO = [[1], [1], [1], [1], [1,-1,-1,-1,-1]]

      # AdjEigToObs (Inner Loop)
      names_AEO     = ['AdjModelToObs', 'AdjCtrlToModel', 'AdjSpecToCtrl', 'Global Comm.', 'Other']
      nkeys_AEO     = [1, 1, 2, 1, 4]
      keys_AEO      = [['Run AdjModelToObs 4'],['Run AdjCtrlToModel 4'],['Run AdjSpecToCtrl 4','Wrap_Repro_Sum_fromAdjSpecToCtrl 4'], ['Wrap_Repro_Sum_fromAdjSpecToCtrl 4'],['AdjEigToObs_Make_Ax 4','Run AdjModelToObs 4','Run AdjCtrlToModel 4','Run AdjSpecToCtrl 4']]
      keys_sign_AEO = [[1], [1], [1,-1], [1], [1,-1,-1,-1]]

      # AdjSpecToCtrl (all)
      names_aAEO    = ['AdjSpecToCtrl']
      nkeys_aAEO    = [4]
      keys_aAEO     = [['Run AdjSpecToCtrl 1', 'Run AdjSpecToCtrl 2', 'Run AdjSpecToCtrl 3', 'Run AdjSpecToCtrl 4']]
      keys_sign_aAEO= [[1,1,1,1]]



      # Dimensions
#      nModels = len(names_model)
#      nRuns   = len(names_otl)
#      nDyns   = len(names_inl)
#      nTEOs   = len(names_TEO)
#      nAEOs   = len(names_AEO)


      # Components
      self.Model    = Component('VAR',        names_model)
      self.Outer    = Component('Outer Loop', names_otl)
      self.Inner    = Component('Inner Loop', names_inl)
      self.LHS      = Component('LHS (Ax)',   names_lhs)
      self.TEO      = Component('TlmEigToObs',names_TEO)
      self.AEO      = Component('AdjEigToObs',names_AEO)
      self.aAEO     = Component('AdjSpecToCtrl',names_aAEO)


      # Set Keys
      self.Model.SetKeys(nkeys_mod, keys_mod, keys_sign_mod)
      self.Outer.SetKeys(nkeys_otl, keys_otl, keys_sign_otl)
      self.Inner.SetKeys(nkeys_inl, keys_inl, keys_sign_inl)
      self.LHS.SetKeys(nkeys_lhs, keys_lhs, keys_sign_lhs)
      self.TEO.SetKeys(nkeys_TEO, keys_TEO, keys_sign_TEO)
      self.AEO.SetKeys(nkeys_AEO, keys_AEO, keys_sign_AEO)
      self.aAEO.SetKeys(nkeys_aAEO, keys_aAEO, keys_sign_aAEO)


      if self.isKeys:
         print 'WORNING :: rewrite keys....'
      self.isKeys   = True

      print ''



   def Initialize(self, nCores, file_dirs, isDiff=False):

      # Set Cores
      if not isDiff:
         self.nCores   = nCores
         self.FileDir  = file_dirs

      self.Model.IniCores(nCores, isDiff=isDiff)
      self.Outer.IniCores(nCores, isDiff=isDiff)
      self.Inner.IniCores(nCores, isDiff=isDiff)
      self.LHS.IniCores(nCores, isDiff=isDiff)
      self.TEO.IniCores(nCores, isDiff=isDiff)
      self.AEO.IniCores(nCores, isDiff=isDiff)
      self.aAEO.IniCores(nCores, isDiff=isDiff)

      # Ini from Files
      self.Model.ReadData(file_dirs, self.FileHead, self.FileTail, isDiff=isDiff, is1File=True)
      self.Outer.ReadData(file_dirs, self.FileHead, self.FileTail, isDiff=isDiff, is1File=True)
      self.Inner.ReadData(file_dirs, self.FileHead, self.FileTail, isDiff=isDiff, is1File=True)
      self.LHS.ReadData(file_dirs, self.FileHead, self.FileTail, isDiff=isDiff, is1File=True)
      self.TEO.ReadData(file_dirs, self.FileHead, self.FileTail, isDiff=isDiff, is1File=True)
      self.aAEO.ReadData(file_dirs, self.FileHead, self.FileTail, isDiff=isDiff, is1File=True)




   def Initialize_DB(self, file_dir=None):
      # Ini from Files
      self.Model.ReadFile(file_dir)
      self.Outer.ReadFile(file_dir)
      self.Inner.ReadFile(file_dir)
      self.LHS.ReadFile(file_dir)
      self.TEO.ReadFile(file_dir)
      self.aAEO.ReadFile(file_dir)



   def Initialize_DIR(self, file_dir, isDiff=False):
      if isinstance(file_dir, list):
         print 'Error: input was list...'
         quit()
      ncores, dirs = self.Dir2CoresDirs(file_dir)
      self.Initialize(ncores, dirs, isDiff)



   def Dir2CoresDirs(self, file_dir):
      files  = os.listdir(file_dir)
      nfiles = len(files)
      
      ncores = []
      dirs   = []
      
      for i in range(nfiles):
         filename = file_dir+'/'+files[i]+'/profile'
         if os.path.isdir(filename):
            try:
               nproc  = int(files[i])
            except:
               continue
            ncores.append(nproc)
            dirs.append(filename)
      
      ncores.sort()
      dirs.sort()
      return ncores, dirs



   def PrintScreen(self):
      # Print Screen
      self.Model.PrintScreen()
      self.Outer.PrintScreen()
      self.Inner.PrintScreen()
      self.LHS.PrintScreen()
      self.TEO.PrintScreen()
      self.AEO.PrintScreen()
      self.aAEO.PrintScreen()


   def PrintWiki(self):
      # Print Wiki
      self.Model.PrintWiki()
      self.Outer.PrintWiki()
      self.Inner.PrintWiki()
      self.LHS.PrintWiki()
      self.TEO.PrintWiki()
      self.AEO.PrintWiki()
      self.aAEO.PrintWiki()


   def PrintWiki_Diff(self):
      # Print Wiki
      self.Model.PrintWiki_Diff()
      self.Outer.PrintWiki_Diff()
      self.Inner.PrintWiki_Diff()
      self.TEO.PrintWiki_Diff()
      self.AEO.PrintWiki_Diff()
      self.aAEO.PrintWiki_Diff()


   def PrintKeysInfo(self):
      print '==== Keys Info. ===='
      print ' '
      print '* Model'
      print '{{{'
      self.Model.PrintKeysInfo()
      print '}}}'
      print ' '
      print '* Run'
      print '{{{'
      self.Outer.PrintKeysInfo()
      print '}}}'
      print ' '
      print '* Inner'
      print '{{{'
      self.Inner.PrintKeysInfo()
      print '}}}'
      print ' '
      print '* LHS'
      print '{{{'
      self.LHS.PrintKeysInfo()
      print '}}}'
      print ' '
      print '* TEO'
      print '{{{'
      self.TEO.PrintKeysInfo()
      print '}}}'
      print '* AEO'
      print '{{{'
      self.AEO.PrintKeysInfo()
      print '}}}'
      print '* aAEO'
      print '{{{'
      self.aAEO.PrintKeysInfo()
      print '}}}'


   def WriteFile(self,filedir=''):
      # Print Screen
      self.Model.WriteFile(filedir)
      self.Outer.WriteFile(filedir)
      self.Inner.WriteFile(filedir)
      self.LHS.WriteFile(filedir)
      self.TEO.WriteFile(filedir)
      self.AEO.WriteFile(filedir)
      self.aAEO.WriteFile(filedir)




   def DoPlot_Scalability(self):
      import os
      import matplotlib.pyplot as plt
      print '## Plotting         ... '
      figdir = './Figs'
      os.system('mkdir -p '+figdir)
      #####################
      # Plotting
      #####################
      nfigs = 2
      nrows = 1
      ncols = 3


      # Plot Setting
      pltfig = [plt.figure(figsize=(18,5)) for i in range(nfigs)]
      #pltfig = plt.figure(figsize=(18,9))
      axis = [[[pltfig[i].add_subplot(nrows, ncols, j*ncols+k+1) for k in range(ncols)] for j in range(nrows)] for i in range(nfigs)]
      
      titles  = [[['Wall-clock Time', 'Ratio', 'Scalability']],[['Wall-clock Time', 'Ratio', 'Time Reduction']]]
      xlabels = [[['Number of prosesses [#]' for k in range(ncols)] for j in range(nrows)] for i in range(nfigs)]
      ylabels = [[['Wall-clock Time [s]', 'Ratio', 'Speed-up']], [['Wall-clock Time [s]', 'Ratio [%]', 'Time Reduction [%]']]]
      xmins   = [[[self.Model.Cores[0]*0.9 for k in range(ncols)] for j in range(nrows)] for i in range(nfigs)]
      xmaxs   = [[[self.Model.Cores[-1]*1.1 for k in range(ncols)] for j in range(nrows)] for i in range(nfigs)]
      ymins   = [[[0.0, 0.0, self.Model.Cores[0]*0.9]], [[0.0, 0.0, 0.0]]]
      ymaxs   = [[[self.Model.Mean[-1][0]*1.1, 101.0, self.Model.Cores[-1]*1.1]], [[self.Model.Mean[-1][0]*1.1, 101.0, 101.0]]]
      for i in range(nfigs):
        for j in range(nrows):
          for k in range(ncols):
            axis[i][j][k].set_xlim(xmins[i][j][k],xmaxs[i][j][k])
            axis[i][j][k].set_ylim(ymins[i][j][k],ymaxs[i][j][k])
            axis[i][j][k].set_title(titles[i][j][k])
            axis[i][j][k].set_xlabel(xlabels[i][j][k])
            axis[i][j][k].set_ylabel(ylabels[i][j][k])
            axis[i][j][k].grid(True)
      


      # Data Setting
      for i in range(nfigs):
         xval   = [[None for k in range(ncols)] for j in range(nrows)]
         maxY   = 2
         yval   = [[[None for k in range(ncols)] for j in range(nrows)] for y in range(maxY)]
         yvalD  = [[[None for k in range(ncols)] for j in range(nrows)] for y in range(maxY)]
         nCores = self.Model.nCores
         for j in range(nrows):
            for k in range(ncols):
               xval[j][k]  = self.Model.Cores
               if j == 0:
                  if k == 0:
                     yval[0][j][k]  = self.Model.Mean[-1]
                     yval[1][j][k]  = self.aAEO.Mean[-1]
                     if i == 1:
                        yvalD[0][j][k] = self.Model.MeanD[-1]
                        yvalD[1][j][k] = self.aAEO.MeanD[-1]
                  if k == 1: 
                     yval[0][j][k]  = [0.0 for l in range(nCores)]
                     for l in range(nCores):
                        yval[0][j][k][l]  = self.aAEO.Mean[-1][l]/self.Model.Mean[-1][l]*100.0
                     if i == 1:
                        yvalD[0][j][k] = [0.0 for l in range(nCores)]
                        for l in range(nCores):
                           yvalD[0][j][k][l] = self.aAEO.MeanD[-1][l]/self.Model.MeanD[-1][l]*100.0
                        #print 'ratio = '
                        #print yval[0][j][k]
                        #print yvalD[0][j][k]
                  if k == 2:
                     if i == 0:
                        yval[0][j][k]  = self.Model.Speed[-1]
                        yval[1][j][k]  = self.aAEO.Speed[-1]
                     if i == 1:
                        yval[0][j][k]  = [0.0 for l in range(nCores)]
                        yval[1][j][k]  = [0.0 for l in range(nCores)]
                        for l in range(nCores):
                           yval[0][j][k][l] = (self.Model.Mean[-1][l]-self.Model.MeanD[-1][l])/self.Model.Mean[-1][l]*100.0
                           yval[1][j][k][l] = (self.aAEO.Mean[-1][l]-self.aAEO.MeanD[-1][l])/self.aAEO.Mean[-1][l]*100.0
                        #print 'time reduction = '
                        #print yval[0][j][k]
                        #print yval[1][j][k]
   
   
         # Plot
         locs = [1,2,1]
         for j in range(nrows):
            for k in range(ncols):
               if k == 0:
                  axis[i][j][k].plot(xval[j][k], yval[0][j][k], 'r*-', label='Original 3DVAR')
                  axis[i][j][k].plot(xval[j][k], yval[1][j][k], 'b*-', label='Original AdjSpecToCtrl')
                  if i == 1:
                    axis[i][j][k].plot(xval[j][k], yvalD[0][j][k], 'r*--', label='Optimized 3DVAR')
                    axis[i][j][k].plot(xval[j][k], yvalD[1][j][k], 'b*--', label='Optimized AdjSpecToCtrl')
               if k == 1:
                  axis[i][j][k].plot(xval[j][k], yval[0][j][k], 'r*-', label='Original Ratio')
                  if i == 1:
                     axis[i][j][k].plot(xval[j][k], yvalD[0][j][k], 'r*--', label='Optimized Ratio')
               if k == 2:
                  if i == 0:
                    axis[i][j][k].plot(xval[j][k], xval[j][k], 'k-.', label='Ideal')
                    axis[i][j][k].plot(xval[j][k], yval[0][j][k], 'r*-', label='Original 3DVAR')
                    axis[i][j][k].plot(xval[j][k], yval[1][j][k], 'b*-', label='Original AdjSpecToCtrl')
                    #axis[i][j][k].plot(xval[j][k], yvalD[0][j][k], 'r*--', label='Optimized 3DVAR')
                    #axis[i][j][k].plot(xval[j][k], yvalD[1][j][k], 'b*--', label='Optimized AdjSpecToCtrl')
                  if i == 1:
                    axis[i][j][k].plot(xval[j][k], yval[0][j][k], 'r*-', label='3DVAR')
                    axis[i][j][k].plot(xval[j][k], yval[1][j][k], 'b*-', label='AdjSpecToCtrl')

               axis[i][j][k].legend(loc=locs[j], fontsize=12)
   
         outfigname = figdir+'/performance_fig'+str(i)+'.png'
         pltfig[i].savefig(outfigname)
      plt.show()



   def DoPlot(self, tag='tmp'):
      import os
      import matplotlib
      print '## Plotting         ... '
      figdir = './Figs'
      os.system('mkdir -p '+figdir)
      #####################
      # Plotting
      #####################
      nrows = 1
      ncols = 6
      
      titles  = ['3D VAR', 'Outer Loop', 'Inner Loop', 'Cal. LHS(Ax)', 'TlmEigToObs', 'AdjEigToObs']
      xlabels = ['','','','','', '']
      ylabels = ['Ratio','','','','','']
      
      xmins   = [0.0]
      xmaxs   = [4.0]
      ymins   = [0.0]
      ymaxs   = [140.0]

#      self.Model    = Component('VAR',        names_model)
#      self.Outer    = Component('Outer Loop', names_otl)
#      self.Inner    = Component('Inner Loop', names_inl)
#      self.TEO      = Component('TlmEigToObs',names_TEO)
#      self.AEO      = Component('AdjEigToObs',names_AEO)
      

      canvas1 = MatPlotLib(18,6,'Profiling of the 3DVAR',18,nrows,ncols)
      axis    = [None for i in range(ncols)]
      colors  = ['r', 'b', 'g', 'c', 'y', 'k']


      labels = [['Initialize', 'Obs', 'Outer Loop', 'Others'],['Cal. RHS(b)', 'Check Adj', 'LHS(Ax0)', 'Inner Loop', 'TlmToModel', 'Others'],['Cal. LHS(Ax)', 'Cost Func', 'Others'],['TlmEigToObs','AdjeigToObs','Others'],['Broadcast','SpecToCtrl','TlmCtrlToModel','ModelToObs','Others'],['AdjModelToObs','AdjCtrlToModel','AdjSpecToCtrl (-)','Global Comm.','Others']]
      # Plot 1-1

      ####################
      ## Plots
      for ic in range(ncols):
         myaxis = canvas1.IniAxis(0,ic,xmins[0],xmaxs[0],ymins[0],ymaxs[0],titles[ic],xlabels[ic],ylabels[ic])
         axis[ic] = myaxis
         if ic == 0: comp = self.Model
         if ic == 1: comp = self.Outer
         if ic == 2: comp = self.Inner
         if ic == 3: comp = self.LHS
         if ic == 4: comp = self.TEO
         if ic == 5: comp = self.AEO

         bottom = 0.0
         for i in range(comp.nNames-1):
            canvas1.Bar(0,ic,1,comp.Ratio[i][0],bottom=bottom,color=colors[i],edgecolor='k',label=labels[ic][i])
            bottom = bottom + comp.Ratio[i][0]
         canvas1.PlotLegend(0,ic,loc=1,fontsize=10)

#         if ic == 2:
#            y1_TEO = 0.0
#            y2_TEO = self.TEO.Mean[-1][0]/comp.Mean[-1][0]*100.
#            canvas1.Bar(0,ic,1,y2_TEO,bottom=y1_TEO,color=colors[0],edgecolor='k')
#
#            #bottom = yyy
#            #yyy    = self.AEO.Mean[-1][0]/comp.Mean[-1][0]*100.0
#            y1_AEO = y2_TEO
#            y2_AEO = self.AEO.Mean[-1][0]/comp.Mean[-1][0]*100.0 
#            canvas1.Bar(0,ic,1,y2_AEO,bottom=y1_AEO,color=colors[0],edgecolor='k')


         #ticks      = [i+0.5 for i in range(self.Outer.nCores)]
         #ticklabels = [self.Outer.Cores[i] for i in range(self.Outer.nCores)]
         myaxis.set_xticks([])
         myaxis.set_xticklabels([])
         if ic != 0: myaxis.set_yticklabels([])
         myaxis.grid(True)

      ####################
      ## Lines
      trans = canvas1.pltfig.transFigure.inverted()
      xL = 1.8
      xR = 1.0

      for ic in range(ncols-1):
         if ic == 0:
            idx  = 2
            comp = self.Model
         elif ic == 1:
            idx  = 3
            comp = self.Outer
         elif ic == 2:
            idx  = 0
            comp = self.Inner
         elif ic == 3:
            idx  = 0
            comp = self.LHS
         elif ic == 4:
            idx  = 1
            comp = self.LHS
#         elif ic == 5:
#            comp = self.AEO

         sic = ic
         tic = ic + 1
         y1, y2 = 0.0, 0.0
         for i in range(idx):
            y1  = y1 + comp.Ratio[i][0]
         y2 = y1 + comp.Ratio[idx][0]

         if ic == 4:
            sic = ic - 1
            tic = ic + 1
            y1, y2 = 0.0, 0.0
            for i in range(idx):
               y1  = y1 + comp.Ratio[i][0]
            y2 = y1 + comp.Ratio[idx][0]


         colors = ['g', 'b', 'r', 'r', 'b']
         #coord1 = trans.transform(axis[0].transData.transform([xL, self.Model.Ratio[0][0]+self.Model.Ratio[1][0]]))
         coord1 = trans.transform(axis[sic].transData.transform([xL, y1]))
         coord2 = trans.transform(axis[tic].transData.transform([xR, 0.0]))
         line = matplotlib.lines.Line2D((coord1[0], coord2[0]),(coord1[1], coord2[1]), color=colors[ic], transform=canvas1.pltfig.transFigure)
         canvas1.pltfig.lines.append(line)
         #coord1 = trans.transform(axis[0].transData.transform([xL, self.Model.Ratio[0][0]+self.Model.Ratio[1][0]+self.Model.Ratio[2][0]]))
         coord1 = trans.transform(axis[sic].transData.transform([xL, y2]))
         coord2 = trans.transform(axis[tic].transData.transform([xR, 100.0]))
         line = matplotlib.lines.Line2D((coord1[0], coord2[0]),(coord1[1], coord2[1]), color=colors[ic], transform=canvas1.pltfig.transFigure)
         canvas1.pltfig.lines.append(line)

      
      canvas1.SaveFigure(figdir+'/Result_VAR_'+tag+'.eps')
      
      # Show
      canvas1.Show()





   def DoPrediction(self):


      # Total walltime & Global communication
      Total = self.Model.Mean[-1][0]
      Comm  = self.AEO.Mean[3][0]
      Ratio = Comm/Total*100.0
      print 'Total = {0:5.2f} s , Comm = {1:5.2f} [s] '.format(Total, Comm)
      print 'Ratio = {0:5.2f} %'.format(Ratio)
      print ' '

      nn = 11
      Improve = [float(i)*10.0 for i in range(nn)]
      newComm = [0.0 for i in range(nn)]
      newTota = [0.0 for i in range(nn)]
      Speedup = [0.0 for i in range(nn)]
 
      for i in range(nn):
         newComm[i] = Comm*(100.0 - Improve[i])/100.0
         newTota[i] = (Total - Comm) + newComm[i]
         Speedup[i] = Total/newTota[i]
      print Improve
      print newComm
      print Speedup


      #####################
      # Plotting
      #####################
      nrows = 1
      ncols = 1
      
      titles  = ['Speedup']
      xlabels = ['Improvement of Communication [%]']
      ylabels = ['Speedup [#]']
      
      xmins   = [0.0]
      xmaxs   = [100.0]
      ymins   = [0.8]
      ymaxs   = [1.5]
      
      canvas1 = MatPlotLib(8,6,'Prediction',18,nrows,ncols)
      axis11 = canvas1.IniAxis(0,0,xmins[0],xmaxs[0],ymins[0],ymaxs[0],titles[0],xlabels[0],ylabels[0])
      
      fmts   = ['k--','k*-','r*-','b*-','g*-','c*-','y*-']
      labels = ['Ideal', 'Model', 'Dynamics', 'Physics', 'Communication', 'Input/Output','Synchronize']

      # Plot 1-1
      axis11.grid(True)
      canvas1.Plot(0,0,Improve,Speedup,'r*-')


      print '## Plotting         ... '
      figdir = './Figs'
      os.system('mkdir -p '+figdir)
      canvas1.SaveFigure(figdir+'/Result_VAR_Prediction.eps')
