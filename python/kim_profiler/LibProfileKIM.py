# -*- coding: utf-8 -*- import sys
import sys
sys.path.append('/home/jhkim/Study/Library/Shared/Python')
from LibParser import *
from LibPlot import *
sys.path.append('/home/jhkim/Study/Library/Shared/Python/ProfileKIM')
from LibProfComp import *
#######################################################################################
# Description: Profiling of KIM
#
# Date: 13MAR2015
#   - JH KIM: First, for KIM version 0.21.01
#
# Date: 18MAR2015
#   - JH KIM: KIM version 0.21.02
#             Add Subroutine: WrtieFile, ReadFile
#
# Date: 24MAR2015
#   - JH KIM: KIM version 0.21.02
#             Modified for Difference and formats
#
# Date 03APR2015
#   - JH KIM: Initailize with input path (Initialize_DIR, Dir2CoresDirs)
#
#######################################################################################


class ProfileKIM:

   def __init__(self, vert='Eulerian'):
      self.isKeys   = False
      self.isCores  = False
      self.isData   = False
      self.nCores   = -1
      self.FileDir  = ''
      self.FileHead = 'ProfileInfo'
      self.FileTail = '.prf'


      # Component's names and Parsing information (key, dependent on KIM version)
      # Model
      names_model   = ['Set', 'Ini', 'Run', 'Fin']
      nkeys_mod     = [1, 1, 1, 1]
      keys_mod      = [['SetAtmosModel '], ['IniAtmosModel '], ['RunAtmosModel '], ['FinAtmosModel ']]
      keys_sign_mod = [[1], [1], [1], [1]]
 
      # Run
      names_run     = ['Dynamics', 'Physics', 'Communication', 'Input/Output', 'Synchronize']
      nkeys_run     = [12, 2, 6, 7, 5]
      keys_run      = [['RunCoreHOMME ', 'sync_prim_advance_exp ', 'sync_compute_and_apply_rhs ', 'sync_advance_hypervis ', 'sync_prim_advec_tracers_remap_k2 ', \
                              'sync_advance_hypervis_scalar ', 'bndry_exchange(compute_and_apply_rhs) ', 'bndry_exchange(biharmonic_wk) ', \
                              'bndry_exchange(advance_hypervis) ', 'bndry_exchange(euler_step) ', 'bndry_exchange(biharmonic_wk_scalar) ', \
                              'bndry_exchange(advance_hypervis_scalar) '], \
                       ['RunPhysicsPackage ', 'WritePhysicsOutput(RunT) '], \
                       ['bndry_exchange(compute_and_apply_rhs) ', 'bndry_exchange(biharmonic_wk) ', 'bndry_exchange(advance_hypervis) ', \
                        'bndry_exchange(euler_step) ', 'bndry_exchange(biharmonic_wk_scalar) ', 'bndry_exchange(advance_hypervis_scalar) '], \
                       ['WriteDyCoreOutput(Ini) ','WriteHOMMEOutput(Ini) ','WritePhysicsOutput(Ini) ',\
                              'WriteDyCoreOutput(Run) ','WriteHOMMEOutput(Run) ','WritePhysicsOutput(RunT) ', 'WriteRestart '], \
                       ['sync_prim_advance_exp ', 'sync_compute_and_apply_rhs ', 'sync_advance_hypervis ', 'sync_prim_advec_tracers_remap_k2 ', 'sync_advance_hypervis_scalar ']]
      if vert != 'Eulerian':
         keys_run[0][7] = 'bndry_exchange(biharmonic_wk_dp3d) '
         keys_run[0][8] = 'bndry_exchange(advance_hypervis_dp) '
         keys_run[2][1] = 'bndry_exchange(biharmonic_wk_dp3d) '
         keys_run[2][2] = 'bndry_exchange(advance_hypervis_dp) '
      keys_sign_run = [[1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1], [1,-1], [1,1,1,1,1,1], [1,1,1,1,1,1,1], [1,1,1,1,1,1]]

      # Dynamics
      names_dyn     = ['Synchronize (dyn-phy)', 'Governing Equation', 'Hyper-viscosity', 'Tracer Advection', 'Remapping']
      nkeys_dyn     = [1, 1, 1, 1, 1]
      keys_dyn      = [['sync_prim_advance_exp '], ['compute_and_apply_rhs                  '], ['advance_hypervis                        '], ['prim_advec_tracers_remap_rk2             '], ['vertical_remap ']]
      keys_sign_dyn = [[1], [1], [1], [1], [1]]
      if vert != 'Eulerian':
         keys_dyn[2][0] = 'advance_hypervis_dp '

      # Physics
      names_phy     = ['Radiation', 'Sfc', 'Land', 'PBL', 'GWDOrg', 'ConvDeep', 'ConvShal', 'CldMacro', 'CldMicro', 'CldMacroAfMicro', 'GWDNonOro']
      nkeys_phy     = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
      keys_phy      = [['RunRad '], ['RunSfc '], ['RunLand '], ['RunPbl '], ['RunGWDOro '], ['RunConvDeep '], \
                       ['RunConvShal '], ['RunCldMacro '], ['RunCldMicro '], ['RunCldMacroAfMicro '], ['RunGWDNonOro ']]
      keys_sign_phy = [[1], [1], [1], [1], [1], [1], [1], [1], [1], [1], [1]]



      # Dimensions
      nModels = len(names_model)
      nRuns   = len(names_run)
      nDyns   = len(names_dyn)
      nPhys   = len(names_phy)


      # Components
      self.Model    = Component('Model',    names_model)
      self.Run      = Component('Run',      names_run)
      self.Dynamics = Component('Dynamics', names_dyn)
      self.Physics  = Component('Physics',  names_phy)


      # Set Keys
      self.Model.SetKeys(nkeys_mod, keys_mod, keys_sign_mod)
      self.Run.SetKeys(nkeys_run, keys_run, keys_sign_run)
      self.Dynamics.SetKeys(nkeys_dyn, keys_dyn, keys_sign_dyn)
      self.Physics.SetKeys(nkeys_phy, keys_phy, keys_sign_phy)


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
      self.Run.IniCores(nCores, isDiff=isDiff)
      self.Dynamics.IniCores(nCores, isDiff=isDiff)
      self.Physics.IniCores(nCores, isDiff=isDiff)

      # Ini from Files
#      self.Model.ReadData_files(file_dirs, self.FileHead, self.FileTail)
#      self.Run.ReadData_files(file_dirs, self.FileHead, self.FileTail)
#      self.Dynamics.ReadData_files(file_dirs, self.FileHead, self.FileTail)
#      self.Physics.ReadData_files(file_dirs, self.FileHead, self.FileTail)
      self.Model.ReadData(file_dirs, self.FileHead, self.FileTail, isDiff=isDiff)
      self.Run.ReadData(file_dirs, self.FileHead, self.FileTail, isDiff=isDiff)
      self.Dynamics.ReadData(file_dirs, self.FileHead, self.FileTail, isDiff=isDiff)
      self.Physics.ReadData(file_dirs, self.FileHead, self.FileTail, isDiff=isDiff)


   def Initialize_DB(self, file_dir=None):
      # Ini from Files
      self.Model.ReadFile(file_dir)
      self.Run.ReadFile(file_dir)
      self.Dynamics.ReadFile(file_dir)
      self.Physics.ReadFile(file_dir)



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
      self.Run.PrintScreen()
      self.Dynamics.PrintScreen()
      self.Physics.PrintScreen()


   def PrintWiki(self):
      # Print Wiki
      self.Model.PrintWiki()
      self.Run.PrintWiki()
      self.Dynamics.PrintWiki()
      self.Physics.PrintWiki()


   def PrintWiki_Diff(self):
      # Print Wiki
      self.Model.PrintWiki_Diff()
      self.Run.PrintWiki_Diff()
      self.Dynamics.PrintWiki_Diff()
      self.Physics.PrintWiki_Diff()


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
      self.Run.PrintKeysInfo()
      print '}}}'
      print ' '
      print '* Dynamics'
      print '{{{'
      self.Dynamics.PrintKeysInfo()
      print '}}}'
      print ' '
      print '* Physics'
      print '{{{'
      self.Physics.PrintKeysInfo()
      print '}}}'


   def WriteFile(self,filedir=''):
      # Print Screen
      self.Model.WriteFile(filedir)
      self.Run.WriteFile(filedir)
      self.Dynamics.WriteFile(filedir)
      self.Physics.WriteFile(filedir)




   def DoPlot(self, tag='tmp'):
      import os
      print '## Plotting         ... '
      figdir = './Figs'
      os.system('mkdir -p '+figdir)
      #####################
      # Plotting
      #####################
      nrows = 1
      ncols = 3
      
      titles  = ['Wall-clock Time', 'Scalability', 'Ratio']
      xlabels = ['Number of prosesses [#]', 'Number of prosesses [#]', 'Number of prosesses [#]']
      ylabels = ['Wall-clock Time [s]', 'Speed-up', 'Ratio']
      
      xmins   = [0.0, 0.0, 0.0]
      xmaxs   = [self.Model.Cores[-1]*1.1, self.Model.Cores[-1]*1.1, self.Run.nCores]
      ymins   = [0.0, 0.0, 0.0]
      ymaxs   = [self.Model.Mean[-1][0]*1.1, self.Model.Cores[-1]*1.1, 101.0]
      
      canvas1 = MatPlotLib(20,6,'',18,nrows,ncols)
      axis11 = canvas1.IniAxis(0,0,xmins[0],xmaxs[0],ymins[0],ymaxs[0],titles[0],xlabels[0],ylabels[0])
      axis12 = canvas1.IniAxis(0,1,xmins[1],xmaxs[1],ymins[1],ymaxs[1],titles[1],xlabels[1],ylabels[1])
      axis13 = canvas1.IniAxis(0,2,xmins[2],xmaxs[2],ymins[2],ymaxs[2],titles[2],xlabels[2],ylabels[2])
      
      fmts   = ['k--','k*-','r*-','b*-','g*-','c*-','y*-']
      labels = ['Ideal', 'Model', 'Dynamics', 'Physics', 'Communication', 'Input/Output','Synchronize']

      # Plot 1-1
      canvas1.Plot(0,0,self.Model.Cores,self.Model.Mean[-1],fmts[1],label=labels[1])
      for i in range(5):
         canvas1.Plot(0,0,self.Model.Cores,self.Run.Mean[i],fmts[i+2],label=labels[i+2])
      canvas1.PlotLegend(0,0,loc=1,fontsize=12)


      # Plot 1-2
      canvas1.Plot(0,1,self.Model.Cores,self.Model.Cores,fmts[0],label=labels[0])
      canvas1.Plot(0,1,self.Model.Cores,self.Model.Speed[-1],fmts[1],label=labels[1])
      for i in range(5):
         canvas1.Plot(0,1,self.Model.Cores,self.Run.Speed[i],fmts[i+2],label=labels[i+2])
      canvas1.PlotLegend(0,1,loc=2,fontsize=12)

      # Plot 1-3
      labels = ['Dynamics', 'Physics', 'Communication', 'Input/Output','Synchronize']
      colors = ['r', 'b', 'g', 'c', 'y']
      xaxis  = [i   for i in range(self.Run.nCores)]
      bottom = [0.0 for i in range(self.Run.nCores)]
      for i in range(len(labels)):
         canvas1.Bar(0,2,xaxis,self.Run.Ratio[i],bottom=bottom,color=colors[i],edgecolor='k',label=labels[i])
         for j in range(self.Run.nCores):
            bottom[j] = bottom[j] + self.Run.Ratio[i][j]
      ticks      = [i+0.5 for i in range(self.Run.nCores)]
      ticklabels = [self.Run.Cores[i] for i in range(self.Run.nCores)]
      axis13.set_xticks(ticks)
      axis13.set_xticklabels(ticklabels)
      canvas1.PlotLegend(0,2,loc=3,fontsize=12)
      
      canvas1.SaveFigure(figdir+'/Result_Model_'+tag+'.eps')




      xmins   = [0.0, 0.0, 0.0]
      xmaxs   = [self.Model.Cores[-1]*1.1, self.Model.Cores[-1]*1.1, self.Dynamics.nCores]
      ymins   = [0.0, 0.0, 0.0]
      ymaxs   = [self.Dynamics.Mean[-1][2]*1.1, self.Dynamics.Cores[-1]*1.1, 101.0]
      
      canvas2 = MatPlotLib(20,6,'',18,nrows,ncols)
      axis21 = canvas2.IniAxis(0,0,xmins[0],xmaxs[0],ymins[0],ymaxs[0],titles[0],xlabels[0],ylabels[0])
      axis22 = canvas2.IniAxis(0,1,xmins[1],xmaxs[1],ymins[1],ymaxs[1],titles[1],xlabels[1],ylabels[1])
      axis23 = canvas2.IniAxis(0,2,xmins[2],xmaxs[2],ymins[2],ymaxs[2],titles[2],xlabels[2],ylabels[2])
      
      fmts   = ['k--','k*-','r*-','b*-','g*-','c*-','y*-']
      labels = ['Ideal', 'Dynamics', 'Synchronize', 'Governnig', 'Hyper-viscosity', 'Tracer Advection','Remapping']

      # Plot 2-1
      canvas2.Plot(0,0,self.Model.Cores,self.Dynamics.Mean[-1],fmts[1],label=labels[1])
      for i in range(5):
         canvas2.Plot(0,0,self.Model.Cores,self.Dynamics.Mean[i],fmts[i+2],label=labels[i+2])
      canvas2.PlotLegend(0,0,loc=1,fontsize=12)

      # Plot 2-2
      canvas2.Plot(0,1,self.Model.Cores,self.Model.Cores,fmts[0],label=labels[0])
      canvas2.Plot(0,1,self.Model.Cores,self.Dynamics.Speed[-1],fmts[1],label=labels[1])
      for i in range(5):
         canvas2.Plot(0,1,self.Model.Cores,self.Dynamics.Speed[i],fmts[i+2],label=labels[i+2])
      canvas2.PlotLegend(0,1,loc=2,fontsize=12)
      
      # Plot 2-3
      labels = ['Synchronize', 'Governnig', 'Hyper-viscosity', 'Tracer Advection','Remapping']
      colors = ['r', 'b', 'g', 'c', 'y']
      xaxis  = [i   for i in range(self.Dynamics.nCores)]
      bottom = [0.0 for i in range(self.Dynamics.nCores)]
      for i in range(len(labels)):
         canvas2.Bar(0,2,xaxis,self.Dynamics.Ratio[i],bottom=bottom,color=colors[i],edgecolor='k',label=labels[i])
         for j in range(self.Dynamics.nCores):
            bottom[j] = bottom[j] + self.Dynamics.Ratio[i][j]
      ticks      = [i+0.5 for i in range(self.Dynamics.nCores)]
      ticklabels = [self.Dynamics.Cores[i] for i in range(self.Dynamics.nCores)]
      axis23.set_xticks(ticks)
      axis23.set_xticklabels(ticklabels)
      canvas2.PlotLegend(0,2,loc=3,fontsize=12)
      
      canvas2.SaveFigure(figdir+'/Result_Dynamics_'+tag+'.eps')



      # Show
      canvas1.Show()


   def DoPlot_Diff(self):
      pass

