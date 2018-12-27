# -*- coding: utf-8 -*- import sys
import os
#######################################################################################
# Description: Profiling of KIM
#
# Date: 03JUL2015
#   - JH KIM: First, for KIMOUT
#
#######################################################################################


class kim_out:

   def __init__(self, version, ne=-1, compiler='all'):
      self.version  = version
      self.old_ver  = 'none'
      self.nlevs    = 6
      self.ndirs    = [0 for i in range(self.nlevs)]
      self.dirs     = [[] for i in range(self.nlevs)]
      self.basedir  = os.path.abspath('/data/jhkim/TestBed/KIM/Output')

      self.def_ne   = ['ne30', 'ne60', 'ne120']
      self.def_comp = ['gnu', 'intel', 'pgi']
      self.def_int  = ['10h', '02d']
      ndyns = 2
      nphys = 2
      self.def_dyn  = ['SH', 'SW']
      #self.def_phy  = ['G', 'K']
      self.def_phy  = ['G']
      self.def_exp  = ['SW_G', 'SH_G']
      self.def_qrp  = ['101','403']   # tracer, remap, physics
      
      nlevs         = self.nlevs
      
      #self.dirs[0] = [basedir+'/'+str(version)]
      self.dirs[0] = [str(version)]
      if ne <= 0:
        self.dirs[1] = self.def_ne
      else:
        self.dirs[1] = ['ne'+str(ne)]
      
      if compiler == 'all':
        self.dirs[2] = self.def_comp
      else:
        self.dirs[2] = [compiler]
      
      self.dirs[3] = self.def_int
      self.dirs[4] = self.def_exp
      self.dirs[5] = self.def_qrp
      
      for ilev in range(nlevs):
        self.ndirs[ilev] = len(self.dirs[ilev])



   # opt = 0(list dirs), 1(make dirs), 2(link dir), 3(copy dir), 4(check difference), 5(check performance)
   def ListDirs(self, nlev_in, ilev_in, basedir_in, subdirs, opt=-1, tpath=''):
      # START: basic algorithm.
      nsubdirs = len(subdirs[ilev_in])
      for idir in range(nsubdirs):
         full_path = basedir_in+'/'+subdirs[ilev_in][idir]

         # START: Options (depend on this object)
         if opt == 0:
            print full_path
         elif opt == 1:
            if os.path.exists(full_path):
               print full_path+': already existed...'
            else:
               print 'making ... '+full_path
               os.system("mkdir "+full_path)
         elif opt == 2 or opt == 3:
            if ilev_in == nlev_in-1:
               src_basedir = self.SourceDir(full_path)
               #print '## ', full_path, full_path[-8:-4], self.s_exps
               if self.s_exps.count(full_path[-8:-4]) > 0:
                  if os.path.exists(full_path) or os.path.lexists(full_path):
                     command = 'rm -rf '+full_path
                     print command
                     os.system(command)
                  if os.path.islink(src_basedir): src_basedir = os.path.realpath(src_basedir)
                  if opt == 2:
                     command = 'ln -s '+src_basedir+' '+full_path
                  else:
                     command = 'cp -r '+src_basedir+' '+full_path
                  print command
                  os.system(command)
         elif opt == 4 or opt == 5:
            if ilev_in == nlev_in-1 and subdirs[ilev_in][idir] == '101':
              new_path = full_path.replace(self.basedir, self.c_basedir)
              new_path = new_path.replace(self.version, self.c_version)
              if new_path.count("ne30/gnu/") > 0:
                 new_path = new_path.replace("ne30/gnu/","")
              elif new_path.count("ne30/intel/") > 0:
                 new_path = new_path.replace("ne30/intel/","")
              elif new_path.count("ne30/pgi/") > 0:
                 new_path = new_path.replace("ne30/pgi/","")
              else:
                 print 'check compiler in opt 4'
                 quit()
              if tpath != '':
                 if tpath[-1] != '/': tpath = tpath+'/'
                 new_path = tpath+new_path[-8:]
              self.CheckFile(full_path, new_path)
         # END  : Options

         if ilev_in != nlev_in-1:
            self.ListDirs(nlev_in, ilev_in+1, full_path, subdirs, opt, tpath)
      # END  : basic algorithm.


   def PrintDirs(self):
      nlevs = self.nlevs
      self.ListDirs(nlevs, 0, self.basedir, self.dirs, 0)


   def CreateDirs(self):
      nlevs = self.nlevs
      self.ListDirs(nlevs, 0, self.basedir, self.dirs, 1)
      os.system('touch '+self.dirs[0][0]+'/'+'log.txt')
      os.system('mkdir '+self.dirs[0][0]+'/'+'Figs')
      os.system('mkdir '+self.dirs[0][0]+'/'+'Performance')


   def ActionDirs(self, old_version, dynamics, physics, iact):
      nlevs = self.nlevs
      if dynamics == 'all':
         dynamics = self.def_dyn
      else:
         if self.def_dyn.count(dynamics) <= 0:
            print 'FATAL: check dynamics'
            exit()
         else:
            dynamics = [dynamics]
      if physics == 'all':
         physics = self.def_phy
      else:
         if self.def_phy.count(physics) <= 0:
            print 'FATAL: check physics'
            exit()
         else:
            physics = [physics]

      ndyns = len(dynamics)
      nphys = len(physics)
      self.old_ver = old_version
      self.s_nexps = ndyns * nphys
      self.s_exps  = ['' for i in range(self.s_nexps)]
      for id in range(ndyns):
         for ip in range(nphys):
            self.s_exps[id*nphys+ip] = dynamics[id]+'_'+physics[ip]
      #print self.s_exps
      if iact != 2 and iact != 3:
         print 'check action...'
         quit()
      self.ListDirs(nlevs, 0, self.basedir, self.dirs, iact)


   def SourceDir(self, target_in):
      source_out = target_in.replace(self.version, self.old_ver)
      return source_out







      



   # opt = 0(diff), 1(ncdump -h), 2(ncdump), 3(ncdump -v)
   def check_difference(self, path='none', version='none', integration='all', dynamics='all', physics='all', opt=0, var='u', ifile=0, tpath=''):
      nlevs = self.nlevs
#      self.c_comp = ['gnu']
      # path
      if path == 'none':
         self.c_basedir = '/scratch/jhkim/TestBed/Data'
      else:
         self.c_basedir = path
      if not os.path.isdir(self.c_basedir):
         print 'check dir:', self.c_basedir

      # version
      if version == 'none':
         print 'check c_version...', version
         quit()
      else:
         self.c_version = version
      
      # integration
      if integration == 'all':
           self.c_int = self.def_int
      else:
        if self.def_int.count(integration) <= 0:
            print 'FATAL: check integration', integration
            exit()
        elif integration == '10h':
           self.c_int = ['10h']
        elif integration == '2d' or integration == '02d':
           self.c_int = ['02d']
        else:
           print 'check integration ...', integration

      if dynamics == 'all':
         self.c_dyn = self.def_dyn
      else:
         if self.def_dyn.count(dynamics) <= 0:
            print 'FATAL: check dynamics', dynamics
            exit()
         else:
            self.c_dyn = [dynamics]
      if physics == 'all':
         self.c_phy = self.def_phy
      else:
         if self.def_phy.count(physics) <= 0:
            print 'FATAL: check physics', physics
            exit()
         else:
            self.c_phy = [physics]
  
      self.c_opt = opt
      self.c_var = var
      self.c_ifile = ifile
      self.c_files = self.makeCFile()

      print 'SOURCE DIR = ', self.basedir
      print 'TARGET DIR = ', self.c_basedir
      if tpath != '':
        print 'T Path     = ', tpath
      self.ListDirs(nlevs, 0, self.basedir, self.dirs, 4, tpath)



   def CheckFile(self, s_path, t_path):
      isInt = False; isDyn = False; isPhy = False
      for i in range(len(self.c_int)):
         if s_path.count(self.c_int[i]) > 0:
            isInt = True
            iint  = i
      for i in range(len(self.c_dyn)):
         if s_path.count(self.c_dyn[i]+'_') > 0:
            isDyn = True
            idyn  = i
      for i in range(len(self.c_phy)):
         if s_path.count('_'+self.c_phy[i]) > 0:
            isPhy = True
            iphy  = i
      if not isInt or not isDyn or not isPhy: return
      
      if s_path.count('10h') > 0:
        i_int = 0
      elif s_path.count('02d') > 0:
        i_int = 1
      else:
        print 'check int index...'
        quit()


      if self.c_opt != 4:
        src_file = s_path+'/'+self.c_files[i_int][self.c_ifile]
        trg_file = t_path+'/'+self.c_files[i_int][self.c_ifile]
      else:
        src_file = s_path+'/'+'profile/ProfileInfoSummary.prf'
        trg_file = t_path+'/'+'profile/ProfileInfoSummary.prf'
  
#      if not os.path.isfile(src_file) or not os.path.isfile(trg_file):
#         if not os.path.isfile(src_file):
#            print 'file not found...', src_file
#         if not os.path.isfile(trg_file):
#            print 'file not found...', trg_file
#         return
      if not os.path.isfile(src_file):
        new_src_file = src_file.replace('101/', '')
        if os.path.isfile(new_src_file):
          src_file = new_src_file
        else:
          print 'file not found...', src_file
          return

      if not os.path.isfile(trg_file):
        new_trg_file = trg_file.replace('101/', '')
        if os.path.isfile(new_trg_file):
          trg_file = new_trg_file
        else:
          print 'file not found...', trg_file
          return

      print '# START... ', self.c_int[iint], self.c_dyn[idyn], self.c_phy[iphy]
      if self.c_opt == 0:
        print ' - diff source: '+src_file
        print ' - diff target: '+trg_file
        command = 'diff '+src_file+' '+trg_file
        os.system(command)
        print ' '
        print '  => more info :: '
        print ' KIMStdOut.exe '+src_file+' '+trg_file
         
      elif self.c_opt == 1 or self.c_opt == 2 or self.c_opt == 3:
        os.system('rm -rf org.txt new.txt')
        if self.c_opt == 1: ncdump = 'ncdump -h '
        if self.c_opt == 2: ncdump = 'ncdump '
        if self.c_opt == 3: ncdump = 'ncdump -v '+self.c_var+' '
        command = ncdump+src_file+' > '+'org.txt'
        print command
        os.system(command)
        command = ncdump+trg_file+' > '+'new.txt'
        print command
        os.system(command)
        command = 'diff org.txt new.txt | tee diff.txt'
        print command
        os.system(command)
      elif self.c_opt == 4:
        par = Parser()
        pro = os.popen('cat '+src_file)
        par.IniFromProcessor(pro)
        src_walltime = float(par.GetValueString('Total          ', '\n').split()[2])
        pro = os.popen('cat '+trg_file)
        par.IniFromProcessor(pro)
        trg_walltime = float(par.GetValueString('Total          ', '\n').split()[2])
        print ' - source walltime = ', src_walltime
        print ' - target walltime = ', trg_walltime
        print ' - speed-up        = ', src_walltime/trg_walltime
        
      else:
         print 'check c option...', self.c_opt
      print '# END  ... ', self.c_int[iint], self.c_dyn[idyn], self.c_phy[iphy]
      print ' '


   def makeCFile(self):
      nint  = 2
      nfile = [11, 9]
      #start = 2011072512
      files = [['' for j in range(nfile[i])] for i in range(nint)]
      '''
      files[0] = ['UP-20110725120000-000000.nc',                                \
                  'UP-20110725130000-000001.nc', 'UP-20110725140000-000002.nc', \
                  'UP-20110725150000-000003.nc', 'UP-20110725160000-000004.nc', \
                  'UP-20110725170000-000005.nc', 'UP-20110725180000-000006.nc', \
                  'UP-20110725190000-000007.nc', 'UP-20110725200000-000008.nc', \
                  'UP-20110725210000-000009.nc', 'UP-20110725220000-000010.nc']
      files[1] = ['UP-20110725120000-000000.nc',                                \
                  'UP-20110725130000-000001.nc', 'UP-20110725140000-000002.nc', \
                  'UP-20110725150000-000003.nc', 'UP-20110725160000-000004.nc', \
                  'UP-20110725170000-000005.nc', 'UP-20110725180000-000006.nc', \
                  'UP-20110725190000-000007.nc', 'UP-20110725200000-000008.nc', \
                  'UP-20110725210000-000009.nc', 'UP-20110725220000-000010.nc']
      files[0] = ['UP-20170601000000-000000.nc',                                \
                  'UP-20170601010000-000001.nc', 'UP-20170601020000-000002.nc', \
                  'UP-20170601030000-000003.nc', 'UP-20170601040000-000004.nc', \
                  'UP-20170601050000-000005.nc', 'UP-20170601060000-000006.nc', \
                  'UP-20170601070000-000007.nc', 'UP-20170601080000-000008.nc', \
                  'UP-20170601090000-000009.nc', 'UP-20170601100000-000010.nc']
      files[1] = ['UP-20170601000000-000000.nc',                                \
                  'UP-20170601010000-000001.nc', 'UP-20170601020000-000002.nc', \
                  'UP-20170601030000-000003.nc', 'UP-20170601040000-000004.nc', \
                  'UP-20170601050000-000005.nc', 'UP-20170601060000-000006.nc', \
                  'UP-20170601070000-000007.nc', 'UP-20170601080000-000008.nc', \
                  'UP-20170601090000-000009.nc', 'UP-20170601100000-000010.nc']
      '''
      files[0] = ['UP-20171125120000-000000.nc',                                \
                  'UP-20171125130000-000001.nc', 'UP-20171125140000-000002.nc', \
                  'UP-20171125150000-000003.nc', 'UP-20171125160000-000004.nc', \
                  'UP-20171125170000-000005.nc', 'UP-20171125180000-000006.nc', \
                  'UP-20171125190000-000007.nc', 'UP-20171125200000-000008.nc', \
                  'UP-20171125210000-000009.nc', 'UP-20171125220000-000010.nc']
      files[1] = ['UP-20171125120000-000000.nc',                                \
                  'UP-20171125130000-000001.nc', 'UP-20171125140000-000002.nc', \
                  'UP-20171125150000-000003.nc', 'UP-20171125160000-000004.nc', \
                  'UP-20171125170000-000005.nc', 'UP-20171125180000-000006.nc', \
                  'UP-20171125190000-000007.nc', 'UP-20171125200000-000008.nc', \
                  'UP-20171125210000-000009.nc', 'UP-20171125220000-000010.nc']
      return files






   def CheckWallTime(self, s_path, t_path):
      isInt = False; isDyn = False; isPhy = False
      for i in range(len(self.c_int)):
         if s_path.count(self.c_int[i]) > 0:
            isInt = True
            iint  = i
      for i in range(len(self.c_dyn)):
         if s_path.count(self.c_dyn[i]+'_') > 0:
            isDyn = True
            idyn  = i
      for i in range(len(self.c_phy)):
         if s_path.count('_'+self.c_phy[i]) > 0:
            isPhy = True
            iphy  = i
      if not isInt or not isDyn or not isPhy: return
      
      if s_path.count('10h') > 0:
        i_int = 0
      elif s_path.count('02d') > 0:
        i_int = 1
      else:
        print 'check int index...'
        quit()


      src_file = s_path+'/'+'profile/ProfileInfoSummary.prf'
      trg_file = t_path+'/'+'profile/ProfileInfoSummary.prf'
      print src_file
      print trg_file

      if not os.path.isfile(src_file) or not os.path.isfile(trg_file):
         if not os.path.isfile(src_file):
            print 'file not found...', src_file
         if not os.path.isfile(trg_file):
            print 'file not found...', trg_file
         return


#      print '# START... ', self.c_int[iint], self.c_dyn[idyn], self.c_phy[iphy]
#      if self.c_opt == 0:
#         print ' - diff source: '+src_file
#         print ' - diff target: '+trg_file
#         command = 'diff '+src_file+' '+trg_file
#         os.system(command)
#         
#      elif self.c_opt == 1 or self.c_opt == 2 or self.c_opt == 3:
#         os.system('rm -rf org.txt new.txt')
#         if self.c_opt == 1: ncdump = 'ncdump -h '
#         if self.c_opt == 2: ncdump = 'ncdump '
#         if self.c_opt == 3: ncdump = 'ncdump -v '+self.c_var+' '
#         command = ncdump+src_file+' > '+'org.txt'
#         print command
#         os.system(command)
#         command = ncdump+trg_file+' > '+'new.txt'
#         print command
#         os.system(command)
#         command = 'diff org.txt new.txt | tee diff.txt'
#         print command
#         os.system(command)
#      else:
#         print 'check c option...', self.c_opt
#      print '# END  ... ', self.c_int[iint], self.c_dyn[idyn], self.c_phy[iphy]
#      print ' '

