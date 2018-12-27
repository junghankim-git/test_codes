from numpy import *
class KIAPSGM_Prof:

   def __init__(self, filename):
      print 'filename = '+filename
      self.filename  = filename
      self.nprocs    = 0
      self.nthreads  = 0
      self.maxlevels = 10
      self.nlevel    = 0

      infile         = open(filename, 'r')
      self.Parse(infile)
      self.Statistic()
      infile.close()
      self.SetProcThread(-1,-1)
      print ' '
      print '#'*63
      print '{0:s} {1:5d} {2:s} {3:3d} {4:20s}'.format('# Total: nprocs = ', self.nprocs, ',     nthreads = ', self.nthreads, ' '*15+'#')
      print '#'*63
      print ' '
       


   # public
   def SetProcThread(self, iproc, ithrd):
      if iproc < 0 and ithrd < 0:
        self.profname_out  = self.tprofname
        self.proflevel_out = self.tproflevel
        self.walltime_out  = self.twalltime
        self.proftime_out  = self.tproftime
        self.walltimeD_out = self.twalltimeD
        self.proftimeD_out = self.tproftimeD
        self.walltimeM_out = self.twalltimeM
        self.proftimeM_out = self.tproftimeM
      elif iproc >= 0 and ithrd < 0:
        self.profname_out  = self.tprofname
        self.proflevel_out = self.tproflevel
        self.walltime_out  = self.pwalltime[iproc]
        self.proftime_out  = self.pproftime[iproc]
        self.walltimeD_out = self.pwalltimeD[iproc]
        self.proftimeD_out = self.pproftimeD[iproc]
      else:
        self.profname_out  = self.profname[iproc][ithrd]
        self.proflevel_out = self.proflevel[iproc][ithrd]
        self.walltime_out  = self.walltime[iproc][ithrd]
        self.proftime_out  = self.proftime[iproc][ithrd]
        self.walltimeD_out = self.walltimeD[iproc][ithrd]
        self.proftimeD_out = self.proftimeD[iproc][ithrd]


   # public
   def ViewProfLevelTree(self, ilev_in=-1):
      print ' ### START:: Profiling Levels ### '

      nlevels = [0 for ilev in range(self.maxlevels)]
      nitems  = self.nitems
      bilev   = -1

      print '-'*96
      print '{0:80s} {1:15.5f} '.format(' (0)0.: Total', self.walltime_out)
      print '-'*96
      for item in range(nitems):
         ilev = self.proflevel_out[item]
         if ilev > bilev:
            for i in range(ilev,self.maxlevels): nlevels[i]=0
         nlevels[ilev] = nlevels[ilev] + 1
         levstr  = ' ('+str(ilev)+')'
         for i in range(1,ilev+1):
            levstr = levstr+str(nlevels[i])+'.'
         levstr = levstr+' '*ilev
         if ilev_in != -1:
            if ilev == ilev_in:
               print '{0:80s} {1:15.5f} '.format(levstr+': '+self.profname_out[item], self.proftime_out[item])
         else:
            print '{0:80s} {1:15.5f} '.format(levstr+': '+self.profname_out[item], self.proftime_out[item])
         bilev = ilev
      print '-'*96
      print ' ###  END :: Profiling Levels ### \n\n'



   # public
   def GetProfInfo(self, ilev, pos, isDev=False):
      # lev 0: Total 
      # lev 1: IniGrid, AtmosModelSet, AtmosModelRun, ...
      if ilev == 0:
         if isDev:
            return 1, 'Total', self.walltime_out, self.walltimeD_out
         else:
            return 1, 'Total', self.walltime_out
      else:
         name   = []
         wtime  = []
         wtimeD = []
         wtimeM = []

         nitems  = self.nitems

         n_uplev = [0 for i in range(ilev-1)]
         onoff   = [False for i in range(ilev-1)]
         for iitem in range(nitems):
            clev = self.proflevel_out[iitem]
            if clev < ilev:
               for ulev in range(ilev-1):
                  if ulev == clev-1:
                     if ulev > 0:
                        if onoff[ulev-1] == True:
                           n_uplev[ulev] = n_uplev[ulev] + 1
                     elif ulev == 0:
                           n_uplev[ulev] = n_uplev[ulev] + 1
                  if n_uplev[ulev] == pos[ulev]:
                     onoff[ulev] = True
                  else:
                     onoff[ulev] = False
            #print '\n'
            #print iitem+1, ilev, n_uplev[:], pos[:], onoff[:]

            fonoff = True
            for ulev in range(ilev-1):
               if onoff[ulev] == False: fonoff = False
   
            if fonoff == True:
               if self.proflevel_out[iitem] == ilev:
                  name.append(self.profname_out[iitem])
                  wtime.append(self.proftime_out[iitem])
                  wtimeD.append(self.proftimeD_out[iitem])
      if isDev:
         return len(wtime), array(name), array(wtime), array(wtimeD)
      else:
         return len(wtime), array(name), array(wtime)


   # public
   def GetProfInfo_Name(self, name, isLev=False, isDev=False, isMax=False):
      if name == 'Total':
         if isLev:
            if isDev==True and isMax==False:
               return self.walltime_out, self.walltimeD_out, 1
            if isDev==True and isMax==True:
               return self.walltime_out, self.walltimeD_out, self.walltimeM_out, 1
            else:
               return self.walltime_out, 1
         else:
            if isDev==True and isMax==False:
               return self.walltime_out, self.walltimeD_out
            if isDev==True and isMax==True:
               return self.walltime_out, self.walltimeD_out, self.walltimeM_out
            else:
               return self.walltime_out
      else:
         nitems  = self.nitems
         isVal = False
         for iitem in range(nitems):
            if self.profname_out[iitem] == name:
               isVal = True
               if isLev:
                  if isDev==True and isMax==False:
                     return self.proftime_out[iitem], self.proftimeD_out[iitem], self.proflevel_out[iitem]
                  if isDev==True and isMax==True:
                     return self.proftime_out[iitem], self.proftimeD_out[iitem], self.proftimeM_out[iitem], self.proflevel_out[iitem]
                  else:
                     return self.proftime_out[iitem], self.proflevel_out[iitem]
                  break
               else:
                  if isDev==True and isMax==False:
                     return self.proftime_out[iitem], self.proftimeD_out[iitem]
                  if isDev==True and isMax==True:
                     return self.proftime_out[iitem], self.proftimeD_out[iitem], self.proftimeM_out[iitem]
                  else:
                     return self.proftime_out[iitem]
                  break
         if isVal == False:
            print 'cannot find variable:', name
            quit()


   # private
   def Parse(self, file):
      self.nprocs, self.nthreads = self.GetNumberOfProcsThreads(file)
      # walltime [nprocs][nthreads]
      # profname [nprocs][nthreads][nitems]
      # proflevel[nprocs][nthreads][nitems]
      # proftime [nprocs][nthreads][nitems]
      self.walltime   = [[0.0 for ithrd in range(self.nthreads)] for iprocs in range(self.nprocs)]
      self.profname   = [[[] for ithrd in range(self.nthreads)] for iprocs in range(self.nprocs)]
      self.proflevel  = [[[] for ithrd in range(self.nthreads)] for iprocs in range(self.nprocs)]
      self.proftime   = [[[] for ithrd in range(self.nthreads)] for iprocs in range(self.nprocs)]

      for iproc in range(self.nprocs):
         for ithrd in range(self.nthreads):
            i, j = self.GetProcThreadNumber(file)
            self.walltime[iproc][ithrd] = self.GetWallTime(file)
            self.GetLevelNameTime(file, iproc, ithrd)



   # private
   def Statistic(self):
   
      self.nitems     = len(self.proflevel[0][0])
      self.tprofname  = self.profname[0][0]
      self.tproflevel = self.proflevel[0][0]
      # tprofname [nitems]
      # tproflevel[nitems]

      # Average
      self.twalltime  = 0.0
      self.tproftime  = [0.0 for iitem in range(self.nitems)]

      self.pwalltime  = [0.0 for iproc in range(self.nprocs)]
      self.pproftime  = [[0.0 for iitem in range(self.nitems)] for iproc in range(self.nprocs)]


      # Standard Deviation
      self.twalltimeD = 0.0
      self.tproftimeD = [0.0 for iitem in range(self.nitems)]

      self.pwalltimeD = [0.0 for iproc in range(self.nprocs)]
      self.pproftimeD = [[0.0 for iitem in range(self.nitems)] for iproc in range(self.nprocs)]

      # Maximium
      self.twalltimeM = 0.0
      self.tproftimeM = [0.0 for iitem in range(self.nitems)]

      # Get Average
      for iproc in range(self.nprocs):
         for ithrd in range(self.nthreads):
            self.twalltime        = self.twalltime + self.walltime[iproc][ithrd]
            self.pwalltime[iproc] = self.pwalltime[iproc]+self.walltime[iproc][ithrd]
         self.pwalltime[iproc]    = self.pwalltime[iproc]/float(self.nthreads)
      self.twalltime = self.twalltime/float(self.nprocs*self.nthreads)

      for item in range(self.nitems):
         for iproc in range(self.nprocs):
            for ithrd in range(self.nthreads):
               self.tproftime[item] = self.tproftime[item] + self.proftime[iproc][ithrd][item]
               self.pproftime[iproc][item] = self.pproftime[iproc][item] + self.proftime[iproc][ithrd][item]
            self.pproftime[iproc][item] = self.pproftime[iproc][item]/float(self.nthreads)
         self.tproftime[item] = self.tproftime[item]/float(self.nprocs*self.nthreads)


      # Get Standard Deviation
      for iproc in range(self.nprocs):
         for ithrd in range(self.nthreads):
            self.twalltimeD        = self.twalltimeD+(self.twalltime - self.walltime[iproc][ithrd])**2.0
            self.pwalltimeD[iproc] = self.pwalltimeD[iproc]+(self.walltime[iproc][ithrd]-self.pwalltime[iproc])**2.0
         self.pwalltimeD[iproc]    = (self.pwalltimeD[iproc]/float(self.nthreads))**0.5
      self.twalltimeD = (self.twalltimeD/float(self.nprocs*self.nthreads))**0.5

      for item in range(self.nitems):
         for iproc in range(self.nprocs):
            for ithrd in range(self.nthreads):
               self.tproftimeD[item] = self.tproftimeD[item] + (self.proftime[iproc][ithrd][item]-self.tproftime[item])**2.0
         self.tproftimeD[item] = (self.tproftimeD[item]/float(self.nprocs*self.nthreads))**0.5

      # Get Maximium
      maxval = 0.0
      for iproc in range(self.nprocs):
         for ithrd in range(self.nthreads):
            self.twalltimeM       = max(self.walltime[iproc][ithrd], maxval)
            maxval = self.twalltimeM

      for item in range(self.nitems):
         maxval = 0.0
         for iproc in range(self.nprocs):
            for ithrd in range(self.nthreads):
               self.tproftimeM[item] = max(self.proftime[iproc][ithrd][item], maxval)
               maxval = self.tproftimeM[item]


   # private
   def GetNumberOfProcsThreads(self, file):
      # fine # of processes and # of threads
      nprocs   = 0
      nthreads = 0
      for iline in file:
         line   = iline.strip()
         splitline = str(line).split()
         if len(splitline) > 5 and splitline[0] == 'Total':
           nprocs   = int(splitline[1])
           nthreads = int(splitline[2])/nprocs
           break
      return nprocs, nthreads


   # private
   def GetProcThreadNumber(self, file):
      iproc = -1
      ithrd = -1
      for iline in file:
         line   = iline.strip()
         splitline = str(line).split()
         if len(splitline) > 5 and splitline[1] == 'PROCESS':
            iproc = splitline[2]
         if len(splitline) > 3 and splitline[0] == 'Stats' and splitline[2] == 'thread':
            ithrd = splitline[3][0]
            break
      return iproc, ithrd


   # private
   def GetWallTime(self, file):
      for iline in file:
         line   = iline.strip()
         splitline = str(line).split()
         if iline[0:7] == '  Total':
            return float(splitline[4])
            break


   # private
   def GetLevelNameTime(self, file, ip, it):
      for iline in file:
         if self.IsThreadEnd(iline):
            break
         else:
            ilev, isstar = self.IsLevel(self.maxlevels, iline)
            (self.proflevel[ip][it]).append(ilev)
            line = iline.strip()
            splitline = str(line).split()
            if isstar:
               (self.profname[ip][it]).append(splitline[1])
               (self.proftime[ip][it]).append(float(splitline[5]))
            else:
               (self.profname[ip][it]).append(splitline[0])
               (self.proftime[ip][it]).append(float(splitline[4]))

   # private
   def IsThreadEnd(self, line):
      if len(line) > 10 and line[0:7] == '  t_prf':
         return True
      else:
         return False


   # private
   def IsLevel(self, maxlev, line):
      for i in range(maxlev):
         ilev = maxlev - i
         if len(line) > 10 and (line[0:2*(ilev+1)] == (ilev+1)*'  ' or line[0:2*(ilev+1)] == '* '+ilev*'  '):
            if line[0:2*(ilev+1)] == (ilev+1)*'  ':
              return ilev, False
              break
            elif line[0:2*(ilev+1)] == '* '+ilev*'  ':
              return ilev, True
              break


