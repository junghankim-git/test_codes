class TimeIntegration:
   def __init__(self, method, dt, nstep, uptype, ndims, arg, halo=0):
      self.ntl    = 3

      self.ntls   = 3
      self.nm1    = 0
      self.n0     = 1
      self.np1    = 2

      self.method = method
      self.dt     = dt
      self.istep  = 0
      self.nstep  = nstep
      self.uptype = uptype
      self.ndims  = ndims
      self.halo   = halo

      if method == 'Lagrangian':
         self.nparts  = arg
      elif method == 'Eulerian':
         self.dimsize = arg
      else:
         quit()

   def UpdateLevel(self, lup=True):
      tmp      = self.nm1
      self.nm1 = self.n0
      self.n0  = self.np1
      self.np1 = tmp
      if lup == True:
         self.istep = self.istep + 1


   def ForwardEuler(self, psi_io, RHS_i, dt):
      method = self.method
      if method == 'Lagrangian':
         nparts = self.nparts
         ndims  = self.ndims
         for ipart in range(nparts):
            for idim in range(ndims):
               psi_io[self.np1][ipart][idim] = psi_io[self.n0][ipart][idim]+dt*RHS_i[ipart][idim]
      elif method == 'Eulerian':
         ndims   = self.ndims
         dimsize = self.dimsize
         halo    = self.halo
         if ndims == 1:
            for i in range(halo,dimsize+halo):
               psi_io[self.np1][i] = psi_io[self.n0][i]+dt*RHS_i[i]
         elif ndims == 2:
            for i in range(halo,dimsize[0]+halo):
               for j in range(halo,dimsize[1]+halo):
                  psi_io[self.np1][i][j] = psi_io[self.n0][i][j]+dt*RHS_i[i][j]
         elif ndims == 3:
            for i in range(halo,dimsize[0]+halo):
               for j in range(halodimsize[1]+halo):
                  for k in range(halodimsize[2]+halo):
                     psi_io[self.np1][i][j][k] = psi_io[self.n0][i][j][k]+dt*RHS_i[i][j][k]

   def Leapfrog(self, psi_io, RHS_i, dt):
      method = self.method
      if method == 'Lagrangian':
         nparts = self.nparts
         ndims  = self.ndims
         for ipart in range(nparts):
            for idim in range(ndims):
               psi_io[self.np1][ipart][idim] = psi_io[self.nm1][ipart][idim]+2.0*dt*RHS_i[ipart][idim]
      elif method == 'Eulerian':
         ndims   = self.ndims
         dimsize = self.dimsize
         halo    = self.halo
         if ndims == 1:
            for i in range(halo,dimsize+halo):
               psi_io[self.np1][i] = psi_io[self.nm1][i]+2.0*dt*RHS_i[i]
         elif ndims == 2:
            for i in range(halo,dimsize[0]+halo):
               for j in range(halo,dimsize[1]+halo):
                  psi_io[self.np1][i][j] = psi_io[self.nm1][i][j]+2.0*dt*RHS_i[i][j]
         elif ndims == 3:
            for i in range(halo,dimsize[0]+halo):
               for j in range(halo,dimsize[1]+halo):
                  for k in range(halo,dimsize[2]+halo):
                     psi_io[self.np1][i][j][k] = psi_io[self.nm1][i][j][k]+2.0*dt*RHS_i[i][j][k]

   def Update(self, psi_io, RHS_i, dt=None):
      dt = self.dt if dt == None else dt

      if self.uptype== 'ForwardEuler':
         self.ForwardEuler(psi_io, RHS_i, dt)
      elif self.uptype == 'Leapfrog':
         if self.istep == 0:
            self.ForwardEuler(psi_io, RHS_i, dt)
         else:
            self.Leapfrog(psi_io, RHS_i, dt)
      elif self.uptype == 'RK2':
         self.ForwardEuler(psi_io, RHS_i, 0.5*dt)
         #self.UpdateLevel(False)
         self.UpdateLevel() #why?
         self.Leapfrog(psi_io, RHS_i, dt)


