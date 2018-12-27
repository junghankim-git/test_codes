class Mesh:
   
   def __init__(self):
      # direction
      self.w  = 0
      self.e  = 1
      self.s  = 2
      self.n  = 3
      self.ws = 4
      self.es = 5 # changed
      self.wn = 6 # changed
      self.en = 7
      self.dirstr = ['w', 'e', 's', 'n']
      # [axis=0: x-axis], [axis=1: y-axis], [axis=2: x.y-axis], [axis=3: (-x).y or x.(-y)-axis], [axis=4: (0,0)]
      self.x   = 0
      self.y   = 1
      self.xy  = 2
      self.mxy = 3
      self.org = 4


   def GetDirectionIndex(self):
      return self.w, self.e, self.s, self.n, self.ws, self.wn, self.es, self.en
   

   def GetAxisIndex(self):
      return self.x, self.y, self.xy, self.mxy, self.org
   

   def SymmetryMesh(self, n, axis, mesh):
   
      result = [[0 for i in range(n)] for j in range(n)]
   
      for i in range(n):
         for j in range(n):
            if axis==self.x:
              result[i][j] = mesh[i][n-j-1]
            elif axis==self.y:
              result[i][j] = mesh[n-i-1][j]
            elif axis==self.xy:
              result[i][j] = mesh[j][i]
            elif axis==self.mxy:
              result[i][j] = mesh[n-j-1][n-i-1]
            elif axis==self.org:
              result[i][j] = mesh[n-i-1][n-j-1]
            else:
              print '## -- Not support a axis value: '+str(axis)
      return result
   
   
   def SymmetryDirMesh(self, n, axis, mesh):
      #   ir=0 (w,left), il=1 (e,right), id=2 (s,down), iu=3(n,up)
      tmp = self.SymmetryMesh(n, axis, mesh)
      result = [[0 for i in range(n)] for j in range(n)]
      for i in range(n):
         for j in range(n):
            result[i][j] = tmp[i][j]
            if axis==self.x:
               if tmp[i][j] == self.s: result[i][j] = self.n
               if tmp[i][j] == self.n: result[i][j] = self.s
            if axis==self.y:
               if tmp[i][j] == self.w: result[i][j] = self.e
               if tmp[i][j] == self.e: result[i][j] = self.w
            if axis==self.xy:
               if tmp[i][j] == self.w: result[i][j] = self.s
               if tmp[i][j] == self.e: result[i][j] = self.n
               if tmp[i][j] == self.s: result[i][j] = self.w
               if tmp[i][j] == self.n: result[i][j] = self.e
            if axis==self.mxy:
               if tmp[i][j] == self.w: result[i][j] = self.n
               if tmp[i][j] == self.e: result[i][j] = self.s
               if tmp[i][j] == self.s: result[i][j] = self.e
               if tmp[i][j] == self.n: result[i][j] = self.w
            if axis==self.org:
               if tmp[i][j] == self.w: result[i][j] = self.e
               if tmp[i][j] == self.e: result[i][j] = self.w
               if tmp[i][j] == self.s: result[i][j] = self.n
               if tmp[i][j] == self.n: result[i][j] = self.s
      return result
              
   
   
   def VectorToMatrix2D(self, nx_in, ny_in, vector_in):
      matrix_out = [[[0.0 for j in range(ny_in)] for i in range(nx_in)] for idim in range(2)]
      for i in range(nx_in):
         for j in range(ny_in):
            for idim in range(2):
               matrix_out[idim][i][j] = vector_in[i][j][idim]
      return matrix_out


   
   
   def PrintMesh2D(self, nx, ny, mesh, message='', dir=False, lmatrix=False):
      print ' '
      if message != '':
         print '# Mesh: '+message
   
   
      if lmatrix == False:
         for j in range(ny):
            string = ' j = '+str(ny-j)+': ['
            for i in range(nx):
               if dir:
                  if mesh[i][ny-j-1] < 0:
                     dirstr = 'x'
                  else:
                     dirstr = str(self.dirstr[mesh[i][ny-j-1]])
                  string = string + dirstr + ', '
               else:
                  string = string + str(mesh[i][ny-j-1]) + ', '
            string = string[0:-2]+']'
            print string
         print '    i =  1, 2, ...'
      else:
         for i in range(nx):
            string = ' i = '+str(i)+': ['
            for j in range(ny):
               string = string + str(mesh[i][j]) + ', '
            string = string[0:-2]+']'
            print string
         print '    j =  1, 2, ...'
      print ' '
