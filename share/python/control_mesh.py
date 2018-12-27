

# direction
ww_ = 0
ee_ = 1
ss_ = 2
nn_ = 3
ws_ = 4
es_ = 5 # changed
wn_ = 6 # changed
en_ = 7
dirstr_ = ['w', 'e', 's', 'n']
# [axis=0: x-axis], [axis=1: y-axis], [axis=2: x.y-axis], [axis=3: (-x).y or x.(-y)-axis], [axis=4: (0,0)]
xx_  = 0
yy_  = 1
xy_  = 2
mxy_ = 3
org_ = 4


def get_direction_index():
   return ww_, ee_, ss_, nn_, ws_, wn_, es_, en_


def get_axis_index():
   return xx_, yy_, xy_, mxy_, org_


def symmetry_mesh(n, axis, mesh):

   result = [[0 for i in range(n)] for j in range(n)]

   for i in range(n):
      for j in range(n):
         if axis==xx_:
           result[i][j] = mesh[i][n-j-1]
         elif axis==yy_:
           result[i][j] = mesh[n-i-1][j]
         elif axis==xy_:
           result[i][j] = mesh[j][i]
         elif axis==mxy_:
           result[i][j] = mesh[n-j-1][n-i-1]
         elif axis==org_:
           result[i][j] = mesh[n-i-1][n-j-1]
         else:
           print '## -- Not support a axis value: '+str(axis)
   return result


def symmetry_dir_mesh(n, axis, mesh):
   #   ir=0 (w,left), il=1 (e,right), id=2 (s,down), iu=3(n,up)
   tmp = symmetry_mesh(n, axis, mesh)
   result = [[0 for i in range(n)] for j in range(n)]
   for i in range(n):
      for j in range(n):
         result[i][j] = tmp[i][j]
         if axis==xx_:
            if tmp[i][j] == ss_: result[i][j] = nn_
            if tmp[i][j] == nn_: result[i][j] = ss_
         if axis==yy_:
            if tmp[i][j] == ww_: result[i][j] = ee_
            if tmp[i][j] == ee_: result[i][j] = ww_
         if axis==xy_:
            if tmp[i][j] == ww_: result[i][j] = ss_
            if tmp[i][j] == ee_: result[i][j] = nn_
            if tmp[i][j] == ss_: result[i][j] = ww_
            if tmp[i][j] == nn_: result[i][j] = ee_
         if axis==mxy_:
            if tmp[i][j] == ww_: result[i][j] = nn_
            if tmp[i][j] == ee_: result[i][j] = ss_
            if tmp[i][j] == ss_: result[i][j] = ee_
            if tmp[i][j] == nn_: result[i][j] = ww_
         if axis==org_:
            if tmp[i][j] == ww_: result[i][j] = ee_
            if tmp[i][j] == ee_: result[i][j] = ww_
            if tmp[i][j] == ss_: result[i][j] = nn_
            if tmp[i][j] == nn_: result[i][j] = ss_
   return result
           


def vector_to_matrix_2d(nx_in, ny_in, vector_in):
   matrix_out = [[[0.0 for j in range(ny_in)] for i in range(nx_in)] for idim in range(2)]
   for i in range(nx_in):
      for j in range(ny_in):
         for idim in range(2):
            matrix_out[idim][i][j] = vector_in[i][j][idim]
   return matrix_out




def print_mesh_2d(nx, ny, mesh, message='', dir=False, lmatrix=False):
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
                  dirstr = str(dirstr_[mesh[i][ny-j-1]])
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
