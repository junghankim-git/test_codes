

def check_dimension(var):
  ndim = 0

  if isinstance(var, list):
    ndim = 1
  else:
    return ndim

  if isinstance(var[0], list):
    ndim = 2
  else:
    return ndim

  if isinstance(var[0][0], list):
    ndim = 3
  else:
    return ndim

  if isinstance(var[0][0][0], list):
    ndim = 4
  else:
    return ndim


def get_lnorm(var, pnorm):
#  import math

  if pnorm != 'L1' and pnorm != 'L2' and pnorm != 'Linf':
    print 'check pnorm...'
    return -99

  ndims = check_dimension(var)
  if ndims > 3:
    print 'check dimensions of var...'
    return -99

  norm  = 0.0

  if ndims == 0:
    norm = abs(var)
  elif ndims == 1:
    nx = len(var)
    if pnorm == 'L1':
      norm = 1.0*sum([abs(var[i]) for i in range(nx)])
    elif pnorm == 'L2':
      norm = 1.0*(sum([var[i]**2 for i in range(nx)]))**0.5
    elif pnorm == 'Linf':
      norm = 1.0*max([abs(var[i]) for i in range(nx)])
  elif ndims == 2:
    nx = len(var)
    ny = len(var[0])
    if pnorm == 'L1':
      norm = 1.0*sum([abs(var[i][j]) for j in range(ny) for i in range(nx)])
    elif pnorm == 'L2':
      norm = 1.0*(sum([var[i][j]**2 for j in range(ny) for i in range(nx)]))**0.5
    elif pnorm == 'Linf':
      norm = 1.0*max([abs(var[i][j]) for j in range(ny) for i in range(nx)])
  elif ndims == 3:
    nx = len(var)
    ny = len(var[0])
    nz = len(var[0][0])
    if pnorm == 'L1':
      norm = 1.0*sum([abs(var[i][j][k]) for k in range(nz) for j in range(ny) for i in range(nx)])
    elif pnorm == 'L2':
      norm = 1.0*(sum([var[i][j][k]**2 for k in range(nz) for j in range(ny) for i in range(nx)]))**0.5
    elif pnorm == 'Linf':
      norm = 1.0*max([abs(var[i][j][k]) for k in range(nz) for j in range(ny) for i in range(nx)])

  return norm




def get_lerror(var1, var2, pnorm):
#  import math

  if pnorm != 'L1' and pnorm != 'L2' and pnorm != 'Linf':
    print 'check pnorm...'
    return -99

  ndims1 = check_dimension(var1)
  ndims2 = check_dimension(var2)
  if ndims1 != ndims2:
    print 'check dimensions of vars...'
    return -99
  if ndims1 > 3:
    print 'check dimensions of var...'
    return -99

  error = 0.0
  if ndims1 == 0:
    diff = var2 - var1
  elif ndims1 == 1:
    nx   = len(var1)
    diff = [var2[i]-var1[i] for i in range(nx)]
  elif ndims1 == 2:
    nx   = len(var1)
    ny   = len(var1[0])
    diff = [[var2[i][j]-var1[i][j] for j in range(ny)] for i in range(nx)]
  elif ndims1 == 3:
    nx   = len(var1)
    ny   = len(var1[0])
    nz   = len(var1[0][0])
    diff = [[[var2[i][j][k]-var1[i][j][k] for k in range(nz)] for j in range(ny)] for i in range(nx)]


  basenorm = get_lnorm(var1, pnorm)
  diffnorm = get_lnorm(diff, pnorm)

  if basenorm == 0:
    error = 0.0
  else:
    error = diffnorm/basenorm

  return error


