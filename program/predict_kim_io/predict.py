#!/usr/bin/env python


# # of components (ini, dyn, phy, io, tot)
ncs = 5

base_proc = 10008.

wall_2 = [160., 14595., 2588., 1637., 0.]
wall_3 = [130., 14841., 3030., 1256., 0.]
wall_4 = [124., 15266., 3449., 1128., 0.]
wall_2[4] = wall_2[0]+wall_2[1]+wall_2[2]+wall_2[3]
wall_3[4] = wall_3[0]+wall_3[1]+wall_3[2]+wall_3[3]
wall_4[4] = wall_4[0]+wall_4[1]+wall_4[2]+wall_4[3]
#wall_2[4] = sum(wall_2[0:3])
#wall_3[4] = sum(wall_3[0:3])
#wall_4[4] = sum(wall_4[0:3])

def io_alpha():
    '''
    a1 = 1920.
    a2 = 19200.
    b1 = 68.
    b2 = 403.
    '''
    a1 = 15360.
    b1 = 319.
    a2 = 30720.
    b2 = 714.
    alpha = (b2-b1)/(a2-a1)
    beta  = b2-(alpha*a2)
    return alpha

def io_wall(proc_base, wall_base, proc):
    return wall_base-io_alpha()*(proc_base-proc)


def get_walls(wall_base, proc_base, proc):
    result = [0. for i in range(ncs)]
    frac = proc/proc_base
    #result[0] = wall_base[0]*frac
    result[0] = io_wall(proc_base, wall_base[0], proc)
    result[1] = wall_base[1]/frac
    result[2] = wall_base[2]/frac
    #result[3] = wall_base[3]*frac
    result[3] = io_wall(proc_base, wall_base[3], proc)
    result[4] = result[0]+result[1]+result[2]+result[3]
    io_ratio  = (result[0]+result[3])/result[4]*100.
    return result, io_ratio


nps   = 12
fracs = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5]
procs = [0. for i in range(nps)]
for i in range(nps):
  procs[i] = fracs[i]*base_proc

'''
print ' 1920 ', io_wall(1920.)
print ' 7680 ', io_wall(7680.)
print '19200 ', io_wall(19200.)
print '24000 ', io_wall(24000.)
print '30720 ', io_wall(30720.)
'''

print(procs)
print(wall_2)
print(wall_3)
print(wall_4)

walls_2t = [[0. for i in range(ncs)] for j in range(nps)]
walls_3t = [[0. for i in range(ncs)] for j in range(nps)]
walls_4t = [[0. for i in range(ncs)] for j in range(nps)]

for i in range(nps):
    walls_2t[i], io_2t = get_walls(wall_2, base_proc, procs[i])
    walls_3t[i], io_3t = get_walls(wall_3, base_proc, procs[i])
    walls_4t[i], io_4t = get_walls(wall_4, base_proc, procs[i])
    print('# nprocs = {:05d}'.format(int(procs[i])))
    #print(' - thread 2: dyn={:7.1f}, phy={:7.1f}, io={:7.1f}'.format(walls_2t[i][1], walls_2t[i][2], walls_2t[i][0]+walls_2t[i][3]))
    #print(' - thread 3: dyn={:7.1f}, phy={:7.1f}, io={:7.1f}'.format(walls_3t[i][1], walls_3t[i][2], walls_3t[i][0]+walls_3t[i][3]))
    #print(' - thread 4: dyn={:7.1f}, phy={:7.1f}, io={:7.1f}'.format(walls_4t[i][1], walls_4t[i][2], walls_4t[i][0]+walls_4t[i][3]))
    print(' - elapse(2t)={:7.1f} ({:4.1f}%), elapse(3t)={:7.1f} ({:4.1f}%), elapse(4t)={:7.1f} ({:4.1f}%)'. \
          format(walls_2t[i][-1], io_2t, walls_3t[i][-1], io_3t, walls_4t[i][-1], io_4t))





