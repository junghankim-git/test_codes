#!/usr/bin/env python




def Numerical(dt_in, x0, v0, a):
   x1 = x0 + dt_in*v0
   v1 = v0 + dt_in*a
   return x1, v1

def Analytic(t_in, x0, v0, a):
   x_out = x0 + v0*t_in + 0.5*a*t_in*t_in
   v_out = v0 + a*t_in
   return x_out, v_out


nstep = 1000

h0 = 2.0*9.8
g  = -9.8

ft = 2.0
dt = ft/nstep

h = [0.0 for i in range(nstep+1)]
v = [0.0 for i in range(nstep+1)]

h[0] = h0
v[0] = 0.0


t = 0.0
xa, va = Analytic(t, h[0], v[0], g)
print '%000007.3f: | %000007.3f, %000007.3f | %000007.3f, %000007.3f | %000007.3f'%(t, xa, va, h[0], v[0], xa-h[0])
for istep in range(nstep):
   h[istep+1], v[istep+1] = Numerical(dt, h[istep], v[istep], g)
   t = dt*(istep+1)
   xa, va = Analytic(t, h[0], v[0], g)
   print '%000007.3f: | %000007.3f, %000007.3f | %000007.3f, %000007.3f | %000007.3f'%(t, xa, va, h[istep+1], v[istep+1], xa-h[istep+1])
