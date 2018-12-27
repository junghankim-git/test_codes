#!/usr/bin/env python

p0 = 100000.0
kapa = 2.8571658640413355e-1
def center(p1,p2):
    p12_org = (p1+p2)/2.
    pi1 = (p1/p0)**kapa
    pi2 = (p2/p0)**kapa
    pi12 = (pi1+pi2)/2.
    p12 = p0*pi12**(1./kapa)
    return p12_org, p12

def exner(p):
    return (p/p0)**kapa

pi1 = 900000.
pi2 = 800000.
pim = (pi1+pi2)/2.0
ex1 = (pi1/p0)**kapa
ex2 = (pi2/p0)**kapa
exm = (pim/p0)**kapa
print pi1, pim, pi2
print ex1, exm, ex2
print (ex1+ex2)/2.0
print center(pi1,pi2)


p = [100000.0,50000.0,30000.0,10000.0]

for i in range(len(p)):
    print p[i], exner(p[i])
