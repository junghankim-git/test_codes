load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "/home/jhkim/Study/Library/Shared/NCL/Lib_NCL.ncl"

begin

a = new((/4/), double)
b = new((/4/), double)

a(0) = 1.0
a(1) = 2.0
a(2) = 3.0
a(3) = 4.0

b(0) = 1.1
b(1) = 2.1
b(2) = 3.1
b(3) = 4.1


print("a = "+a)
print("b = "+b)

ave = avg(a)
dev = stddev(a)

print("ave = "+ave)
print("dev = "+dev)

output = AveDev(a)
print("output = "+output)
print(" ")

L1Norm = GetL1Norm(a,b)
L2Norm = GetL2Norm(a,b)
LiNorm = GetLinfNorm(a,b)
print("L1   norm = "+L1Norm)
print("L2   norm = "+L2Norm)
print("Linf norm = "+LiNorm)
L1Norm = GetLNorm(a,b,"L1")
L2Norm = GetLNorm(a,b,"L2")
LiNorm = GetLNorm(a,b,"Linf")
print("L1   norm = "+L1Norm)
print("L2   norm = "+L2Norm)
print("Linf norm = "+LiNorm)

L1Error = GetL1Error(a,b)
L2Error = GetL2Error(a,b)
LiError = GetLinfError(a,b)
print("L1   err. = "+L1Error)
print("L2   err. = "+L2Error)
print("Linf err. = "+LiError)

L1Error = GetLError(a,b,"L1")
L2Error = GetLError(a,b,"L2")
LiError = GetLError(a,b,"Linf")
print("L1   err. = "+L1Error)
print("L2   err. = "+L2Error)
print("Linf err. = "+LiError)


c = new((/2,2/),double)
d = new((/2,2/),double)
c(0,0) = 1.0
c(0,1) = 2.0
c(1,0) = 3.0
c(1,1) = 4.0
d(0,0) = 1.1
d(0,1) = 2.1
d(1,0) = 3.1
d(1,1) = 4.1
print(" ")
L1Norm = GetLNorm2D(c,d,"L1")
L2Norm = GetLNorm2D(c,d,"L2")
LiNorm = GetLNorm2D(c,d,"Linf")
print("L1   norm = "+L1Norm)
print("L2   norm = "+L2Norm)
print("Linf norm = "+LiNorm)
L1Error = GetLError2D(c,d,"L1")
L2Error = GetLError2D(c,d,"L2")
LiError = GetLError2D(c,d,"Linf")
print("L1   err. = "+L1Error)
print("L2   err. = "+L2Error)
print("Linf err. = "+LiError)


end
