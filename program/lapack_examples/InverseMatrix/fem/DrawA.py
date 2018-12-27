from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
#ax = fig.gca(projection='3d')
ax = fig.add_subplot(1,1,1)
X = np.arange(-5, 5, 1.0)
Y = np.arange(-5, 5, 1.0)
X, Y = np.meshgrid(X, Y)
#R = np.sqrt(X**2 + Y**2)
#Z = np.sin(R)

'''
Z = [[0.100000e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00], \
 [0.100000e+01, 0.111111e+00, 0.123457e-01, 0.137174e-02, 0.152416e-03, 0.169351e-04, 0.188168e-05, 0.209075e-06, 0.232306e-07, 0.258117e-08], \
 [0.100000e+01, 0.222222e+00, 0.493827e-01, 0.109739e-01, 0.243865e-02, 0.541923e-03, 0.120427e-03, 0.267616e-04, 0.594703e-05, 0.132156e-05], \
 [0.100000e+01, 0.333333e+00, 0.111111e+00, 0.370370e-01, 0.123457e-01, 0.411523e-02, 0.137174e-02, 0.457247e-03, 0.152416e-03, 0.508053e-04], \
 [0.100000e+01, 0.444444e+00, 0.197531e+00, 0.877915e-01, 0.390184e-01, 0.173415e-01, 0.770735e-02, 0.342549e-02, 0.152244e-02, 0.676639e-03], \
 [0.100000e+01, 0.555556e+00, 0.308642e+00, 0.171468e+00, 0.952599e-01, 0.529221e-01, 0.294012e-01, 0.163340e-01, 0.907444e-02, 0.504136e-02], \
 [0.100000e+01, 0.666667e+00, 0.444444e+00, 0.296296e+00, 0.197531e+00, 0.131687e+00, 0.877915e-01, 0.585277e-01, 0.390184e-01, 0.260123e-01], \
 [0.100000e+01, 0.777778e+00, 0.604938e+00, 0.470508e+00, 0.365950e+00, 0.284628e+00, 0.221377e+00, 0.172182e+00, 0.133920e+00, 0.104160e+00], \
 [0.100000e+01, 0.888889e+00, 0.790123e+00, 0.702332e+00, 0.624295e+00, 0.554929e+00, 0.493270e+00, 0.438462e+00, 0.389744e+00, 0.346439e+00], \
 [0.100000e+01, 0.100000e+01, 0.100000e+01, 0.100000e+01, 0.100000e+01, 0.100000e+01, 0.100000e+01, 0.100000e+01, 0.100000e+01, 0.100000e+01]]
'''

Z = [[0.100000e+01, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,-0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00], \
 [-0.254607e+02, 0.810000e+02,-0.162000e+03, 0.252000e+03,-0.283500e+03, 0.226800e+03,-0.126000e+03, 0.462857e+02,-0.101250e+02, 0.100000e+01], \
 [ 0.261763e+03,-0.133332e+04, 0.339564e+04,-0.566010e+04, 0.658024e+04,-0.536625e+04, 0.301905e+04,-0.111896e+04, 0.246399e+03,-0.244607e+02], \
 [-0.145382e+04, 0.920297e+04,-0.271253e+05, 0.489841e+05,-0.594044e+05, 0.497087e+05,-0.284536e+05, 0.106772e+05,-0.237316e+04, 0.237303e+03], \
 [ 0.486949e+04,-0.349328e+05, 0.113455e+06,-0.219411e+06, 0.278499e+06,-0.240251e+06, 0.140501e+06,-0.535633e+05, 0.120501e+05,-0.121652e+04], \
 [-0.102960e+05, 0.800339e+05,-0.278309e+06, 0.568880e+06,-0.753879e+06, 0.671949e+06,-0.402804e+06, 0.156521e+06,-0.357472e+05, 0.365297e+04], \
 [ 0.138396e+05,-0.113669e+06, 0.415557e+06,-0.887949e+06, 0.122268e+07,-0.112562e+07, 0.693088e+06,-0.275316e+06, 0.640313e+05,-0.664301e+04], \
 [-0.114671e+05, 0.979844e+05,-0.372009e+06, 0.823734e+06,-0.117249e+07, 0.111270e+07,-0.704159e+06, 0.286599e+06,-0.680909e+05, 0.719660e+04], \
 [ 0.533814e+04,-0.469756e+05, 0.183632e+06,-0.418510e+06, 0.612818e+06,-0.597871e+06, 0.388616e+06,-0.162279e+06, 0.395022e+05,-0.427051e+04], \
 [-0.106763e+04, 0.960864e+04,-0.384346e+05, 0.896807e+05,-0.134521e+06, 0.134521e+06,-0.896807e+05, 0.384346e+05,-0.960864e+04, 0.106763e+04]]



#surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
surf = ax.imshow(Z, cmap=cm.coolwarm)
#ax.set_ylim(0.00, 1.00)

#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()

outfilename = 'A.png'
fig.savefig(outfilename)