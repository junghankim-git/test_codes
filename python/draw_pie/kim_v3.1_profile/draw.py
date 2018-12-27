#!/usr/bin/env python

import matplotlib.pyplot as plt

# total 14960.713
# dyn:  11279.746
# pys:  1607.902
# I/O:   826.717

#print 826.717/14960.713
ncomp = 4
total = 14960.713
walltime = [11279.746,1607.902,826.717,14960.713-11279.746-1607.902-826.717]
sizes = [walltime[i]/total*100. for i in range(ncomp)]

#labels = 'Dynamics','Physics','Write','Others (initialize, read, ...)'
labels = ['Dynamics','Physics','Write','Others']
#sizes = [15,30,45,10]
explode = (0,0,0,0)  # only "explode" the 2nd slice (i.e. 'Hogs')

'''
fig1,ax1 = plt.subplots()
ax1.pie(sizes,explode=explode,labels=labels,autopct='%1.1f%%',shadow=True,startangle=90)
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
#ax1.set_clip_on(False)
plt.show()

'''
pltfig = plt.figure(figsize=(8,7))
#pltfig.subplots_adjust(left=0,right=1,bottom=0,top=1)
axis   = pltfig.add_subplot(1,1,0)

#axis.pie(sizes,explode=explode,labels=labels,colors=['lightcoral','c','g','y'],  \
#         autopct='%1.1f%%',shadow=0,startangle=90)#,labeldistance=0.7)
#patche,text,autotext = axis.pie(sizes,explode=explode,labels=labels,colors=['lightcoral','c','g','y'],  \
patche,text,autotext = axis.pie(sizes,explode=explode,colors=['lightcoral','c','g','y'],  \
         autopct='%1.1f%%',shadow=0,startangle=90,labeldistance=0.60,pctdistance=0.5)
for t in text:
    t.set_size(20)
for t in autotext:
    t.set_size(0)
    #t.set_color('y')
axis.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
pltfig.savefig('pie.png')

plt.show()
