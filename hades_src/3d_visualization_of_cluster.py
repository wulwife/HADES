import numpy as num
import matplotlib.pyplot as plt

path='/home/francesco/Seistools/HADES/hades_src/RIDGECREST'

filename='/hades_ridgecrestcartesian.txt'
f=open(path+filename,'r')
f.readline()
f.readline()
events=[]
for line in f:
   toks=line.split()
   events.append([eval(toks[1]),eval(toks[2]),eval(toks[3])])
events=num.array(events)
f.close()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(events[:,0],events[:,1],events[:,2])
ax.set_xlim([-30,30])
ax.set_ylim([-30,30])
ax.set_zlim([-10,10])
plt.show()