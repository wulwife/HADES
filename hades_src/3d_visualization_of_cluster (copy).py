import numpy as num
import matplotlib.pyplot as plt

path='/home/francesco/Seistools/HADES/hades_src/RIDGECREST'

filename1='/hades_ridgecrestcartesian_relative.txt'
f=open(path+filename1,'r')
f.readline()
f.readline()
events=[]
for line in f:
   toks=line.split()
   events.append([eval(toks[1]),eval(toks[2]),eval(toks[3])])
events=num.array(events)
f.close()

filename2='/hades_ridgecrestcartesian.txt'
f=open(path+filename2,'r')
f.readline()
f.readline()
events2=[]
for line in f:
   toks=line.split()
   events2.append([eval(toks[1]),eval(toks[2]),eval(toks[3])])
events2=num.array(events2)
f.close()

print(events[5,:],events2[5,:])
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(events[:,0],events[:,1],events[:,2])
ax.scatter(events2[:,0],events2[:,1],events2[:,2])
ax.set_xlim([-30,30])
ax.set_ylim([-30,30])
ax.set_zlim([5,15])
plt.show()