from hades_input import hades_input
from hades_location import hades_location
from numpy import sqrt

data_path='/Users/francesco/Codes/HADES_NEW/Datasets/Synth'
input_file='Synth.dat'
sta_file='Stations_synth.txt'
out_file='hades_joint'
stations=['STA1','STA2']
Vp=6000
Vs=Vp/sqrt(3)
hobj=hades_input(data_path,input_file,sta_file)
output_folder='./SYNTH/'
fixed_dist=True
joint=True

hobj.distance_calculation(Vp,Vs,stations,refcat='tsp')
hloc=hades_location(hobj,output_folder)
hloc.location(out_file,joint,fixed_dist)