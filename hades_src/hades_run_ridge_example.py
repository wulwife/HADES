from hades_input import hades_input
from hades_location import hades_location

data_path='./'
input_file='ridgecrest_gji.dat'
sta_file='stations_ridge.txt'
out_file='hades_ridgecrest'
stations=['WBS','WMF']
Vp=6000
Vs=Vp/1.73
hobj=hades_input(data_path,input_file,sta_file)
hobj.distance_calculation(Vp,Vs,stations)
hloc=hades_location(hobj,'./NAPA/')
hloc.location(out_file)
