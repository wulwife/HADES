from hades_input import hades_input
from hades_location import hades_location

data_path='./'
input_file='napa_gji.dat'
sta_file='stations_napa.txt'
out_file='hades_napa'
stations='CVS'
Vp=6000
Vs=Vp/1.73
hobj=hades_input(data_path,input_file,sta_file)
hobj.distance_calculation(Vp,Vs,stations)
hloc=hades_location(hobj,'./NAPA/')
hloc.location(out_file)
