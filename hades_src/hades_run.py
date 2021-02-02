from hades_input import hades_input
from hades_location import hades_location

data_path='./'
input_file='napa_gji.dat'#'ridgecrest_gji.dat'
sta_file='stations_napa.txt'#'stations_ridge.txt'
out_file='hades_napa'#'hades_ridgecrest'
stations='CVS'#['WBS','WMF']
Vp=6000
Vs=Vp/1.73
hobj=hades_input(data_path,input_file,sta_file)
hobj.distance_calculation(Vp,Vs,stations)
hloc=hades_location(hobj,'./')
hloc.location(out_file)
