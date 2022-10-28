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
output_folder='./RIDGECREST/'
mode='classic'

if mode=='classic':
    hobj.distance_calculation(Vp,Vs,stations)
    hloc=hades_location(hobj,output_folder)
    hloc.location(out_file,master=False,fixed=False)
elif mode=='master':
    hobj.relative_frame(Vp,Vs,stations,y_ref=-1,z_ref=-1,fixed_depth=1)
    hobj.distance_calculation(Vp,Vs,stations)
    hloc=hades_location(hobj,output_folder)
    hloc.location(out_file,master=False,fixed=False)
else:
    print('mode must be "classic" or "master"!!!')