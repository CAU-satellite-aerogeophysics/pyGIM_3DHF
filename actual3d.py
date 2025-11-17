"""
Georg HÃ¼ttner - September 2024

This script performs three dimensional thermal modelling from provided lithospheric
interfaces. Based on this a surface geothermal heat flow can be calculated.

Look on my works, ye Mighty, and despair!

There are a bunch of really stupid things in this script and the entire method needs improvement.
Some points to keep in mind: 
 - world size and extent of data should align with data resolution, as
   otherwise the gaps at the edges could lead "misapplied" markers
 - there are a lot of transposed 2d arrays in here, really make sure that your data is oriented correctly
 - mesh creation should probably be done with an external program, as the pygimli internal stuff
   takes unbearably long for large regions
 - marker application could probably be overhauled to allow for faster result (who would have thunk that
   a for loop over a couple of million cells takes long)
 - the solver itself is incredibly ram hungry and for larger models (cell count > 3 mil) a normal pc is
   gonna die
"""

import pygimli as pg
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from matplotlib import ticker
import cartopy.crs as ccrs
from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter)
from scipy.spatial import KDTree
import datetime
import pickle
import gepy

#%%___General parameters______________________________

#_____Coordinates_____________________________________

# like profile coord, but take x0,x1 and y0,y1 of area, area width/length should be divisible by 10000
area_coords = [1300000,1400000,-275000,-375000] # Lake Vostok "south-western" edge
# area_coord_s = [1355000,1435000,-855000,-985000] # Dome C
# area_coord_s = [1800000,1900000,-650000,-750000] # ASB

#_____Layers & thermal parameters_____________________

layers = 4
do_sediments = True
do_hp = True

k0 = 1.5
k1 = 2.7
k2 = 3.5

tc = [k0,k1,k2]

hp_0 = 0.000001
# If you have input HP data: do not use hp_0 in any function.
# If you do not have input HP data: set a hp_0 value and give it in 

#_____Data Input______________________________________

data_input = []

data_input.append(np.load("data_example/topography.npz"))
data_input.append(np.load("data_example/Sediments.npz"))
data_input.append(np.load("data_example/Moho.npz"))
data_input.append(np.load("data_example/LAB.npz"))
data_input.append(np.load("data_example/HP.npz"))

data_resolution = [1000,9250,10000,10000,12500]

#%%___Thermal parameters______________________________

k0 = 1.5
k1 = 2.7
k2 = 3.5
hp_0 = 0.000001

#%%___cut data to size & interpolate into region______

data_list = []

for i,layer in enumerate(data_input):
    data = gepy.cut_and_interpolate_3d(layer, data_resolution, i, layers,area_coords,
                                        do_sediments=do_sediments, do_hp=do_hp, hp_0=hp_0)
    data_list.append(data)

"""
#%% ja/nein/vielleicht?

notes = '3d_LV_A_1.5_norm'
print(notes)

grid_topo = np.full(np.shape(grid_topo),np.mean(grid_topo))
# grid_topo = grid_topo/1000
grid_sed = np.full(np.shape(grid_sed),np.mean(grid_sed))
grid_moho = np.full(np.shape(grid_moho),np.mean(grid_moho))
grid_lab = np.full(np.shape(grid_lab),np.mean(grid_lab))
# grid_A_for_topo = np.full(np.shape(grid_A_for_topo),np.mean(grid_A_for_topo))

grid_A_norm = grid_A_for_topo/np.max(grid_A_for_topo)
grid_A_for_topo = grid_A_for_topo*grid_A_norm*1.5
"""

#%%___build the world_________________________________

mesh,layer_list = gepy.build_world_3d(data_list, area_coords, area=2.5,
                                      layers = layers, do_sediments=do_sediments, border = 10)

# print(geometry)
# pg.show(geometry,showMesh=True,alpha=0.7)

# print(mesh)
# pg.show(mesh,showMesh=True,alpha=0.7)

#%%___apply correct markers to cells and nodes________



mesh,force = gepy.assign_markers_3d(data_list, mesh,
                         layers=layers, do_sediments=do_sediments, do_hp=do_hp, hp_0=hp_0, tc=tc)

# pg.show(mesh,showMesh=True)

#%% run calc and get shf

T = gepy.calc_temp(mesh, "3D", force, layers,do_sediments = True,tc = tc)

# pg.show(mesh,data=T,showMesh=True, label='Temperature in °C', cMap="inferno")

shf_3D = gepy.calc_ghf(T, mesh, data_list, "3D", tc)

# pg.show(mesh,data=gradient[:,2],filter={'clip':{'origin':(1050, 0, 0)},})


#%% Plot data and compare to 1D result

plt.figure(dpi=300)
plt.pcolormesh(x_int_t,y_int_t,shf_3D,vmin=40,vmax=80)
# plt.pcolormesh(x_int_t,y_int_t,shf_3D_norm,vmin=30,vmax=100)
# plt.pcolormesh(gradientTop_rs)
# plt.pcolormesh(x_int_t,y_int_t,geology_rs)
# plt.pcolormesh(x_int_t,y_int_t,grid_topo)
# plt.pcolormesh(x_int_t,y_int_t,grid_sed,vmin=0,vmax=70)
# plt.pcolormesh(x_int_HP,y_int_HP,np.transpose(grid_A))
# plt.pcolormesh(grid_A_for_topo)
plt.xlabel('x in m')
plt.ylabel('y in m')
plt.colorbar(label='mW/m$^2$')

grid_moho_t = grid_moho_for_topo - grid_sed + grid_topo
grid_lab_t = grid_lab_for_topo - grid_moho_for_topo

A_1D = grid_A_for_topo
# A_1D = 0.000001

shf_1D = (1315-0+(A_1D*grid_moho_t*grid_lab_t)/k2+(1/2*A_1D*grid_moho_t**2)/k1)/(grid_lab_t/k2+grid_moho_t/k1+grid_sed/k0)*1000

plt.figure(dpi=300)
plt.pcolormesh(x_int_t,y_int_t,shf_1D,vmin=40,vmax=80)#
plt.xlabel('x in m')
plt.ylabel('y in m')
plt.colorbar(label='mW/m$^2$')

plt.figure(dpi=300)
plt.pcolormesh(x_int_t,y_int_t,shf_3D-shf_1D,cmap='bwr',vmin=-15,vmax=15)
plt.xlabel('x in m')
plt.ylabel('y in m')
plt.colorbar(label='mW/m$^2$')


#%% save 

print(notes)
save_dict = {'shf_3D': shf_3D,
             'topcoord': topcoord,
             'T': T_list,
             'shf_1D': shf_1D,
             'grid_topo': grid_topo,
             'x_int_t': x_int_t,
             'y_int_t': y_int_t,
             'grid_sed': grid_sed,
             'grid_A': grid_A,
             'grid_A_for_topo': grid_A_for_topo,
             'x_int_HP': x_int_HP,
             'y_int_HP': y_int_HP,
             'grid_moho': grid_moho,
             'grid_moho_for_topo': grid_moho_for_topo,
             'x_int_m': x_int_m,
             'y_int_m': y_int_m,
             'grid_lab': grid_lab,
             'grid_lab_for_topo': grid_lab_for_topo,
             'x_int_l': x_int_l,
             'y_int_l': y_int_l
             }

now = datetime.datetime.now()
mesh.exportVTK('outputVTK/'+str(now.strftime('%Y-%m-%d_%H-%M-%S'))+'_'+str(notes))
geometry.exportPLC('outputVTK/'+str(now.strftime('%Y-%m-%d_%H-%M-%S'))+'_'+str(notes))
# geometry.exportSTL('outputVTK/'+str(now.strftime('%Y-%m-%d_%H-%M-%S'))+'_'+str(notes))
with open('output/'+str(now.strftime('%Y-%m-%d_%H-%M-%S'))+'_'+str(notes)+'.pkl','wb+') as file:
    pickle.dump(save_dict,file)


