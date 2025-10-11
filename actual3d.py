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
import pygimli.meshtools as mt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from matplotlib import ticker
import cartopy.crs as ccrs
from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter)
#import polygon_setup as ps
from scipy.interpolate import griddata,LinearNDInterpolator, RegularGridInterpolator
from scipy.spatial import KDTree
import datetime
import pickle
from scipy.stats import norm
import gepy
from pyproj import CRS, Transformer
km = 1000
# from skimage.util.shape import view_as_windows

# data_calc,data_resolution = gepy.load_data()

#%%  Parameters

#_____Coordinates_____

# like profile coord, but take x0,x1 and y0,y1 of area, area width/length should be divisible by 10000
area_coords = [1300000,1400000,-275000,-375000] # Lake Vostok "south-western" edge
# area_coord_s = [1355000,1435000,-855000,-985000] # Dome C
# area_coord_s = [1800000,1900000,-650000,-750000] # ASB

#_____Layers & thermal parameters_____

layers = 4
do_sediments = True
do_hp = True

k0 = 1.5
k1 = 2.7
k2 = 3.5
hp_0 = 0.000001

#_____Data Input____

data_input = []

data_input.append(np.load("data_example/topography.npz"))
data_input.append(np.load("data_example/Sediments.npz"))
data_input.append(np.load("data_example/Moho.npz"))
data_input.append(np.load("data_example/LAB.npz"))
data_input.append(np.load("data_example/HP.npz"))

"""
data_list = []
for i in data:
    data_list.append(i["data"],i["x"],i["y"])
"""

data_resolution = [1000,9250,10000,10000,12500]

#%% cut data to size & interpolate into region


data_list = []

for i,layer in enumerate(data_input):
    data = gepy.load_and_interpolate_3d(layer, data_resolution, i, layers,area_coords,
                                        do_sediments=do_sediments, do_hp=do_hp, hp_0=hp_0)

"""
gepy.load_and_interpolate(topo["data"],topo["x"],topo["y"],
                          sed["data"],sed["x"],sed["y"],
                          moho["data"],moho["x"],moho["y"],
                          LAB["data"],LAB["x"],LAB["y"],
                          HP["data"],HP["x"],HP["y"],
                          resolution,
                          area_coords)
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

#%% build the world

# create 3d mesh from the input data

# the +/- 10 are to create a slight buffer around the input data and world edge, there are issues
# in the mesh creation when edges of boundaries touch other boundaries
border = 10
start = [(area_coords[0]-border)/km,(area_coords[2]+border)/km,4] 
end = [(area_coords[1]+border)/km,(area_coords[3]-border)/km,-220]

world = mt.createWorld(start=start,end=end, worldMarker=False)

# create layer polygons and shift them to the correct height
mesh_moho = mt.createMesh2D(x_int_m/km,y_int_m/km)
mesh_lab  = mt.createMesh2D(x_int_l/km,y_int_l/km)

surface_moho = mt.createSurface(mesh_moho)
surface_lab  = mt.createSurface(mesh_lab)

surface_moho = gepy.fix_surface_height(surface_moho, x_int_m/km, y_int_m/km, -grid_moho/km)
surface_lab  = gepy.fix_surface_height(surface_lab, x_int_l/km, y_int_l/km, -grid_lab/km)

for boundary in surface_moho.boundaries():
    boundary.setMarker(6)
for boundary in surface_lab.boundaries():
    boundary.setMarker(7)

surface_topo,surface_sed = gepy.create_sed_interface(x_int_t, y_int_t, grid_topo, grid_sed)

geometry = world + surface_topo + surface_sed + surface_moho + surface_lab

mesh = mt.createMesh(geometry,quality=34,area=2.5)#,area=2.5

# print(geometry)
# pg.show(geometry,showMesh=True,alpha=0.7)

# print(mesh)
# pg.show(mesh,showMesh=True,alpha=0.7)

#%% apply correct markers to cells and nodes, and create temp and force vectors

grid_xy_ts = np.array(list(zip(xi_int_t.flatten()/km,yi_int_t.flatten()/km)))
tree_ts = KDTree(grid_xy_ts)
grid_xy_m = np.array(list(zip(xi_int_m.flatten()/km,yi_int_m.flatten()/km)))
tree_m = KDTree(grid_xy_m)
grid_xy_l = np.array(list(zip(xi_int_l.flatten()/km,yi_int_l.flatten()/km)))
tree_l = KDTree(grid_xy_l)

# create lists for temperature and heat production that are to be applied to the nodes

force = np.zeros(mesh.nodeCount())
# Tnode = []
for i, node in enumerate(mesh.nodes()):
    x,y,z = node.pos()
    z_topo = gepy.compare_to_plane(x,y,grid_topo,grid_xy_ts,tree_ts)
    z_sed = gepy.compare_to_plane(x,y,(grid_topo-grid_sed),grid_xy_ts,tree_ts)
    z_moho = gepy.compare_to_plane(x,y,-grid_moho,grid_xy_m,tree_m)
    z_lab = gepy.compare_to_plane(x,y,-grid_lab,grid_xy_l,tree_l)
    
    idx_A,idy_A = np.argmin(np.abs(x_int_t/km-x)),np.argmin(np.abs(y_int_t/km-y))
    
    if z >= z_moho and z <= z_sed:
        force[i] = grid_A_for_topo[idy_A,idx_A]
        # force[i] = 0.000001
    else:
        force[i] = 0
    
# apply correct marker to the cells
for i,cell in enumerate(mesh.cells()):
    center = cell.center()
    x, y, z = center.x(), center.y(), center.z()
    
    z_topo = gepy.compare_to_plane(x,y,grid_topo,grid_xy_ts,tree_ts)
    z_sed = gepy.compare_to_plane(x,y,(grid_topo-grid_sed),grid_xy_ts,tree_ts)
    z_moho = gepy.compare_to_plane(x,y,-grid_moho,grid_xy_m,tree_m)
    z_lab = gepy.compare_to_plane(x,y,-grid_lab,grid_xy_l,tree_l)
    
    if z < z_topo:
        cell.setMarker(4)
    if z < z_sed:
        cell.setMarker(5)
    if z < z_moho:
        cell.setMarker(6)
    if z < z_lab:
        cell.setMarker(7)

# pg.show(mesh,showMesh=True)

#%% run calc and get shf

T = pg.solver.solveFiniteElements(mesh,
                                a={1: 1.0*km, 4: 1.5*km, 5: 2.7*km, 6: 3.5*km, 7: 4.0*km},
                                f=force*km*km*km,
                                bc={'Dirichlet': {7: 1315, 4: 0}},verbose=True)#{'Node': Tnode}

T_list = [T[i] for i in range(len(T))]

# pg.show(mesh,data=T,showMesh=True, label='Temperature in °C', cMap="inferno")

# gradient = pg.solver.grad(mesh, T)

# pg.show(mesh,data=gradient[:,2],filter={'clip':{'origin':(1050, 0, 0)},})

topcoord = []
for i in range(len(x_int_t)):
    for j in range(len(y_int_t)):
        topcoord.append([x_int_t[i]/km,y_int_t[j]/km,(grid_topo[j,i]-5)/km])

gradientTop = pg.solver.grad(mesh, T, topcoord)

geology_list = []

shf = np.zeros(len(gradientTop))
for i in range(len(gradientTop)):
    point = pg.RVector3(topcoord[i][0],topcoord[i][1],topcoord[i][2])
    cell = mesh.findCell(point)
    geology = cell.marker()
    geology_list.append(geology)
    if geology == 5:   
        shf[i] = np.linalg.norm(gradientTop[i,:])*2.7
    elif geology == 4:
        shf[i] = np.linalg.norm(gradientTop[i,:])*1.5

shf_3D = np.transpose(np.reshape(shf,(len(x_int_t),len(y_int_t))))
gradientTop_rs = np.reshape(gradientTop[:,2],(len(x_int_t),len(y_int_t)))

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


