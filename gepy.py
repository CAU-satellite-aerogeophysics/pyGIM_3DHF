"""
Georg HÃ¼ttner - October 2024

One of probably a bunch of scripts to import dumbass functions
"""

import pygimli as pg
import pygimli.meshtools as mt
import numpy as np
import pandas as pd
from scipy.interpolate import griddata, RegularGridInterpolator
from pyproj import Proj

km = 1000

# def load_and_interpolate_2d()
# def load_and_interpolate_3d()


# class build_world ? 
#   would maybe allow for things like setting up differen thermal parameters for the same mesh
# def build_world()


# def assign t_params

# def calc_ghf()

def build_world(data_list           ,
                area_coords         ,
                layers = None       ,
                do_sediments = False,
                do_hp = False       ,
                border = None       
                ):
    
    start = [(area_coords[0]-border)/km,(area_coords[2]+border)/km,4] 
    end = [(area_coords[1]+border)/km,(area_coords[3]-border)/km,-220]
    
    world = mt.createWorld(start=start,end=end, worldMarker=False)
    
    # create layer polygons and shift them to the correct height
    
    if layers == 2:
        # topo + LAB
    elif layers == 2 and do_hp:
        # topo + LAB + HP
    elif layers == 3 and not do_sediments and not do_hp:    
        # topo + Moho + LAB
    elif layers == 3 and not do_sediments and do_hp:
        # topo + Moho + LAB + HP
    elif layers == 3 and do_sediments and not do_hp:
        # topo + Sed + LAB
    elif layers == 3 and do_sediments and do_hp:
        # topo + Sed + LAB + HP
    elif layers == 4 and do_sediments and not do_hp:
        # topo + Sed + Moho + LAB
    elif layers == 4 and do_sediments and do_hp:   
        # topo + Sed + Moho + LAB + HP
        

    mesh_moho = mt.createMesh2D(x_int_m/km,y_int_m/km)
    mesh_lab  = mt.createMesh2D(x_int_l/km,y_int_l/km)
    
return mesh
    
    
    
"""
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
"""




















"""
def load_and_interpolate_3d():
    xi_topo_s, yi_topo_s, grid_topo_s = gepy.in_area_s(area_coord_s,data_calc[0][0],data_calc[0][1],data_calc[0][2])
    xi_sed_s, yi_sed_s, grid_sed_s    = gepy.in_area_s(area_coord_s,data_calc[1][0],data_calc[1][1],data_calc[1][2])
    xi_moho_s, yi_moho_s, grid_moho_s = gepy.in_area_s(area_coord_s,data_calc[2][0],data_calc[2][1],data_calc[2][2])
    xi_lab_s, yi_lab_s, grid_lab_s    = gepy.in_area_s(area_coord_s,data_calc[3][0],data_calc[3][1],data_calc[3][2])

    xii_topo_s,yii_topo_s = np.meshgrid(xi_topo_s,yi_topo_s)
    xii_sed_s,yii_sed_s   = np.meshgrid(xi_sed_s,yi_sed_s)
    xii_moho_s,yii_moho_s = np.meshgrid(xi_moho_s,yi_moho_s)
    xii_lab_s,yii_lab_s   = np.meshgrid(xi_lab_s,yi_lab_s)

    # 

    res_topo,res_moho,res_lab = (1000,10000,10000)

    x_int_t,y_int_t = np.arange(area_coord_s[0],area_coord_s[1]+res_topo,res_topo),np.arange(area_coord_s[3],area_coord_s[2]+res_topo,res_topo)
    x_int_m,y_int_m = np.arange(area_coord_s[0],area_coord_s[1]+res_moho,res_moho),np.arange(area_coord_s[3],area_coord_s[2]+res_moho,res_moho)
    x_int_l,y_int_l = np.arange(area_coord_s[0],area_coord_s[1]+res_lab,res_lab),np.arange(area_coord_s[3],area_coord_s[2]+res_lab,res_lab)

    xi_int_t,yi_int_t = np.meshgrid(x_int_t,y_int_t)
    xi_int_m,yi_int_m = np.meshgrid(x_int_m,y_int_m)
    xi_int_l,yi_int_l = np.meshgrid(x_int_l,y_int_l)

    interp_topo = RegularGridInterpolator((np.unique(xi_topo_s),np.unique(yi_topo_s)),grid_topo_s)
    grid_topo = interp_topo((xi_int_t,yi_int_t))

    interp_sed = RegularGridInterpolator((np.unique(xi_sed_s),np.unique(yi_sed_s)),grid_sed_s)
    grid_sed = interp_sed((xi_int_t,yi_int_t),method='linear')

    for i in range(np.shape(grid_sed)[0]):
        for j in range(np.shape(grid_sed)[1]):
            grid_sed[i,j] = round(grid_sed[i,j]/35)*35

    interp_moho = RegularGridInterpolator((np.unique(xi_moho_s),np.unique(yi_moho_s)),grid_moho_s)
    grid_moho = interp_moho((xi_int_m,yi_int_m))
    grid_moho_for_topo = interp_moho((xi_int_t,yi_int_t))

    interp_lab = RegularGridInterpolator((np.unique(xi_lab_s),np.unique(yi_lab_s)),grid_lab_s)
    grid_lab = interp_lab((xi_int_l,yi_int_l))
    grid_lab_for_topo = interp_lab((xi_int_t,yi_int_t))

    # heat production

    xi_smos_s, yi_smos_s, A_s = gepy.in_area_s(area_coord_s,xi_smos,yi_smos,A)
    xii_smos_s,yii_smos_s   = np.meshgrid(xi_smos_s,yi_smos_s)
        
    res_HP = 12500
        
    x_int_HP,y_int_HP = np.arange(area_coord_s[0],area_coord_s[1]+res_HP,res_HP),np.arange(area_coord_s[3],area_coord_s[2]+res_HP,res_HP)
    xi_int_HP,yi_int_HP = np.meshgrid(x_int_HP,y_int_HP)

    grid_A  = griddata((xii_smos_s.flatten(), yii_smos_s.flatten()), np.transpose(A_s).flatten(), (xi_int_HP,yi_int_HP),fill_value=0.000001)
    interp_A = RegularGridInterpolator((x_int_HP,y_int_HP),np.transpose(grid_A))
    grid_A_for_topo = interp_A((xi_int_t,yi_int_t))
     
    return 

    del A_s,xi_smos,yi_smos,xi_smos_s,yi_smos_s,xii_smos_s,yii_smos_s    
    del grid_topo_s,xi_topo_s,yi_topo_s,xii_topo_s,yii_topo_s
    del grid_sed_s,xi_sed_s,yi_sed_s,xii_sed_s,yii_sed_s
    del grid_moho_s,xi_moho_s,yi_moho_s,xii_moho_s,yii_moho_s
    del grid_lab_s,xi_lab_s,yi_lab_s,xii_lab_s,yii_lab_s

"""

def cut_and_interpolate_3d(layer               ,
                           resolution          ,
                           i                   ,
                           layers              ,
                           area_coords         ,
                           do_sediments = False,
                           do_hp = False       ,
                           hp_0 = None
                           ):
    
    data,xi,yi = layer["data"],layer["x"],layer["y"]
    
    xi_s,yi_s,data_s = in_area_s(area_coords,xi,yi,data)
    xii_s,yii_s = np.meshgrid(xi_s,yi_s)
    x_int,y_int = np.arange(area_coords[0],area_coords[1]+resolution[i],resolution[i]),np.arange(area_coords[3],area_coords[2]+resolution[i],resolution[i])
    xi_int,yi_int = np.meshgrid(x_int,y_int)

    if do_hp and layers == i:
        grid_temp  = griddata((xii_s.flatten(), yii_s.flatten()), data.flatten(), (xi_int,yi_int),fill_value=hp_0)
        interp = RegularGridInterpolator((x_int,y_int),grid_temp)
        x_int_t,y_int_t = np.arange(area_coords[0],area_coords[1]+resolution[0],resolution[0]),np.arange(area_coords[3],area_coords[2]+resolution[0],resolution[0])
        xi_int_t,yi_int_t = np.meshgrid(x_int_t,y_int_t)
        grid = interp((xi_int_t,yi_int_t))
    else:
        interp = RegularGridInterpolator((np.unique(xi_s),np.unique(yi_s)),data_s)
        grid = interp((xi_int,yi_int))

    if do_sediments and i==1:
        for i in range(np.shape(grid)[0]):
            for j in range(np.shape(grid)[1]):
                grid[i,j] = round(grid[i,j]/35)*35

    return x_int,y_int,grid

"""
data_calc,data_resolution = gepy.load_data()

x_profile,y_profile = (profile_coord[0],profile_coord[1])

profile_length = np.round(np.sqrt((x_profile[1]-x_profile[0])**2 + (y_profile[1]-y_profile[0])**2))

x_int_t, y_int_t, dist_prof_t = gepy.define_profile(x_profile, y_profile, data_resolution[0])
x_int_s, y_int_s, dist_prof_s = gepy.define_profile(x_profile, y_profile, data_resolution[1])
x_int_m, y_int_m, dist_prof_m = gepy.define_profile(x_profile, y_profile, data_resolution[2])
x_int_l, y_int_l, dist_prof_l = gepy.define_profile(x_profile, y_profile, data_resolution[3])

dist_prof_t = np.round(dist_prof_t)
dist_prof_s = np.round(dist_prof_s)
dist_prof_m = np.round(dist_prof_m)
dist_prof_l = np.round(dist_prof_l)

xvals_t = np.unique(np.asarray(data_calc[0][0]))
yvals_t = np.unique(np.asarray(data_calc[0][1]))

xvals_s = np.unique(np.asarray(data_calc[1][0]))
yvals_s = np.unique(np.asarray(data_calc[1][1]))

xvals_m = np.unique(np.asarray(data_calc[2][0]))
yvals_m = np.unique(np.asarray(data_calc[2][1]))

xvals_l = np.unique(np.asarray(data_calc[3][0]))
yvals_l = np.unique(np.asarray(data_calc[3][1]))

interp_func = RegularGridInterpolator((yvals_t,xvals_t),np.asarray(data_calc[0][2]))
profile_topo = interp_func((y_int_t,x_int_t),method="linear").tolist()
profile_topo = np.round(profile_topo,decimals=3)

interp_func = RegularGridInterpolator((yvals_s,xvals_s),np.asarray(data_calc[1][2]))
profile_sed = interp_func((y_int_s,x_int_s),method="linear").tolist()
profile_sed = np.round(profile_sed,decimals=3)

interp_func = RegularGridInterpolator((yvals_m,xvals_m),np.asarray(data_calc[2][2]))
profile_moho = interp_func((y_int_m,x_int_m),method="linear").tolist()
profile_moho = np.round(profile_moho,decimals=3)

interp_func = RegularGridInterpolator((yvals_l,xvals_l),np.asarray(data_calc[3][2]))
profile_lab = interp_func((y_int_l,x_int_l),method="linear").tolist()
profile_lab = np.round(profile_lab,decimals=3)

profile_topo_for_sed = np.interp(dist_prof_s,dist_prof_t,profile_topo)
profile_sed_for_topo = np.interp(dist_prof_t,dist_prof_s,profile_sed)
profile_sed_for_topo = np.round(profile_sed_for_topo,decimals=-1)

profile_moho_for_topo = np.interp(dist_prof_t,dist_prof_m,profile_moho)
profile_lab_for_topo = np.interp(dist_prof_t,dist_prof_l,profile_lab)

# heat production:
data_HP = np.load("data/SMOS_hp.npz")
grid_A_ws = data_HP['A_ws']
grid_A_wo = data_HP['A_wo']

xi_smos,yi_smos = data_HP['xi_smos'],data_HP['yi_smos']
x_int_A, y_int_A, dist_prof_A = gepy.define_profile(x_profile, y_profile, 12500)
dist_prof_A = np.round(dist_prof_A)
xvals_A = np.unique(xi_smos)
yvals_A = np.unique(yi_smos)
interp_func = RegularGridInterpolator((yvals_A,xvals_A),np.asarray(grid_A_ws))
profile_A = interp_func((y_int_A,x_int_A),method="linear").tolist()
# profile_A = np.round(profile_A,decimals=3)
profile_A_for_topo = np.interp(dist_prof_t,dist_prof_A,profile_A)

# profile_A_norm = profile_A_for_topo/np.max(profile_A_for_topo)
# profile_A_for_topo = profile_A_for_topo*profile_A_norm*1.5


# single layer scheiße

xi_topo_s_a2, yi_topo_s_a2, grid_topo_s_a2 = gepy.in_area_s(area_coord_s,data_calc[0][0],data_calc[0][1],data_calc[0][2])
xi_sed_s_a2, yi_sed_s_a2, grid_sed_s_a2    = gepy.in_area_s(area_coord_s,data_calc[1][0],data_calc[1][1],data_calc[1][2])
xi_moho_s_a2, yi_moho_s_a2, grid_moho_s_a2 = gepy.in_area_s(area_coord_s,data_calc[2][0],data_calc[2][1],data_calc[2][2])
xi_lab_s_a2, yi_lab_s_a2, grid_lab_s_a2    = gepy.in_area_s(area_coord_s,data_calc[3][0],data_calc[3][1],data_calc[3][2])
xii_topo_s_a2,yii_topo_s_a2 = np.meshgrid(xi_topo_s_a2,yi_topo_s_a2)
xii_sed_s_a2,yii_sed_s_a2   = np.meshgrid(xi_sed_s_a2,yi_sed_s_a2)
xii_moho_s_a2,yii_moho_s_a2 = np.meshgrid(xi_moho_s_a2,yi_moho_s_a2)
xii_lab_s_a2,yii_lab_s_a2   = np.meshgrid(xi_lab_s_a2,yi_lab_s_a2)
res_topo,res_moho,res_lab = (1000,10000,10000)
x_int_t_a2,y_int_t_a2 = np.arange(area_coord_s[0],area_coord_s[1]+res_topo,res_topo),np.arange(area_coord_s[3],area_coord_s[2]+res_topo,res_topo)
x_int_m_a2,y_int_m_a2 = np.arange(area_coord_s[0],area_coord_s[1]+res_moho,res_moho),np.arange(area_coord_s[3],area_coord_s[2]+res_moho,res_moho)
x_int_l_a2,y_int_l_a2 = np.arange(area_coord_s[0],area_coord_s[1]+res_lab,res_lab),np.arange(area_coord_s[3],area_coord_s[2]+res_lab,res_lab)
xi_int_t_a2,yi_int_t_a2 = np.meshgrid(x_int_t_a2,y_int_t_a2)
xi_int_m_a2,yi_int_m_a2 = np.meshgrid(x_int_m_a2,y_int_m_a2)
xi_int_l_a2,yi_int_l_a2 = np.meshgrid(x_int_l_a2,y_int_l_a2)
interp_topo_a2 = RegularGridInterpolator((np.unique(xi_topo_s_a2),np.unique(yi_topo_s_a2)),grid_topo_s_a2)
grid_topo_a2 = interp_topo_a2((xi_int_t_a2,yi_int_t_a2))
interp_sed_a2 = RegularGridInterpolator((np.unique(xi_sed_s_a2),np.unique(yi_sed_s_a2)),grid_sed_s_a2)
grid_sed_a2 = interp_sed_a2((xi_int_t_a2,yi_int_t_a2),method='linear')
interp_moho_a2 = RegularGridInterpolator((np.unique(xi_moho_s_a2),np.unique(yi_moho_s_a2)),grid_moho_s_a2)
grid_moho_a2 = interp_moho_a2((xi_int_m_a2,yi_int_m_a2))
grid_moho_for_topo_a2 = interp_moho_a2((xi_int_t_a2,yi_int_t_a2))
interp_lab_a2 = RegularGridInterpolator((np.unique(xi_lab_s_a2),np.unique(yi_lab_s_a2)),grid_lab_s_a2)
grid_lab_a2 = interp_lab_a2((xi_int_l_a2,yi_int_l_a2))
grid_lab_for_topo_a2 = interp_lab_a2((xi_int_t_a2,yi_int_t_a2))
xi_smos_s_a2, yi_smos_s_a2, A_s_a2 = gepy.in_area_s(area_coord_s,xi_smos,yi_smos,grid_A_ws)
xii_smos_s_a2,yii_smos_s_a2   = np.meshgrid(xi_smos_s_a2,yi_smos_s_a2)    
res_HP = 12500
x_int_HP_a2,y_int_HP_a2 = np.arange(area_coord_s[0],area_coord_s[1]+res_HP,res_HP),  np.arange(area_coord_s[3],area_coord_s[2]+res_HP,res_HP)
xi_int_HP_a2,yi_int_HP_a2 = np.meshgrid(x_int_HP_a2,y_int_HP_a2)
grid_A_a2  = griddata((xii_smos_s_a2.flatten(), yii_smos_s_a2.flatten()), np.transpose(A_s_a2).flatten(), (xi_int_HP_a2,yi_int_HP_a2),fill_value=0.000001)
# interp_A = RegularGridInterpolator((x_int_HP,y_int_HP),np.transpose(grid_A))
interp_A_a2 = RegularGridInterpolator((x_int_HP_a2,y_int_HP_a2),np.transpose(grid_A_a2))
grid_A_for_topo_a2 = interp_A_a2((xi_int_t_a2,yi_int_t_a2))

# select which interfaces are part of the world

# profile_topo = np.full(len(profile_topo),np.mean(grid_topo_a2))
# profile_sed = np.full(len(profile_sed),np.mean(grid_sed_a2))
# profile_sed_for_topo = np.full(len(profile_sed_for_topo),np.mean(grid_sed_a2))
# profile_moho = np.full(len(profile_moho),np.mean(grid_moho_a2))
# profile_moho_for_topo = np.full(len(profile_topo),np.mean(grid_moho_for_topo_a2))
# profile_lab = np.full(len(profile_lab),np.mean(grid_lab_a2))
# profile_lab_for_topo = np.full(len(profile_topo),np.mean(grid_lab_for_topo_a2))
# profile_A = np.full(len(profile_A),np.mean(grid_A_a2))
# profile_A_for_topo = np.full(len(profile_A_for_topo),np.mean(grid_A_for_topo_a2))

profile_sed_depth = profile_topo-profile_sed_for_topo
"""

def cut_and_interpolate_2d(layer               ,
                           resolution          ,
                           i                   ,
                           layers              ,
                           area_coords         ,
                           profile_coord       ,
                           do_sediments = False,
                           do_hp = False       ,
                           hp_0 = None
                           ):
    
    data,xi,yi = layer["data"],layer["x"],layer["y"]
    
    x_profile,y_profile = (profile_coord[0],profile_coord[1])
    x_int, y_int, dist_prof = define_profile(x_profile, y_profile, resolution[i])
    dist_prof = np.round(dist_prof)
    
    xvals = np.unique(np.asarray(xi))
    yvals = np.unique(np.asarray(yi))
    
    interp_func = RegularGridInterpolator((yvals,xvals),np.asarray(data))
    profile = interp_func((y_int,x_int),method="linear").tolist()
    profile = np.round(profile,decimals=3)
    
    if do_sediments and i==1:
        x_int_t, y_int_t, dist_prof_t = define_profile(x_profile, y_profile, resolution[0])
        dist_prof_t = np.round(dist_prof_t)
        profile = np.interp(dist_prof_t,dist_prof,profile)
        profile = np.round(profile,decimals=-1)
    elif do_sediments and i==layers:
        x_int_t, y_int_t, dist_prof_t = define_profile(x_profile, y_profile, resolution[0])
        interp_func = RegularGridInterpolator((yvals,xvals),np.asarray(data))
        profile = interp_func((y_int,x_int),method="linear").tolist()
        profile = np.interp(dist_prof_t,dist_prof,profile)
    else:
        x_int_t, y_int_t, dist_prof_t = define_profile(x_profile, y_profile, resolution[0])
        dist_prof = np.round(dist_prof)
        profile = np.interp(dist_prof_t,dist_prof,profile)
        profile = np.round(profile,decimals=-1)
    
    return dist_prof,profile




def in_area_s(acs,x,y,g):
    '''
    Limit the size of the originial dataset to make interpolation easier

    Parameters
    ----------
    acs : research area coordinates
    x : xi of grid
    y : yi of grid
    g : data grid

    Returns
    -------
    x_s : smaller xi of grid
    y_s : smaller yi of grid
    grid_s : smaller data grid
    '''
    
    x = np.unique(x)
    y = np.unique(y)
    
    x_s = x[((acs[0]-50000)<=x) & (x<=(acs[1]+50000))]
    y_s = y[((acs[2]+50000)>=y) & (y>=(acs[3]-50000))]
    
    grid_s = g[np.ix_(((acs[2]+50000>=y)) & (y>=(acs[3]-50000)),((acs[0]-50000)<=x) & (x<=(acs[1]+50000)))]
    grid_s = np.transpose(grid_s)
    
    return x_s,y_s,grid_s

def fix_surface_height(surface,x_data,y_data,data):
    '''
    Set the height of the 3d surface to the correct interface depth

    Parameters
    ----------
    surface : the flat 3d mesh
    x_data : unique x coordinates of mesh
    y_data : unique y coordinates of mesh
    data : input data in area

    Returns
    -------
    surface : height corrected 3d surface
    '''
    
    for node in surface.nodes():
        x,y,z = node.pos()
        idx,idy = np.argmin(np.abs(x_data-x)),np.argmin(np.abs(y_data-y))
        node.setPos(pg.RVector3(x,y,data[idy,idx]))
    
    return surface

def compare_to_plane(x,y,data,grid_xy,tree):
    '''
    check the height of the plane of the three neighbouring points for the cell center at x,y(,z)
    this is rather simple algebra but took me like 20 attempts to get right

    Parameters
    ----------
    x : x coord of node
    y : y coord of node
    x_int : x_int of data grid
    y_int : y_int of data grid
    data : grid of topo or bottom sed

    Returns
    -------
    F : height of plane at x,y
    '''
    
    km = 1000
    
    distances,indices = tree.query([x,y],k=3)
    nearest_points = grid_xy[indices]
    height_data = np.transpose(data).flatten('F')[indices]/km

    A = [nearest_points[0,0],nearest_points[0,1],height_data[0]]
    B = [nearest_points[1,0],nearest_points[1,1],height_data[1]]
    C = [nearest_points[2,0],nearest_points[2,1],height_data[2]]
    
    AB = [-A[0]+B[0],-A[1]+B[1],-A[2]+B[2]]
    AC = [-A[0]+C[0],-A[1]+C[1],-A[2]+C[2]]
    
    cross = np.cross(AB,AC)
    
    norm = cross/(np.sqrt(cross[0]**2+cross[1]**2+cross[2]**2))
    
    D = -(A[0]*norm[0]+A[1]*norm[1]+A[2]*norm[2])
    
    F = (-x*norm[0] - y*norm[1] - D)/norm[2]

    return F


def define_profile(x_profile,y_profile,dist_int):
    '''
    generate the coordinate profile to a given resolution

    Parameters
    ----------
    x_profile : x0 and x1
    y_profile : y0 and y1 - both in south polar azimuthal equidistant projection
    dist_int : input data resolution

    Returns
    -------
    x_int : x coordinates of new profile
    y_int : y coordinates of new profile
    dist_prof : profile coordinates in resolution (as close as possible, usual discrepancies around .5%)
    '''
    
    # get direction that profile is facing
    if x_profile[0] < x_profile[1]:
        facdx = 1
    elif x_profile[0] > x_profile[1]:
        facdx = -1
    if y_profile[0] < y_profile[1]:
        facdy = 1
    elif y_profile[0] > y_profile[1]:
        facdy = -1
    # case fpr horizontal/vertical profile
    if x_profile[0] == x_profile[1]:
        dy = dist_int*facdy
        y_int=np.arange(y_profile[0],y_profile[1]+dy,dy)
        # y_int=np.arange(y_profile[0]-y_profile[0],y_profile[1]-y_profile[0]+dy,dy)
        x_int=np.full(len(y_int),x_profile[0])
        dist_prof=np.copy(y_int-y_int[0])
        return x_int,y_int,dist_prof
    if y_profile[0] == y_profile[1]:
        dx = dist_int*facdx
        x_int=np.arange(x_profile[0],x_profile[1]+dx,dx)
        # x_int=np.arange(x_profile[0]-x_profile[0],x_profile[1]-x_profile[0]+dx,dx)
        y_int=np.full(len(x_int),y_profile[0])
        dist_prof=np.copy(x_int-x_int[0])
        return x_int,y_int,dist_prof
    # normal case
    m=(y_profile[1]-y_profile[0])/(x_profile[1]-x_profile[0])
    dx=np.sqrt(dist_int**2/(1+m**2))*facdx
    dy=np.sqrt(dist_int**2-dx**2)*facdy
    prof_len = np.sqrt((x_profile[1]-x_profile[0])**2 + (y_profile[1]-y_profile[0])**2)
    x_int = np.linspace(x_profile[0],x_profile[1],np.round(prof_len/dist_int).astype(int))
    y_int = np.linspace(y_profile[0],y_profile[1],np.round(prof_len/dist_int).astype(int))
    
    # x_int=np.arange(x_profile[0],x_profile[1]+dx,dx)
    # y_int=np.arange(y_profile[0],y_profile[1]+dy,dy)
    # assure arrays are same lengths, for some dist_int this is not true
    if len(x_int) < len(y_int):
        x_int.extend(x_int[-1]+dx)
    if len(x_int) > len(y_int):
        y_int.extend(y_int[-1]+dy)
    # create distance profile in here because otherwise its a headache
    dist_prof=np.linspace(0,prof_len,np.round(prof_len/dist_int).astype(int))
    # dist_prof=np.arange(0,np.sqrt((x_int[-1]-x_int[0])**2+(y_int[-1]-y_int[0])**2)+dist_int,dist_int)
    if len(dist_prof) != len(x_int):
        dist_prof = dist_prof[:-1]
    return x_int,y_int,dist_prof

def create_sed_interface(x_int_t,y_int_t,grid_topo,grid_sed):
    '''
    create the costum interface for the sediment layer, which needs to lie under the topography

    Parameters
    ----------
    x_int_t : gridded x coordinates of topography and sediment thickness
    y_int_t : gridded y coordinates of topography and sediment thickness
    grid_topo : gridded topography
    grid_sed : gridded sediment_thickness

    Returns
    -------
    surface_sed : pyGIMLi 3D surface
    missing : the list of coordinates of edge nodes where sed = 0

    '''
    km = 1000
    triangles_topo = []
    triangles_sed = []
    grid_topo = np.transpose(grid_topo)
    grid_sed = np.transpose(grid_sed)
    for j in range(np.shape(grid_sed)[1]-1):
        for i in range(np.shape(grid_sed)[0]-1):
            zero_count = sum(1 for x in [grid_sed[i,j],grid_sed[i+1,j],grid_sed[i,j+1],grid_sed[i+1,j+1]] if x == 0)
            if zero_count == 4:
                # sed
                # none
                
                # topo
                triangle = pg.Mesh(3,isGeometry=True)
                triangle.createNode(x_int_t[i]/km,y_int_t[j]/km,grid_topo[i,j]/km,marker=4)
                triangle.createNode(x_int_t[i+1]/km,y_int_t[j]/km,grid_topo[i+1,j]/km,marker=4)
                triangle.createNode(x_int_t[i]/km,y_int_t[j+1]/km,grid_topo[i,j+1]/km,marker=4)
                surf = [0,1,2]
                triangle.createPolygonFace(triangle.nodes(surf),marker=4)
                triangles_topo.append(triangle)
                
                triangle = pg.Mesh(3,isGeometry=True)
                triangle.createNode(x_int_t[i+1]/km,y_int_t[j]/km,grid_topo[i+1,j]/km,marker=4)
                triangle.createNode(x_int_t[i]/km,y_int_t[j+1]/km,grid_topo[i,j+1]/km,marker=4)
                triangle.createNode(x_int_t[i+1]/km,y_int_t[j+1]/km,grid_topo[i+1,j+1]/km,marker=4)
                surf = [0,1,2]
                triangle.createPolygonFace(triangle.nodes(surf),marker=4)
                triangles_topo.append(triangle)
                continue
            elif zero_count == 3:
                if grid_sed[i,j] != 0:
                    # sed
                    triangle = pg.Mesh(3,isGeometry=True)
                    triangle.createNode([x_int_t[i]/km,y_int_t[j]/km,(grid_topo[i,j]-grid_sed[i,j])/km],marker=5)
                    triangle.createNode([x_int_t[i+1]/km,y_int_t[j]/km,(grid_topo[i+1,j]-grid_sed[i+1,j])/km],marker=9)
                    triangle.createNode([x_int_t[i]/km,y_int_t[j+1]/km,(grid_topo[i,j+1]-grid_sed[i,j+1])/km],marker=9)
                    surf = [0,1,2]
                    triangle.createPolygonFace(triangle.nodes(surf),marker=5)
                    triangles_sed.append(triangle)
                    
                    # topo
                    triangle = pg.Mesh(3,isGeometry=True)
                    triangle.createNode(x_int_t[i]/km,y_int_t[j]/km,grid_topo[i,j]/km,marker=4)
                    triangle.createNode(x_int_t[i+1]/km,y_int_t[j]/km,grid_topo[i+1,j]/km,marker=9)
                    triangle.createNode(x_int_t[i]/km,y_int_t[j+1]/km,grid_topo[i,j+1]/km,marker=9)
                    surf = [0,1,2]
                    triangle.createPolygonFace(triangle.nodes(surf),marker=4)
                    triangles_topo.append(triangle)
                    
                    triangle = pg.Mesh(3,isGeometry=True)
                    triangle.createNode(x_int_t[i+1]/km,y_int_t[j]/km,grid_topo[i+1,j]/km,marker=9)
                    triangle.createNode(x_int_t[i]/km,y_int_t[j+1]/km,grid_topo[i,j+1]/km,marker=9)
                    triangle.createNode(x_int_t[i+1]/km,y_int_t[j+1]/km,grid_topo[i+1,j+1]/km,marker=9)
                    surf = [0,1,2]
                    triangle.createPolygonFace(triangle.nodes(surf),marker=4)
                    triangles_topo.append(triangle)
                    continue
                elif grid_sed[i+1,j] != 0:
                    # sed
                    triangle = pg.Mesh(3,isGeometry=True)
                    triangle.createNode([x_int_t[i+1]/km,y_int_t[j]/km,(grid_topo[i+1,j]-grid_sed[i+1,j])/km],marker=5)
                    triangle.createNode([x_int_t[i]/km,y_int_t[j]/km,(grid_topo[i,j]-grid_sed[i,j])/km],marker=9)
                    triangle.createNode([x_int_t[i+1]/km,y_int_t[j+1]/km,(grid_topo[i+1,j+1]-grid_sed[i+1,j+1])/km],marker=9)
                    surf = [0,1,2]
                    triangle.createPolygonFace(triangle.nodes(surf),marker=5)
                    triangles_sed.append(triangle)
                    
                    # topo
                    triangle = pg.Mesh(3,isGeometry=True)
                    triangle.createNode(x_int_t[i]/km,y_int_t[j]/km,grid_topo[i,j]/km,marker=9)
                    triangle.createNode(x_int_t[i+1]/km,y_int_t[j]/km,grid_topo[i+1,j]/km,marker=4)
                    triangle.createNode(x_int_t[i+1]/km,y_int_t[j+1]/km,grid_topo[i+1,j+1]/km,marker=9)
                    surf = [0,1,2]
                    triangle.createPolygonFace(triangle.nodes(surf),marker=4)
                    triangles_topo.append(triangle)
                    
                    triangle = pg.Mesh(3,isGeometry=True)
                    triangle.createNode(x_int_t[i]/km,y_int_t[j]/km,grid_topo[i,j]/km,marker=9)
                    triangle.createNode(x_int_t[i]/km,y_int_t[j+1]/km,grid_topo[i,j+1]/km,marker=9)
                    triangle.createNode(x_int_t[i+1]/km,y_int_t[j+1]/km,grid_topo[i+1,j+1]/km,marker=9)
                    surf = [0,1,2]
                    triangle.createPolygonFace(triangle.nodes(surf),marker=4)
                    triangles_topo.append(triangle)
                    continue
                elif grid_sed[i,j+1] != 0:
                    # sed
                    triangle = pg.Mesh(3,isGeometry=True)
                    triangle.createNode([x_int_t[i]/km,y_int_t[j+1]/km,(grid_topo[i,j+1]-grid_sed[i,j+1])/km],marker=5)
                    triangle.createNode([x_int_t[i]/km,y_int_t[j]/km,(grid_topo[i,j]-grid_sed[i,j])/km],marker=9)
                    triangle.createNode([x_int_t[i+1]/km,y_int_t[j+1]/km,(grid_topo[i+1,j+1]-grid_sed[i+1,j+1])/km],marker=9)
                    surf = [0,1,2]
                    triangle.createPolygonFace(triangle.nodes(surf),marker=5)
                    triangles_sed.append(triangle)
                    
                    # topo
                    triangle = pg.Mesh(3,isGeometry=True)
                    triangle.createNode(x_int_t[i]/km,y_int_t[j+1]/km,grid_topo[i,j+1]/km,marker=4)
                    triangle.createNode(x_int_t[i]/km,y_int_t[j]/km,grid_topo[i,j]/km,marker=9)
                    triangle.createNode(x_int_t[i+1]/km,y_int_t[j+1]/km,grid_topo[i+1,j+1]/km,marker=9)
                    surf = [0,1,2]
                    triangle.createPolygonFace(triangle.nodes(surf),marker=4)
                    triangles_topo.append(triangle)
                    
                    triangle = pg.Mesh(3,isGeometry=True)
                    triangle.createNode(x_int_t[i+1]/km,y_int_t[j]/km,grid_topo[i+1,j]/km,marker=9)
                    triangle.createNode(x_int_t[i]/km,y_int_t[j]/km,grid_topo[i,j]/km,marker=9)
                    triangle.createNode(x_int_t[i+1]/km,y_int_t[j+1]/km,grid_topo[i+1,j+1]/km,marker=9)
                    surf = [0,1,2]
                    triangle.createPolygonFace(triangle.nodes(surf),marker=4)
                    triangles_topo.append(triangle)
                    continue
                elif grid_sed[i+1,j+1] != 0:
                    # sed
                    triangle = pg.Mesh(3,isGeometry=True)
                    triangle.createNode([x_int_t[i+1]/km,y_int_t[j+1]/km,(grid_topo[i+1,j+1]-grid_sed[i+1,j+1])/km],marker=5)
                    triangle.createNode([x_int_t[i+1]/km,y_int_t[j]/km,(grid_topo[i+1,j]-grid_sed[i+1,j])/km],marker=9)
                    triangle.createNode([x_int_t[i]/km,y_int_t[j+1]/km,(grid_topo[i,j+1]-grid_sed[i,j+1])/km],marker=9)
                    surf = [0,1,2]
                    triangle.createPolygonFace(triangle.nodes(surf),marker=5)
                    triangles_sed.append(triangle)
                    
                    # topo
                    triangle = pg.Mesh(3,isGeometry=True)
                    triangle.createNode(x_int_t[i]/km,y_int_t[j]/km,grid_topo[i,j]/km,marker=9)
                    triangle.createNode(x_int_t[i+1]/km,y_int_t[j]/km,grid_topo[i+1,j]/km,marker=9)
                    triangle.createNode(x_int_t[i]/km,y_int_t[j+1]/km,grid_topo[i,j+1]/km,marker=9)
                    surf = [0,1,2]
                    triangle.createPolygonFace(triangle.nodes(surf),marker=4)
                    triangles_topo.append(triangle)
                    
                    triangle = pg.Mesh(3,isGeometry=True)
                    triangle.createNode(x_int_t[i+1]/km,y_int_t[j]/km,grid_topo[i+1,j]/km,marker=9)
                    triangle.createNode(x_int_t[i]/km,y_int_t[j+1]/km,grid_topo[i,j+1]/km,marker=9)
                    triangle.createNode(x_int_t[i+1]/km,y_int_t[j+1]/km,grid_topo[i+1,j+1]/km,marker=4)
                    surf = [0,1,2]
                    triangle.createPolygonFace(triangle.nodes(surf),marker=4)
                    triangles_topo.append(triangle)
                    continue
            else:
                # sed
                triangle = pg.Mesh(3,isGeometry=True)
                triangle.createNode(x_int_t[i]/km,y_int_t[j]/km,(grid_topo[i,j]-grid_sed[i,j])/km,marker=5)
                triangle.createNode(x_int_t[i+1]/km,y_int_t[j]/km,(grid_topo[i+1,j]-grid_sed[i+1,j])/km,marker=5)
                triangle.createNode(x_int_t[i]/km,y_int_t[j+1]/km,(grid_topo[i,j+1]-grid_sed[i,j+1])/km,marker=5)
                surf = [0,1,2]
                triangle.createPolygonFace(triangle.nodes(surf),marker=5)
                triangles_sed.append(triangle)
                
                triangle = pg.Mesh(3,isGeometry=True)
                triangle.createNode(x_int_t[i+1]/km,y_int_t[j]/km,(grid_topo[i+1,j]-grid_sed[i+1,j])/km,marker=5)
                triangle.createNode(x_int_t[i]/km,y_int_t[j+1]/km,(grid_topo[i,j+1]-grid_sed[i,j+1])/km,marker=5)
                triangle.createNode(x_int_t[i+1]/km,y_int_t[j+1]/km,(grid_topo[i+1,j+1]-grid_sed[i+1,j+1])/km,marker=5)
                surf = [0,1,2]
                triangle.createPolygonFace(triangle.nodes(surf),marker=5)
                triangles_sed.append(triangle)
                
                # topo
                triangle = pg.Mesh(3,isGeometry=True)
                triangle.createNode(x_int_t[i]/km,y_int_t[j]/km,grid_topo[i,j]/km,marker=4)
                triangle.createNode(x_int_t[i+1]/km,y_int_t[j]/km,grid_topo[i+1,j]/km,marker=4)
                triangle.createNode(x_int_t[i]/km,y_int_t[j+1]/km,grid_topo[i,j+1]/km,marker=4)
                surf = [0,1,2]
                triangle.createPolygonFace(triangle.nodes(surf),marker=4)
                triangles_topo.append(triangle)
                
                triangle = pg.Mesh(3,isGeometry=True)
                triangle.createNode(x_int_t[i+1]/km,y_int_t[j]/km,grid_topo[i+1,j]/km,marker=4)
                triangle.createNode(x_int_t[i]/km,y_int_t[j+1]/km,grid_topo[i,j+1]/km,marker=4)
                triangle.createNode(x_int_t[i+1]/km,y_int_t[j+1]/km,grid_topo[i+1,j+1]/km,marker=4)
                surf = [0,1,2]
                triangle.createPolygonFace(triangle.nodes(surf),marker=4)
                triangles_topo.append(triangle)
                continue
    
    # merge all
    surface_topo = mt.mergePLC(triangles_topo)
    surface_sed = mt.mergePLC(triangles_sed)
    return surface_topo, surface_sed

#%% top be discarded:
    
def load_data():
    '''
    load the input data with this external function to make the actual scripts slightly smaller

    Returns
    -------
    list gridded x and y, and the data grid for topo, sed thickness, moho depth, and lab depth
    list of data grid resolutions
    '''
    # bedmap
    
    data_topo = pd.read_csv('data/bedmap2_bed.txt', skiprows=6, sep=' ',header=None)
    data_topo = data_topo.iloc[:,:-1]
    data_topo = data_topo.iloc[::-1].reset_index(drop=True)
    grid_topo = np.asarray(data_topo)

    # x_sp_topo = np.arange(-3333500,3333500,1000)
    # y_sp_topo = np.arange(-3333500,3333500,1000)

    # xi_sp_topo,yi_sp_topo = np.meshgrid(x_sp_topo,y_sp_topo)

    # transformer = Transformer.from_crs("EPSG:3031","ESRI:102019")
    # xi_topo = np.zeros(np.shape(xi_sp_topo))
    # yi_topo = np.zeros(np.shape(yi_sp_topo))
    # for i in range(np.shape(xi_topo)[0]):
    #     for j in range(np.shape(xi_topo)[1]):
    #         xi_topo[i,j],yi_topo[i,j] = transformer.transform(xi_sp_topo[i,j],yi_sp_topo[i,j])
            
    x_ae_topo = np.linspace(-3348244,3348244,6667)
    y_ae_topo = np.linspace(-3348244,3348244,6667)
    
    xi_topo,yi_topo = np.meshgrid(x_ae_topo,y_ae_topo)
            
    # sediment

    data_sed = pd.read_csv('data/sedimentary_layers_in_equidistant_projection.dat', delimiter=' ')
    y_sed = data_sed.phy*111e3
    x_sed = data_sed.theta*111e3
    
    x_ae_sed = np.arange(np.min(x_sed),np.max(x_sed)+9250,9250)
    y_ae_sed = np.arange(np.min(y_sed),np.max(y_sed)+9250,9250)

    xi_sed, yi_sed = np.meshgrid(x_ae_sed,y_ae_sed)

    grid_sed = griddata((x_sed,y_sed),data_sed.total_thickness,(xi_sed,yi_sed),method='linear')*1000
    
    # moho
    
    data_H_Moho = pd.read_csv("data/Haeger_Moho.csv", sep=',', comment="#")
    data_H_Moho['Lat'] = data_H_Moho['Lat'].apply(lambda x: -x if x > 0 else x)
    
    myProj = Proj("+proj=aeqd +lat_0=-90")
    x_moho, y_moho  = myProj(data_H_Moho.Lon, data_H_Moho.Lat)
    
    x_ae_moho = np.arange(np.min(x_moho),np.max(x_moho)+10000,10000)
    y_ae_moho = np.arange(np.min(y_moho),np.max(y_moho)+10000,10000)
     
    xi_moho, yi_moho = np.meshgrid(x_ae_moho,y_ae_moho)
    
    grid_moho = griddata((x_moho,y_moho),data_H_Moho.Moho,(xi_moho,yi_moho))*1000
    
    # lab
    
    data_H_LAB = pd.read_csv("data/Haeger_LAB.csv", sep=',', comment="#")
    data_H_LAB['Lat'] = data_H_LAB['Lat'].apply(lambda x: -x if x > 0 else x)
    
    x_lab, y_lab = myProj(data_H_LAB.Lon, data_H_LAB.Lat)
    
    x_ae_lab = np.arange(np.min(x_lab),np.max(x_lab)+10000,10000)
    y_ae_lab = np.arange(np.min(y_lab),np.max(y_lab)+10000,10000)
    
    xi_lab, yi_lab = np.meshgrid(x_ae_lab,y_ae_lab)
    
    grid_lab = griddata((x_lab,y_lab),data_H_LAB.LAB,(xi_lab,yi_lab))*1000
    
    return [[xi_topo,yi_topo,grid_topo],[xi_sed,yi_sed,grid_sed],[xi_moho,yi_moho,grid_moho],[xi_lab,yi_lab,grid_lab]],[1000,9250,10000,10000]









