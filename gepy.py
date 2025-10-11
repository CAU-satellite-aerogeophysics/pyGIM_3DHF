"""
Georg HÃ¼ttner - October 2024

One of probably a bunch of scripts to import dumbass functions
"""

import pygimli as pg
import pygimli.meshtools as mt
import numpy as np
from scipy.interpolate import griddata, RegularGridInterpolator

km = 1000

# def load_and_interpolate_2d()
# def load_and_interpolate_3d()


# class build_world ? 
#   would maybe allow for things like setting up differen thermal parameters for the same mesh
# def build_world()


# def assign t_params

# def calc_ghf()

def assign_markers_3d(data_list,
                      mesh,
                      layers = None,
                      do_sediments = False,
                      do_hp = False,
                      hp_0 = None,
                      tc = None
                      ):
    '''
    take the previously created 3d pygimli mesh and assign thermal parameters to cells and nodes

    Parameters
    ----------
    data_list :    data input list of [xi, yi, grid]
    mesh :         pygimli mesh from build_world
    layers :       number of layers
    do_sediments : is a sediment layer present?
    do_hp :        is a heat production distribution provided?
    hp_0 :         which value to take in case no hp is provided
    tc :           list of thermal conductivities 

    Returns
    -------
    mesh : mesh with assigned markers and parameters
    '''
    
    return mesh
    
def assign_markers_2d(data_list,
                      mesh,
                      layers = None,
                      do_sediments = False,
                      do_hp = False,
                      hp_0 = None,
                      tc = None
                      ):
    '''
    take the previously created 2d pygimli mesh and assign thermal parameters to cells and nodes

    Parameters
    ----------
    data_list :    data input list of [xi, yi, grid]
    mesh :         pygimli mesh from build_world
    layers :       number of layers
    do_sediments : is a sediment layer present?
    do_hp :        is a heat production distribution provided?
    hp_0 :         which value to take in case no hp is provided
    tc :           list of thermal conductivities 

    Returns
    -------
    mesh : mesh with assigned markers and parameters
    '''
    
    return mesh

def build_world_2d(data_list          ,
                  profile_coord       ,
                  area = None         ,
                  layers = None       ,
                  do_sediments = False          
                  ):
    '''
    Build the 2d pygimli mesh from the provided dataset.
    
    Parameters
    ----------
    data_list :     data input list of [xi, yi, grid]
    profile_coord : profile coordinates
    area :          maximum size of mesh triangles
    layers :        number of layers 
    do_sediments :  is a sediment layer present?

    Returns
    -------
    mesh : the pygimli mesh of the world with all layers
    line_list : pygimli geometries of all layers for inspection
    '''
    profile_length = np.round(np.sqrt((profile_coord[0][1]-profile_coord[0][0])**2 + (profile_coord[1][1]-profile_coord[1][0])**2))
    
    world_start = [0,4*km] 
    world_end = [profile_length,-220*km]

    world = mt.createWorld(start=world_start,end=world_end, worldMarker=False)
    
    line_list = []
    
    topo_line = mt.createPolygon([[a,b] for a,b in zip(data_list[0][0],data_list[0][1])],marker = 5,boundaryMarker = 5,isClosed=False)
    line_list.append(topo_line)
    if layers == 2:
        lab_line = mt.createPolygon([[a,b] for a,b in zip(data_list[1][0],-data_list[1][1])],marker = 8,boundaryMarker = 8,isClosed=False)
        line_list.append(lab_line)
    elif do_sediments and layers == 3:
        sed_line = mt.createPolygon([[a,b] for a,b in zip(data_list[1][0],data_list[1][1])],marker = 6,boundaryMarker = 6,isClosed=False) 
        lab_line = mt.createPolygon([[a,b] for a,b in zip(data_list[2][0],-data_list[2][1])],marker = 8,boundaryMarker = 8,isClosed=False)
        line_list.append(sed_line)
        line_list.append(lab_line)
    elif not do_sediments and layers == 3:
        moho_line = mt.createPolygon([[a,b] for a,b in zip(data_list[1][0],data_list[1][1])],marker = 7,boundaryMarker = 6,isClosed=False) 
        lab_line = mt.createPolygon([[a,b] for a,b in zip(data_list[2][0],-data_list[2][1])],marker = 8,boundaryMarker = 8,isClosed=False)
        line_list.append(moho_line)
        line_list.append(lab_line)
    elif do_sediments and layers == 3:
        sed_line = mt.createPolygon([[a,b] for a,b in zip(data_list[1][0],data_list[1][1])],marker = 6,boundaryMarker = 6,isClosed=False)
        moho_line = mt.createPolygon([[a,b] for a,b in zip(data_list[2][0],-data_list[2][1])],marker = 7,boundaryMarker = 7,isClosed=False)
        lab_line = mt.createPolygon([[a,b] for a,b in zip(data_list[3][0],-data_list[3][1])],marker = 8,boundaryMarker = 8,isClosed=False)
        line_list.append(sed_line)
        line_list.append(moho_line)
        line_list.append(lab_line)
    
    for i in line_list:
        world = world + i
    
    if area != None:
        mesh = mt.createMesh(world,quality=34,area=area)
    else:
        mesh = mt.createMesh(world,quality=34,area=profile_length*4)
    
    return mesh, line_list

def build_world_3d(data_list          ,
                  area_coords         ,
                  area = None         ,
                  layers = None       ,
                  do_sediments = False,
                  border = None       
                  ):
    '''
    Build the 3d pygimli mesh from the provided dataset.
    
    Parameters
    ----------
    data_list :     data input list of [xi, yi, grid]
    area_coords :   research area coordinates
    area :          maximum size of mesh triangles
    layers :        number of layers 
    do_sediments :  is a sediment layer present?
    border : buffer between input data and world edge, can cause issues otherwise

    Returns
    -------
    mesh : the pygimli mesh of the world with all layers
    layer_list : pygimli geometries of all layers for inspection
    '''
    
    start = [(area_coords[0]-border)/km,(area_coords[2]+border)/km,4] 
    end = [(area_coords[1]+border)/km,(area_coords[3]-border)/km,-220]
    
    world = mt.createWorld(start=start,end=end, worldMarker=False)
    
    # create layer polygons and shift them to the correct height
    
    layer_list = []
    
    if do_sediments:
        surface_topo,surface_sed = create_sed_interface(data_list[0][0], data_list[0][1], data_list[0][2], data_list[1][2])
        layer_list.append(surface_topo)
        layer_list.append(surface_sed)
    elif not do_sediments:
        mesh_topo = mt.createMesh2D(data_list[0][0]/km,data_list[0][1]/km)
        surface_topo = mt.createSurface(mesh_topo)
        surface_topo = fix_surface_height(surface_topo, data_list[0][0]/km,data_list[0][1]/km, data_list[0][2]/km)
        layer_list.append(surface_topo)
    else:
        print("wrong do_sediments input")
    
    if layers == 3 and not do_sediments:
        mesh_moho = mt.createMesh2D(data_list[1][0]/km,data_list[1][1]/km)
        surface_moho = mt.createSurface(mesh_moho)
        surface_moho = fix_surface_height(surface_moho, data_list[1][0]/km,data_list[1][1]/km, -data_list[1][2]/km)
        layer_list.append(surface_moho)
    
        mesh_lab = mt.createMesh2D(data_list[2][0]/km,data_list[2][1]/km)
        surface_lab = mt.createSurface(mesh_lab)
        surface_lab = fix_surface_height(surface_lab, data_list[2][0]/km,data_list[2][1]/km, -data_list[2][2]/km)
        layer_list.append(surface_lab)
    elif layers == 3 and do_sediments:
        mesh_lab = mt.createMesh2D(data_list[2][0]/km,data_list[2][1]/km)
        surface_lab = mt.createSurface(mesh_lab)
        surface_lab = fix_surface_height(surface_lab, data_list[2][0]/km,data_list[2][1]/km, -data_list[2][2]/km)
        layer_list.append(surface_lab)
    else:
        print("wrong do_sediments input")
    
    if layers == 4:
        mesh_moho = mt.createMesh2D(data_list[2][0]/km,data_list[2][1]/km)
        surface_moho = mt.createSurface(mesh_moho)
        surface_moho = fix_surface_height(surface_moho, data_list[2][0]/km,data_list[2][1]/km, -data_list[2][2]/km)
        layer_list.append(surface_moho)
    
        mesh_lab = mt.createMesh2D(data_list[3][0]/km,data_list[3][1]/km)
        surface_lab = mt.createSurface(mesh_lab)
        surface_lab = fix_surface_height(surface_lab, data_list[3][0]/km,data_list[3][1]/km, -data_list[3][2]/km)
        layer_list.append(surface_lab)
    
    geometry = world
    for i in layer_list:
        geometry = geometry + i
    
    mesh = mt.createMesh(geometry,quality=34,area=area)

    return mesh,layer_list
    
def cut_and_interpolate_3d(layer               ,
                           resolution          ,
                           i                   ,
                           layers              ,
                           area_coords         ,
                           do_sediments = False,
                           do_hp = False       ,
                           hp_0 = None
                           ):
    '''
    Take input format-correct input data, cut it to the research area and interpolate all layers through research area

    Parameters
    ----------
    layer :         data input with xi, yi, grid
    resolution :    list of data resolution in m
    i :             loop counter for layer selection
    layers :        number of layers 
    area_coords :   research area coordinates
    profile_coord : profile coordinates
    do_sediments :  is a sediment layer present?
    do_hp :         is a heat production distribution provided?
    hp_0 :          which value to take in case no hp is provided

    Returns
    -------
    x_int : gridded x coordinate of data grid
    y_int : gridded y coordinate of data grid
    grid : data grid
    '''
    
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
    '''
    Take input format-correct input data, cut it to the research area and interpolate all layers onto their respective profiles

    Parameters
    ----------
    layer :         data input with xi, yi, grid
    resolution :    list of data resolution in m
    i :             loop counter for layer selection
    layers :        number of layers 
    area_coords :   research area coordinates
    profile_coord : profile coordinates
    do_sediments :  is a sediment layer present?
    do_hp :         is a heat production distribution provided
    hp_0 :          which value to take in case no hp is provided

    Returns
    -------
    dist_prof : profile coordinates of specific layer
    profile : data along profile
    '''
    
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
    
    #if not do_hp:
     # implement the no hp case   
    
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
