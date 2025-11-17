"""
Georg Hüttner - November 2025

pyGIM_3DGHF functions that rely on pyGIMLi and its meshtools
"""

import pygimli as pg
import pygimli.meshtools as mt
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import KDTree

km = 1000

#%% calc_ghf
def calc_ghf(T,
             mesh,
             data_list,
             dimension,
             tc,
             do_sediments = False,
             ):

    if dimension == '3D':
        topcoord = []
        x_int = np.unique(data_list[0][0])
        y_int = np.unique(data_list[0][1])
        for i in range(len(x_int)):
            for j in range(len(y_int)):
                topcoord.append([x_int[i]/km,y_int[j]/km,(data_list[0][2][j,i]-5)/km])
                
        gradientTop = pg.solver.grad(mesh, T, topcoord)
        shf = np.zeros(len(gradientTop))
        for i in range(len(gradientTop)):
            point = pg.RVector3(topcoord[i][0],topcoord[i][1],topcoord[i][2])
            cell = mesh.findCell(point)
            geology = cell.marker()
            if geology == 5:   
                shf[i] = np.linalg.norm(gradientTop[i,:])*tc[1]
            elif geology == 4:
                shf[i] = np.linalg.norm(gradientTop[i,:])*tc[0]
        shf = np.transpose(np.reshape(shf,(len(x_int),len(y_int))))
        
    elif dimension == '2D':
        topcoord = []
        for i in range(len(data_list[0][0])):
            topcoord.append([data_list[0][0][i],data_list[0][1][i]-1,0])
            
        gradientTop = pg.solver.grad(mesh, T, topcoord)
        shf = np.zeros(len(gradientTop))
        for i in range(len(gradientTop)):
            point = pg.RVector3(topcoord[i][0],topcoord[i][1],topcoord[i][2])
            cell = mesh.findCell(point)
            geology = cell.marker()
            if geology == 6:   
                shf[i] = np.linalg.norm(gradientTop[i,:])*tc[1]
            elif geology == 5:
                shf[i] = np.linalg.norm(gradientTop[i,:])*tc[0]
        shf = shf*1000 
    else:
        print('Please give dimension')
    
    return shf
    
#%% calc_temp
def calc_temp(mesh,
              dimension,
              layers,
              force = None,
              do_sediments = False,
              do_hp = False,
              tc = None,
              hp_0 = None,
              Tbound = None
              ):
    '''
    Calculate the temperature field on the supplied mesh

    Parameters
    ----------
    mesh :         pygimli mesh from build_world, with assigned markers
    dimension :    '3D' or '2D' ?
    force :        nodal heat production values from assign_markers_
    layers :       number of layers
    do_sediments : is a sediment layer present?
    do_hp :        REF/CONST/NONE - how to deal with heat production
    tc :           list of thermal conductivities
    hp_0 :         which CONST value to take for heat production
    Tbound :       optional costum boundary conditions for [T_LAB , T_surface]

    Returns
    -------
    T : temperature field on mesh

    '''
    
    if Tbound == None:
        T_LAB = 1315
        T_surface = 0
    else:
        T_LAB = Tbound[0]
        T_surface = Tbound[1]
    
    if dimension == '3D':
        if layers == 2:
            if do_hp == 'REF' or do_hp == 'CONST':
                T = pg.solver.solveFiniteElements(mesh,
                                                  a={1: 1.0*km, 4: tc[0]*km, 7: 4.0*km},
                                                  f=force*km*km*km,
                                                  bc={'Dirichlet': {10: T_LAB, 7: T_surface}},verbose=True)
            elif do_hp == 'NONE':
                T = pg.solver.solveFiniteElements(mesh,
                                                  a={1: 1.0*km, 4: tc[0]*km, 7: 4.0*km},
                                                  bc={'Dirichlet': {10: T_LAB, 7: T_surface}},verbose=True)
        elif do_sediments and layers == 3:
            if do_hp == 'REF' or do_hp == 'CONST':
                T = pg.solver.solveFiniteElements(mesh,
                                                  a={1: 1.0*km, 4: tc[0]*km, 5: tc[1]*km, 7: 4.0*km},
                                                  f=force*km*km*km,
                                                  bc={'Dirichlet': {10: T_LAB, 7: T_surface}},verbose=True)
            elif do_hp == 'NONE':
                T = pg.solver.solveFiniteElements(mesh,
                                                  a={1: 1.0*km, 4: tc[0]*km, 5: tc[1]*km, 7: 4.0*km},
                                                  bc={'Dirichlet': {10: T_LAB, 7: T_surface}},verbose=True)
        elif not do_sediments and layers == 3:
            if do_hp == 'REF' or do_hp == 'CONST':
                T = pg.solver.solveFiniteElements(mesh,
                                                  a={1: 1.0*km, 4: tc[0]*km, 6: tc[1]*km, 7: 4.0*km},
                                                  f=force*km*km*km,
                                                  bc={'Dirichlet': {10: T_LAB, 7: T_surface}},verbose=True)
            elif do_hp == 'NONE':
                T = pg.solver.solveFiniteElements(mesh,
                                                  a={1: 1.0*km, 4: tc[0]*km, 6: tc[1]*km, 7: 4.0*km},
                                                  bc={'Dirichlet': {10: T_LAB, 7: T_surface}},verbose=True)
        elif do_sediments and layers == 4:
            if do_hp == 'REF' or do_hp == 'CONST':
                T = pg.solver.solveFiniteElements(mesh,
                                                  a={1: 1.0*km, 4: tc[0]*km, 5: tc[1]*km, 6: tc[2]*km, 7: 4.0*km},
                                                  f=force*km*km*km,
                                                  bc={'Dirichlet': {10: T_LAB, 7: T_surface}},verbose=True)
            elif do_hp == 'NONE':
                T = pg.solver.solveFiniteElements(mesh,
                                                  a={1: 1.0*km, 4: tc[0]*km, 5: tc[1]*km, 6: tc[2]*km, 7: 4.0*km},
                                                  bc={'Dirichlet': {10: T_LAB, 7: T_surface}},verbose=True)
    elif dimension == '2D':
        if layers == 2:
            if do_hp == 'REF' or do_hp == 'CONST':
                T = pg.solver.solveFiniteElements(mesh,
                                                  a={0: 1.5, 5: tc[0], 8: 4.0},
                                                  f=force,
                                                  bc={'Dirichlet': {8: T_LAB, 5: T_surface}},verbose=True)
            elif do_hp == 'NONE':
                T = pg.solver.solveFiniteElements(mesh,
                                                  a={0: 1.5, 5: tc[0], 8: 4.0},
                                                  bc={'Dirichlet': {8: T_LAB, 5: T_surface}},verbose=True)
        elif do_sediments and layers == 3:
            if do_hp == 'REF' or do_hp == 'CONST':
                T = pg.solver.solveFiniteElements(mesh,
                                                  a={0: 1.5, 5: tc[0], 6: tc[1], 8: 4.0},
                                                  f=force,
                                                  bc={'Dirichlet': {8: T_LAB, 5: T_surface}},verbose=True)
            elif do_hp == 'NONE':
                T = pg.solver.solveFiniteElements(mesh,
                                                  a={0: 1.5, 5: tc[0], 6: tc[1], 8: 4.0},
                                                  bc={'Dirichlet': {8: T_LAB, 5: T_surface}},verbose=True)
        elif not do_sediments and layers == 3:
            if do_hp == 'REF' or do_hp == 'CONST':
                T = pg.solver.solveFiniteElements(mesh,
                                                  a={0: 1.5, 5: tc[0], 7: tc[1], 8: 4.0},
                                                  f=force,
                                                  bc={'Dirichlet': {8: T_LAB, 5: T_surface}},verbose=True)
            elif do_hp == 'NONE':
                T = pg.solver.solveFiniteElements(mesh,
                                                  a={0: 1.5, 5: tc[0], 7: tc[1], 8: 4.0},
                                                  bc={'Dirichlet': {8: T_LAB, 5: T_surface}},verbose=True)
        elif do_sediments and layers == 4:
            if do_hp == 'REF' or do_hp == 'CONST':
                T = pg.solver.solveFiniteElements(mesh,
                                                  a={0: 1.5, 5: tc[0], 6: tc[1], 7: tc[2], 8: 4.0},
                                                  f=force,
                                                  bc={'Dirichlet': {8: T_LAB, 5: T_surface}},verbose=True)
            elif do_hp == 'NONE':T = pg.solver.solveFiniteElements(mesh,
                                              a={0: 1.5, 5: tc[0], 6: tc[1], 7: tc[2], 8: 4.0},
                                              bc={'Dirichlet': {8: T_LAB, 5: T_surface}},verbose=True)
    else:
        print('Please give correct dimension')
    
    return T

#%% assing_markers_3d
def assign_markers_3d(data_list,
                      mesh,
                      layers,
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
    do_hp :        REF/CONST/NONE - how to deal with heat production
    tc :           list of thermal conductivities
    hp_0 :         which CONST value to take for heat production

    Returns
    -------
    mesh : mesh with assigned markers and parameters
    '''
    
    x_int_t, y_int_t, grid = data_list[0]
    xi_t, yi_t = np.meshgrid(x_int_t/km, y_int_t/km)
    grid_xy_t = np.column_stack((xi_t.ravel(), yi_t.ravel()))
    tree_t = KDTree(grid_xy_t)
    
    if layers == 2:
        x_int_l, y_int_l, grid = data_list[1]
        xi_l, yi_l = np.meshgrid(x_int_l/km, y_int_l/km)
        grid_xy_l = np.column_stack((xi_l.ravel(), yi_l.ravel()))
        tree_l = KDTree(grid_xy_l)
        
        force = np.zeros(mesh.nodeCount())
        for i, node in enumerate(mesh.nodes()): # apply heat production values to nodes
            if do_hp == 'NONE':
                break
            
            x,y,z = node.pos()
            z_topo = compare_to_plane(x,y,data_list[0][2],grid_xy_t,tree_t)
            z_lab = compare_to_plane(x,y,-data_list[1][2],grid_xy_l,tree_l)
            
            idx_A,idy_A = np.argmin(np.abs(np.unique(data_list[0][0])/km-x)),np.argmin(np.abs(np.unique(data_list[0][1])/km-y))
            if z >= z_lab and z <= z_topo:
                if do_hp == 'REF':
                    force[i] = data_list[2][2][idy_A,idx_A]
                elif do_hp == 'CONST':
                    force[i] = hp_0
            else:
                force[i] = 0
            
        for i,cell in enumerate(mesh.cells()): # apply geology markers to cells
            center = cell.center()
            x, y, z = center.x(), center.y(), center.z()
            
            z_topo = compare_to_plane(x,y,data_list[0][2],grid_xy_t,tree_t)
            z_lab = compare_to_plane(x,y,-data_list[1][2],grid_xy_l,tree_l)
            
            if z < z_topo:
                cell.setMarker(4)
            if z < z_lab:
                cell.setMarker(7)
        
    elif do_sediments and layers == 3:
        x_int_l, y_int_l, grid = data_list[2]
        xi_l, yi_l = np.meshgrid(x_int_l/km, y_int_l/km)
        grid_xy_l = np.column_stack((xi_l.ravel(), yi_l.ravel()))
        tree_l = KDTree(grid_xy_l)
        
        force = np.zeros(mesh.nodeCount())
        for i, node in enumerate(mesh.nodes()): # apply heat production values to nodes
            if do_hp == 'NONE':
                break
            
            x,y,z = node.pos()
            z_topo = compare_to_plane(x,y,data_list[0][2],grid_xy_t,tree_t)
            z_sed = compare_to_plane(x,y,(data_list[0][2]-data_list[1][2]),grid_xy_t,tree_t)
            z_lab = compare_to_plane(x,y,-data_list[2][2],grid_xy_l,tree_l)
            
            idx_A,idy_A = np.argmin(np.abs(np.unique(data_list[0][0])/km-x)),np.argmin(np.abs(np.unique(data_list[0][1])/km-y))
            if z >= z_lab and z <= z_sed:
                if do_hp == 'REF':
                    force[i] = data_list[3][2][idy_A,idx_A]
                elif do_hp == 'CONST':
                    force[i] = hp_0
            else:
                force[i] = 0
            
        for i,cell in enumerate(mesh.cells()): # apply geology markers to cells
            center = cell.center()
            x, y, z = center.x(), center.y(), center.z()
            
            z_topo = compare_to_plane(x,y,data_list[0][2],grid_xy_t,tree_t)
            z_sed = compare_to_plane(x,y,(data_list[1][2]-data_list[1][2]),grid_xy_t,tree_t)
            z_lab = compare_to_plane(x,y,-data_list[2][2],grid_xy_l,tree_l)
            
            if z < z_topo:
                cell.setMarker(4)
            if z < z_sed:
                cell.setMarker(5)
            if z < z_lab:
                cell.setMarker(7)
        
    elif not do_sediments and layers == 3:
        x_int_m, y_int_m, grid = data_list[1]
        xi_m, yi_m = np.meshgrid(x_int_m/km, y_int_m/km)
        grid_xy_m = np.column_stack((xi_m.ravel(), yi_m.ravel()))
        tree_m = KDTree(grid_xy_m)
        x_int_l, y_int_l, grid = data_list[2]
        xi_l, yi_l = np.meshgrid(x_int_l/km, y_int_l/km)
        grid_xy_l = np.column_stack((xi_l.ravel(), yi_l.ravel()))
        tree_l = KDTree(grid_xy_l)
        
        force = np.zeros(mesh.nodeCount())
        for i, node in enumerate(mesh.nodes()): # apply heat production values to nodes
            if do_hp == 'NONE':
                break
            
            x,y,z = node.pos()
            z_topo = compare_to_plane(x,y,data_list[0][2],grid_xy_t,tree_t)
            z_moho = compare_to_plane(x,y,-data_list[1][2],grid_xy_m,tree_m)
            z_lab = compare_to_plane(x,y,-data_list[2][2],grid_xy_l,tree_l)
            
            idx_A,idy_A = np.argmin(np.abs(np.unique(data_list[0][0])/km-x)),np.argmin(np.abs(np.unique(data_list[0][1])/km-y))
            if z >= z_moho and z <= z_topo:
                if do_hp == 'REF':
                    force[i] = data_list[3][2][idy_A,idx_A]
                elif do_hp == 'CONST':
                    force[i] = hp_0
            else:
                force[i] = 0
            
        for i,cell in enumerate(mesh.cells()): # apply geology markers to cells
            center = cell.center()
            x, y, z = center.x(), center.y(), center.z()
            
            z_topo = compare_to_plane(x,y,data_list[0][2],grid_xy_t,tree_t)
            z_moho = compare_to_plane(x,y,-data_list[1][2],grid_xy_m,tree_m)
            z_lab = compare_to_plane(x,y,-data_list[2][2],grid_xy_l,tree_l)
            
            if z < z_topo:
                cell.setMarker(4)
            if z < z_moho:
                cell.setMarker(6)
            if z < z_lab:
                cell.setMarker(7)
        
    elif do_sediments and layers == 4:
        x_int_m, y_int_m, grid = data_list[2]
        xi_m, yi_m = np.meshgrid(x_int_m/km, y_int_m/km)
        grid_xy_m = np.column_stack((xi_m.ravel(), yi_m.ravel()))
        tree_m = KDTree(grid_xy_m)
        x_int_l, y_int_l, grid = data_list[3]
        xi_l, yi_l = np.meshgrid(x_int_l/km, y_int_l/km)
        grid_xy_l = np.column_stack((xi_l.ravel(), yi_l.ravel()))
        tree_l = KDTree(grid_xy_l)
    
        force = np.zeros(mesh.nodeCount())
        for i, node in enumerate(mesh.nodes()): # apply heat production values to nodes
            if do_hp == 'NONE':
                break
            
            x,y,z = node.pos()
            z_topo = compare_to_plane(x,y,data_list[0][2],grid_xy_t,tree_t)
            z_sed = compare_to_plane(x,y,(data_list[0][2]-data_list[1][2]),grid_xy_t,tree_t)
            z_moho = compare_to_plane(x,y,-data_list[2][2],grid_xy_m,tree_m)
            z_lab = compare_to_plane(x,y,-data_list[3][2],grid_xy_l,tree_l)
            
            idx_A,idy_A = np.argmin(np.abs(np.unique(data_list[0][0])/km-x)),np.argmin(np.abs(np.unique(data_list[0][1])/km-y))
            if z >= z_moho and z <= z_sed:
                if do_hp == 'REF':
                    force[i] = data_list[4][2][idy_A,idx_A]
                elif do_hp == 'CONST':
                    force[i] = hp_0
            else:
                force[i] = 0
        
        for i,cell in enumerate(mesh.cells()): # apply geology markers to cells
            center = cell.center()
            x, y, z = center.x(), center.y(), center.z()
            
            z_topo = compare_to_plane(x,y,data_list[0][2],grid_xy_t,tree_t)
            z_sed = compare_to_plane(x,y,(data_list[0][2]-data_list[1][2]),grid_xy_t,tree_t)
            z_moho = compare_to_plane(x,y,-data_list[2][2],grid_xy_m,tree_m)
            z_lab = compare_to_plane(x,y,-data_list[3][2],grid_xy_l,tree_l)
            
            if z < z_topo:
                cell.setMarker(4)
            if z < z_sed:
                cell.setMarker(5)
            if z < z_moho:
                cell.setMarker(6)
            if z < z_lab:
                cell.setMarker(7)
    if do_hp == 'NONE':
        return mesh,None
    else:         
        return mesh,force

#%% assign_markers_2d
def assign_markers_2d(data_list,
                      mesh,
                      layers,
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
    
    force = np.zeros(mesh.nodeCount())
    
    if layers == 2:
        for i, node in enumerate(mesh.nodes()): # apply heat production values to nodes
            if do_hp == 'NONE':
                break
            x,y,z = node.pos()
            
            idx_t_0,idx_t_1 = np.argsort(np.abs(data_list[0][0]-x))[0],np.argsort(np.abs(data_list[0][0]-x))[1]
            
            m_t = (data_list[0][1][idx_t_0]-data_list[0][1][idx_t_1])/(data_list[0][0][idx_t_0]-data_list[0][0][idx_t_1])
            y_t = m_t * (x-data_list[0][0][idx_t_1]) + data_list[0][1][idx_t_1] 
            
            idx_l = np.argmin(np.abs(data_list[1][0]-x))
            idx_A = np.argmin(np.abs(data_list[2][0]-x))
            
            if y >= -data_list[1][1][idx_l] and y <= y_t:
                if do_hp == 'REF':
                    force[i] = data_list[2][1][idx_A]
                elif do_hp == 'CONST':
                    force[i] = hp_0
            else:
                force[i] = 0
    
        for cell in mesh.cells(): # apply geology markers to cells
            center = cell.center()
            x, y, z = center.x(), center.y(), center.z()
            
            idx_t_0,idx_t_1 = np.argsort(np.abs(data_list[0][0]-x))[0],np.argsort(np.abs(data_list[0][0]-x))[1]
            
            m_t = (data_list[0][1][idx_t_0]-data_list[0][1][idx_t_1])/(data_list[0][0][idx_t_0]-data_list[0][0][idx_t_1])
            y_t = m_t * (x-data_list[0][0][idx_t_1]) + data_list[0][1][idx_t_1]
            
            idx_l = np.argmin(np.abs(data_list[1][0]-x))
            
            if y < y_t:
                cell.setMarker(5)
            if y < -data_list[1][1][idx_l]:
                cell.setMarker(8)
        
    elif do_sediments and layers == 3:
        for i, node in enumerate(mesh.nodes()): # apply heat production values to nodes
            if do_hp == 'NONE':
                break
            x,y,z = node.pos()
            
            idx_t_0,idx_t_1 = np.argsort(np.abs(data_list[0][0]-x))[0],np.argsort(np.abs(data_list[0][0]-x))[1]
            
            m_s = (data_list[1][1][idx_t_0]-data_list[1][1][idx_t_1])/(data_list[0][0][idx_t_0]-data_list[0][0][idx_t_1])
            y_s = m_s * (x-data_list[0][0][idx_t_1]) + data_list[1][1][idx_t_1] 
            
            m_t = (data_list[0][1][idx_t_0]-data_list[0][1][idx_t_1])/(data_list[0][0][idx_t_0]-data_list[0][0][idx_t_1])
            y_t = m_t * (x-data_list[0][0][idx_t_1]) + data_list[0][1][idx_t_1] 
            
            idx_l = np.argmin(np.abs(data_list[2][0]-x))
            idx_A = np.argmin(np.abs(data_list[3][0]-x))
            
            if y >= -data_list[2][1][idx_l] and y <= y_s:
                if do_hp == 'REF':
                    force[i] = data_list[3][1][idx_A]
                elif do_hp == 'CONST':
                    force[i] = hp_0
            else:
                force[i] = 0
    
        for cell in mesh.cells(): # apply geology markers to cells
            center = cell.center()
            x, y, z = center.x(), center.y(), center.z()
            
            idx_t_0,idx_t_1 = np.argsort(np.abs(data_list[0][0]-x))[0],np.argsort(np.abs(data_list[0][0]-x))[1]
            
            m_t = (data_list[0][1][idx_t_0]-data_list[0][1][idx_t_1])/(data_list[0][0][idx_t_0]-data_list[0][0][idx_t_1])
            m_s = (data_list[1][1][idx_t_0]-data_list[1][1][idx_t_1])/(data_list[0][0][idx_t_0]-data_list[0][0][idx_t_1])
            y_t = m_t * (x-data_list[0][0][idx_t_1]) + data_list[0][1][idx_t_1] 
            y_s = m_s * (x-data_list[0][0][idx_t_1]) + data_list[1][1][idx_t_1] 
            
            idx_l = np.argmin(np.abs(data_list[2][0]-x))
            
            if y < y_t and y > y_s:
                cell.setMarker(5)
            if y < y_s:
                cell.setMarker(6)
            if y < -data_list[2][1][idx_l]:
                cell.setMarker(8)
        
    elif not do_sediments and layers == 3:
        for i, node in enumerate(mesh.nodes()): # apply heat production values to nodes
            if do_hp == 'NONE':
                break
            x,y,z = node.pos()
            
            idx_t_0,idx_t_1 = np.argsort(np.abs(data_list[0][0]-x))[0],np.argsort(np.abs(data_list[0][0]-x))[1]
            
            m_t = (data_list[0][1][idx_t_0]-data_list[0][1][idx_t_1])/(data_list[0][0][idx_t_0]-data_list[0][0][idx_t_1])
            y_t = m_t * (x-data_list[0][0][idx_t_1]) + data_list[0][1][idx_t_1] 
            
            idx_bA = np.argmin(np.abs(data_list[0][0]-x))
            idx_l = np.argmin(np.abs(data_list[2][0]-x))
            idx_A = np.argmin(np.abs(data_list[3][0]-x))
            
            if y >= -data_list[1][1][idx_bA] and y <= y_t:
                if do_hp == 'REF':
                    force[i] = data_list[3][1][idx_A]
                elif do_hp == 'CONST':
                    force[i] = hp_0
            else:
                force[i] = 0
    
        for cell in mesh.cells(): # apply geology markers to cells
            center = cell.center()
            x, y, z = center.x(), center.y(), center.z()
            
            idx_t_0,idx_t_1 = np.argsort(np.abs(data_list[0][0]-x))[0],np.argsort(np.abs(data_list[0][0]-x))[1]
            
            m_t = (data_list[0][1][idx_t_0]-data_list[0][1][idx_t_1])/(data_list[0][0][idx_t_0]-data_list[0][0][idx_t_1])
            y_t = m_t * (x-data_list[0][0][idx_t_1]) + data_list[0][1][idx_t_1]
            
            idx_m = np.argmin(np.abs(data_list[1][0]-x))
            idx_l = np.argmin(np.abs(data_list[2][0]-x))
            
            if y < y_t:
                cell.setMarker(5)
            if y < -data_list[1][1][idx_m]:
                cell.setMarker(7)
            if y < -data_list[2][1][idx_l]:
                cell.setMarker(8)
        
    elif do_sediments and layers == 4:
        for i, node in enumerate(mesh.nodes()): # apply heat production values to nodes
            if do_hp == 'NONE':
                break
            x,y,z = node.pos()
            
            idx_t_0,idx_t_1 = np.argsort(np.abs(data_list[0][0]-x))[0],np.argsort(np.abs(data_list[0][0]-x))[1]
            
            m_s = (data_list[1][1][idx_t_0]-data_list[1][1][idx_t_1])/(data_list[0][0][idx_t_0]-data_list[0][0][idx_t_1])
            y_s = m_s * (x-data_list[0][0][idx_t_1]) + data_list[1][1][idx_t_1] 
            
            m_t = (data_list[0][1][idx_t_0]-data_list[0][1][idx_t_1])/(data_list[0][0][idx_t_0]-data_list[0][0][idx_t_1])
            y_t = m_t * (x-data_list[0][0][idx_t_1]) + data_list[0][1][idx_t_1] 
            
            idx_bA = np.argmin(np.abs(data_list[2][0]-x))
            idx_l = np.argmin(np.abs(data_list[3][0]-x))
            idx_A = np.argmin(np.abs(data_list[4][0]-x))
            
            if y >= -data_list[2][1][idx_bA] and y <= y_s:
                if do_hp == 'REF':
                    force[i] = data_list[4][1][idx_A]
                elif do_hp == 'CONST':
                    force[i] = hp_0
            else:
                force[i] = 0
    
        for cell in mesh.cells(): # apply geology markers to cells
            center = cell.center()
            x, y, z = center.x(), center.y(), center.z()
            
            idx_t_0,idx_t_1 = np.argsort(np.abs(data_list[0][0]-x))[0],np.argsort(np.abs(data_list[0][0]-x))[1]
            
            m_t = (data_list[0][1][idx_t_0]-data_list[0][1][idx_t_1])/(data_list[0][0][idx_t_0]-data_list[0][0][idx_t_1])
            m_s = (data_list[1][1][idx_t_0]-data_list[1][1][idx_t_1])/(data_list[0][0][idx_t_0]-data_list[0][0][idx_t_1])
            y_t = m_t * (x-data_list[0][0][idx_t_1]) + data_list[0][1][idx_t_1] 
            y_s = m_s * (x-data_list[0][0][idx_t_1]) + data_list[1][1][idx_t_1] 
            
            idx_m = np.argmin(np.abs(data_list[2][0]-x))
            idx_l = np.argmin(np.abs(data_list[3][0]-x))
            
            if y < y_t and y > y_s:
                cell.setMarker(5)
            if y < y_s:
                cell.setMarker(6)
            if y < -data_list[2][1][idx_m]:
                cell.setMarker(7)
            if y < -data_list[3][1][idx_l]:
                cell.setMarker(8)
    
    if do_hp == 'NONE':
        return mesh,None
    else:         
        return mesh,force

#%% build_world_2d
def build_world_2d(data_list          ,
                  profile_coord       ,
                  layers       ,
                  area = None         ,
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
    if layers == 2:
        lab_line = mt.createPolygon([[a,b] for a,b in zip(data_list[1][0],-data_list[1][1])],marker = 8,boundaryMarker = 8,isClosed=False)
        line_list.append(lab_line)
    elif do_sediments and layers == 3:
        sed_line = mt.createPolygon([[a,b] for a,b in zip(data_list[1][0],data_list[1][1])],marker = 6,boundaryMarker = 6,isClosed=False) 
        lab_line = mt.createPolygon([[a,b] for a,b in zip(data_list[2][0],-data_list[2][1])],marker = 8,boundaryMarker = 8,isClosed=False)
        line_list.append(sed_line)
        line_list.append(lab_line)
    elif not do_sediments and layers == 3:
        moho_line = mt.createPolygon([[a,b] for a,b in zip(data_list[1][0],-data_list[1][1])],marker = 7,boundaryMarker = 6,isClosed=False) 
        lab_line = mt.createPolygon([[a,b] for a,b in zip(data_list[2][0],-data_list[2][1])],marker = 8,boundaryMarker = 8,isClosed=False)
        line_list.append(moho_line)
        line_list.append(lab_line)
    elif do_sediments and layers == 4:
        sed_line = mt.createPolygon([[a,b] for a,b in zip(data_list[0][0],data_list[1][1])],marker = 6,boundaryMarker = 6,isClosed=False)
        moho_line = mt.createPolygon([[a,b] for a,b in zip(data_list[2][0],-data_list[2][1])],marker = 7,boundaryMarker = 7,isClosed=False)
        lab_line = mt.createPolygon([[a,b] for a,b in zip(data_list[3][0],-data_list[3][1])],marker = 8,boundaryMarker = 8,isClosed=False)
        line_list.append(sed_line)
        line_list.append(moho_line)
        line_list.append(lab_line)
    line_list.append(topo_line)
    
    for i in line_list:
        world = world + i
    
    if area != None:
        mesh = mt.createMesh(world,quality=34,area=area)
    else:
        mesh = mt.createMesh(world,quality=34,area=profile_length*4)
    
    return mesh, line_list

#%% build_world_3d
def build_world_3d(data_list          ,
                  area_coords         ,
                  layers       ,
                  area = None         ,
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
        for boundary in surface_topo.boundaries():
            boundary.setMarker(7)
        layer_list.append(surface_topo)
    
    if layers == 2:
        mesh_lab = mt.createMesh2D(data_list[1][0]/km,data_list[1][1]/km)
        surface_lab = mt.createSurface(mesh_lab)
        surface_lab = fix_surface_height(surface_lab, data_list[1][0]/km,data_list[1][1]/km, -data_list[1][2]/km)
        layer_list.append(surface_lab)
        for boundary in surface_lab.boundaries():
            boundary.setMarker(10)
        layer_list.append(surface_lab)
    elif layers == 3 and not do_sediments:
        mesh_moho = mt.createMesh2D(data_list[1][0]/km,data_list[1][1]/km)
        surface_moho = mt.createSurface(mesh_moho)
        surface_moho = fix_surface_height(surface_moho, data_list[1][0]/km,data_list[1][1]/km, -data_list[1][2]/km)
        for boundary in surface_moho.boundaries():
            boundary.setMarker(9)
        layer_list.append(surface_moho)
    
        mesh_lab = mt.createMesh2D(data_list[2][0]/km,data_list[2][1]/km)
        surface_lab = mt.createSurface(mesh_lab)
        surface_lab = fix_surface_height(surface_lab, data_list[2][0]/km,data_list[2][1]/km, -data_list[2][2]/km)
        for boundary in surface_lab.boundaries():
            boundary.setMarker(10)
        layer_list.append(surface_lab)
    elif layers == 3 and do_sediments:
        mesh_lab = mt.createMesh2D(data_list[2][0]/km,data_list[2][1]/km)
        surface_lab = mt.createSurface(mesh_lab)
        surface_lab = fix_surface_height(surface_lab, data_list[2][0]/km,data_list[2][1]/km, -data_list[2][2]/km)
        for boundary in surface_lab.boundaries():
            boundary.setMarker(10)
        layer_list.append(surface_lab)
    elif layers == 4:
        mesh_moho = mt.createMesh2D(data_list[2][0]/km,data_list[2][1]/km)
        surface_moho = mt.createSurface(mesh_moho)
        surface_moho = fix_surface_height(surface_moho, data_list[2][0]/km,data_list[2][1]/km, -data_list[2][2]/km)
        for boundary in surface_moho.boundaries():
            boundary.setMarker(9)
        layer_list.append(surface_moho)
    
        mesh_lab = mt.createMesh2D(data_list[3][0]/km,data_list[3][1]/km)
        surface_lab = mt.createSurface(mesh_lab)
        surface_lab = fix_surface_height(surface_lab, data_list[3][0]/km,data_list[3][1]/km, -data_list[3][2]/km)
        for boundary in surface_lab.boundaries():
            boundary.setMarker(10)
        layer_list.append(surface_lab)
    else:
        print('something went wrong with the input')
    
    geometry = world
    for i in layer_list:
        geometry = geometry + i
    
    if area == None:
        mesh = mt.createMesh(geometry,quality=34)
    else:
        mesh = mt.createMesh(geometry,quality=34,area=area)

    return mesh,layer_list

#%% cut_and_interpolate_3d    
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

    if layer == 'CONST' and i==layers:
        x_int,y_int = np.arange(area_coords[0],area_coords[1]+resolution[0],resolution[0]),np.arange(area_coords[3],area_coords[2]+resolution[0],resolution[0])
        xi_int,yi_int = np.meshgrid(x_int,y_int)
        grid = np.full(np.shape(yi_int),hp_0)
        
        return x_int,y_int,grid

    data,xi,yi = layer["data"],layer["x"],layer["y"]
    
    xi_s,yi_s,data_s = in_area_s(area_coords,xi,yi,data)
    x_int,y_int = np.arange(area_coords[0],area_coords[1]+resolution[i],resolution[i]),np.arange(area_coords[3],area_coords[2]+resolution[i],resolution[i])
    xi_int,yi_int = np.meshgrid(x_int,y_int)

    if do_sediments and i==1:
        interp = RegularGridInterpolator((np.unique(xi_s),np.unique(yi_s)),data_s)
        x_int_t,y_int_t = np.arange(area_coords[0],area_coords[1]+resolution[0],resolution[0]),np.arange(area_coords[3],area_coords[2]+resolution[0],resolution[0])
        xi_int_t,yi_int_t = np.meshgrid(x_int_t,y_int_t)
        grid = interp((xi_int_t,yi_int_t))
    elif (do_hp=='REF' or do_hp=='CONST') and layers == i:
        #interp = RegularGridInterpolator((x_int,y_int),data_s)
        interp = RegularGridInterpolator((np.unique(xi_s),np.unique(yi_s)),data_s)
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

#%% cut_and_interpolate_2d
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
    
    if layer == 'CONST' and i==layers:
        x_profile,y_profile = (profile_coord[0],profile_coord[1])
        x_int, y_int, dist_prof = define_profile(x_profile, y_profile, resolution[0])
        dist_prof = np.round(dist_prof)
        profile = np.full(np.shape(dist_prof),hp_0)
        
        return dist_prof,profile
    
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
    elif do_hp == 'REF' and i==layers:
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

#%% in_area_s
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

#%% fix_surface_height
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

#%% compare_to_plane
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
    height_data = data.flatten()[indices]/km #np.transpose(data).flatten('F')[indices]/km

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

#%% define_profile
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

#%% create_sed_interface
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
                triangle.createPolygonFace(triangle.nodes(surf),marker=7)
                triangles_topo.append(triangle)
                
                triangle = pg.Mesh(3,isGeometry=True)
                triangle.createNode(x_int_t[i+1]/km,y_int_t[j]/km,grid_topo[i+1,j]/km,marker=4)
                triangle.createNode(x_int_t[i]/km,y_int_t[j+1]/km,grid_topo[i,j+1]/km,marker=4)
                triangle.createNode(x_int_t[i+1]/km,y_int_t[j+1]/km,grid_topo[i+1,j+1]/km,marker=4)
                surf = [0,1,2]
                triangle.createPolygonFace(triangle.nodes(surf),marker=7)
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
                    triangle.createPolygonFace(triangle.nodes(surf),marker=8)
                    triangles_sed.append(triangle)
                    
                    # topo
                    triangle = pg.Mesh(3,isGeometry=True)
                    triangle.createNode(x_int_t[i]/km,y_int_t[j]/km,grid_topo[i,j]/km,marker=4)
                    triangle.createNode(x_int_t[i+1]/km,y_int_t[j]/km,grid_topo[i+1,j]/km,marker=9)
                    triangle.createNode(x_int_t[i]/km,y_int_t[j+1]/km,grid_topo[i,j+1]/km,marker=9)
                    surf = [0,1,2]
                    triangle.createPolygonFace(triangle.nodes(surf),marker=7)
                    triangles_topo.append(triangle)
                    
                    triangle = pg.Mesh(3,isGeometry=True)
                    triangle.createNode(x_int_t[i+1]/km,y_int_t[j]/km,grid_topo[i+1,j]/km,marker=9)
                    triangle.createNode(x_int_t[i]/km,y_int_t[j+1]/km,grid_topo[i,j+1]/km,marker=9)
                    triangle.createNode(x_int_t[i+1]/km,y_int_t[j+1]/km,grid_topo[i+1,j+1]/km,marker=9)
                    surf = [0,1,2]
                    triangle.createPolygonFace(triangle.nodes(surf),marker=7)
                    triangles_topo.append(triangle)
                    continue
                elif grid_sed[i+1,j] != 0:
                    # sed
                    triangle = pg.Mesh(3,isGeometry=True)
                    triangle.createNode([x_int_t[i+1]/km,y_int_t[j]/km,(grid_topo[i+1,j]-grid_sed[i+1,j])/km],marker=5)
                    triangle.createNode([x_int_t[i]/km,y_int_t[j]/km,(grid_topo[i,j]-grid_sed[i,j])/km],marker=9)
                    triangle.createNode([x_int_t[i+1]/km,y_int_t[j+1]/km,(grid_topo[i+1,j+1]-grid_sed[i+1,j+1])/km],marker=9)
                    surf = [0,1,2]
                    triangle.createPolygonFace(triangle.nodes(surf),marker=8)
                    triangles_sed.append(triangle)
                    
                    # topo
                    triangle = pg.Mesh(3,isGeometry=True)
                    triangle.createNode(x_int_t[i]/km,y_int_t[j]/km,grid_topo[i,j]/km,marker=9)
                    triangle.createNode(x_int_t[i+1]/km,y_int_t[j]/km,grid_topo[i+1,j]/km,marker=4)
                    triangle.createNode(x_int_t[i+1]/km,y_int_t[j+1]/km,grid_topo[i+1,j+1]/km,marker=9)
                    surf = [0,1,2]
                    triangle.createPolygonFace(triangle.nodes(surf),marker=7)
                    triangles_topo.append(triangle)
                    
                    triangle = pg.Mesh(3,isGeometry=True)
                    triangle.createNode(x_int_t[i]/km,y_int_t[j]/km,grid_topo[i,j]/km,marker=9)
                    triangle.createNode(x_int_t[i]/km,y_int_t[j+1]/km,grid_topo[i,j+1]/km,marker=9)
                    triangle.createNode(x_int_t[i+1]/km,y_int_t[j+1]/km,grid_topo[i+1,j+1]/km,marker=9)
                    surf = [0,1,2]
                    triangle.createPolygonFace(triangle.nodes(surf),marker=7)
                    triangles_topo.append(triangle)
                    continue
                elif grid_sed[i,j+1] != 0:
                    # sed
                    triangle = pg.Mesh(3,isGeometry=True)
                    triangle.createNode([x_int_t[i]/km,y_int_t[j+1]/km,(grid_topo[i,j+1]-grid_sed[i,j+1])/km],marker=5)
                    triangle.createNode([x_int_t[i]/km,y_int_t[j]/km,(grid_topo[i,j]-grid_sed[i,j])/km],marker=9)
                    triangle.createNode([x_int_t[i+1]/km,y_int_t[j+1]/km,(grid_topo[i+1,j+1]-grid_sed[i+1,j+1])/km],marker=9)
                    surf = [0,1,2]
                    triangle.createPolygonFace(triangle.nodes(surf),marker=8)
                    triangles_sed.append(triangle)
                    
                    # topo
                    triangle = pg.Mesh(3,isGeometry=True)
                    triangle.createNode(x_int_t[i]/km,y_int_t[j+1]/km,grid_topo[i,j+1]/km,marker=4)
                    triangle.createNode(x_int_t[i]/km,y_int_t[j]/km,grid_topo[i,j]/km,marker=9)
                    triangle.createNode(x_int_t[i+1]/km,y_int_t[j+1]/km,grid_topo[i+1,j+1]/km,marker=9)
                    surf = [0,1,2]
                    triangle.createPolygonFace(triangle.nodes(surf),marker=7)
                    triangles_topo.append(triangle)
                    
                    triangle = pg.Mesh(3,isGeometry=True)
                    triangle.createNode(x_int_t[i+1]/km,y_int_t[j]/km,grid_topo[i+1,j]/km,marker=9)
                    triangle.createNode(x_int_t[i]/km,y_int_t[j]/km,grid_topo[i,j]/km,marker=9)
                    triangle.createNode(x_int_t[i+1]/km,y_int_t[j+1]/km,grid_topo[i+1,j+1]/km,marker=9)
                    surf = [0,1,2]
                    triangle.createPolygonFace(triangle.nodes(surf),marker=7)
                    triangles_topo.append(triangle)
                    continue
                elif grid_sed[i+1,j+1] != 0:
                    # sed
                    triangle = pg.Mesh(3,isGeometry=True)
                    triangle.createNode([x_int_t[i+1]/km,y_int_t[j+1]/km,(grid_topo[i+1,j+1]-grid_sed[i+1,j+1])/km],marker=5)
                    triangle.createNode([x_int_t[i+1]/km,y_int_t[j]/km,(grid_topo[i+1,j]-grid_sed[i+1,j])/km],marker=9)
                    triangle.createNode([x_int_t[i]/km,y_int_t[j+1]/km,(grid_topo[i,j+1]-grid_sed[i,j+1])/km],marker=9)
                    surf = [0,1,2]
                    triangle.createPolygonFace(triangle.nodes(surf),marker=8)
                    triangles_sed.append(triangle)
                    
                    # topo
                    triangle = pg.Mesh(3,isGeometry=True)
                    triangle.createNode(x_int_t[i]/km,y_int_t[j]/km,grid_topo[i,j]/km,marker=9)
                    triangle.createNode(x_int_t[i+1]/km,y_int_t[j]/km,grid_topo[i+1,j]/km,marker=9)
                    triangle.createNode(x_int_t[i]/km,y_int_t[j+1]/km,grid_topo[i,j+1]/km,marker=9)
                    surf = [0,1,2]
                    triangle.createPolygonFace(triangle.nodes(surf),marker=7)
                    triangles_topo.append(triangle)
                    
                    triangle = pg.Mesh(3,isGeometry=True)
                    triangle.createNode(x_int_t[i+1]/km,y_int_t[j]/km,grid_topo[i+1,j]/km,marker=9)
                    triangle.createNode(x_int_t[i]/km,y_int_t[j+1]/km,grid_topo[i,j+1]/km,marker=9)
                    triangle.createNode(x_int_t[i+1]/km,y_int_t[j+1]/km,grid_topo[i+1,j+1]/km,marker=4)
                    surf = [0,1,2]
                    triangle.createPolygonFace(triangle.nodes(surf),marker=7)
                    triangles_topo.append(triangle)
                    continue
            else:
                # sed
                triangle = pg.Mesh(3,isGeometry=True)
                triangle.createNode(x_int_t[i]/km,y_int_t[j]/km,(grid_topo[i,j]-grid_sed[i,j])/km,marker=5)
                triangle.createNode(x_int_t[i+1]/km,y_int_t[j]/km,(grid_topo[i+1,j]-grid_sed[i+1,j])/km,marker=5)
                triangle.createNode(x_int_t[i]/km,y_int_t[j+1]/km,(grid_topo[i,j+1]-grid_sed[i,j+1])/km,marker=5)
                surf = [0,1,2]
                triangle.createPolygonFace(triangle.nodes(surf),marker=8)
                triangles_sed.append(triangle)
                
                triangle = pg.Mesh(3,isGeometry=True)
                triangle.createNode(x_int_t[i+1]/km,y_int_t[j]/km,(grid_topo[i+1,j]-grid_sed[i+1,j])/km,marker=5)
                triangle.createNode(x_int_t[i]/km,y_int_t[j+1]/km,(grid_topo[i,j+1]-grid_sed[i,j+1])/km,marker=5)
                triangle.createNode(x_int_t[i+1]/km,y_int_t[j+1]/km,(grid_topo[i+1,j+1]-grid_sed[i+1,j+1])/km,marker=5)
                surf = [0,1,2]
                triangle.createPolygonFace(triangle.nodes(surf),marker=8)
                triangles_sed.append(triangle)
                
                # topo
                triangle = pg.Mesh(3,isGeometry=True)
                triangle.createNode(x_int_t[i]/km,y_int_t[j]/km,grid_topo[i,j]/km,marker=4)
                triangle.createNode(x_int_t[i+1]/km,y_int_t[j]/km,grid_topo[i+1,j]/km,marker=4)
                triangle.createNode(x_int_t[i]/km,y_int_t[j+1]/km,grid_topo[i,j+1]/km,marker=4)
                surf = [0,1,2]
                triangle.createPolygonFace(triangle.nodes(surf),marker=7)
                triangles_topo.append(triangle)
                
                triangle = pg.Mesh(3,isGeometry=True)
                triangle.createNode(x_int_t[i+1]/km,y_int_t[j]/km,grid_topo[i+1,j]/km,marker=4)
                triangle.createNode(x_int_t[i]/km,y_int_t[j+1]/km,grid_topo[i,j+1]/km,marker=4)
                triangle.createNode(x_int_t[i+1]/km,y_int_t[j+1]/km,grid_topo[i+1,j+1]/km,marker=4)
                surf = [0,1,2]
                triangle.createPolygonFace(triangle.nodes(surf),marker=7)
                triangles_topo.append(triangle)
                continue
    
    # merge all
    surface_topo = mt.mergePLC(triangles_topo)
    surface_sed = mt.mergePLC(triangles_sed)
    return surface_topo, surface_sed
