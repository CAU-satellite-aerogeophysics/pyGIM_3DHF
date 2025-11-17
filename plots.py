"""
Georg HÃ¼ttner - November 2025

plotting functions to not have to import pygmli in main script
"""

import matplotlib.pyplot as plt
import pygimli as pg

def plot_3D(data_list,
            ghf
            ):
    '''
    plot the 3D GHF results as pcolormesh

    Parameters
    ----------
    data_list :    data input list of [xi, yi, grid]
    ghf :          3D GHF

    Returns
    -------
    Plot
    '''
    
    fig,ax = plt.subplots(figsize=(5,5),dpi=150)
    im = ax.pcolormesh(data_list[0][0][1:-1]/1000,data_list[0][1][1:-1]/1000,ghf[1:-1,1:-1],cmap="afmhot")
    ax.set_xlabel('Easting in km')
    ax.set_ylabel('Northing in km')
    # ax2.set_ylabel('y in m')
    fig.colorbar(im,ax=ax,label='GHF in mW/m$^2$',location='right',shrink=0.65)
    ax.set_aspect('equal')
    ax.tick_params(top=False, labeltop=False, bottom=True, labelbottom=True, left=True, labelleft=True)


def plot_2D(data_list,
            ghf
            ):
    '''
    plot the 2D GHF results as lineplot

    Parameters
    ----------
    data_list : data
    ghf :       2D GHF

    Returns
    -------
    Plot
    '''
    
    fig, ax = plt.subplots(figsize=(7,3),dpi=150)
    ax.plot(data_list[0][0]/1000,ghf,label="1D") 
    ax.set_xlabel("Profile length in km")
    ax.set_ylabel('GHF in mW/m$^2$')
    ax.grid(alpha=0.5)

def plot_2D_mesh(mesh,
                 T=None
                 ):
    '''
    plot the 2D mesh, either with markers or with temperature results
    
    Parameters
    ----------
    mesh :  pygimli mesh from build_world
    T :     Temperature on mesh

    Returns
    -------
    Plot
    '''
    
    if T==None:
        fig, ax = plt.subplots(figsize=(7,7),dpi=300)
        pg.show(mesh,showMesh=True,markers=True,ax=ax)
        # ax.set_aspect('equal','box')
        # ax.set_xlim([0,20000])
        ax.set_xlabel('x in m')
        ax.set_ylabel('depth in m')
    else:
        fig, ax = plt.subplots(figsize=(7,7),dpi=300)
        pg.show(mesh,data=T,showMesh=True, label='Temperature in °C', cMap="inferno",nCols=265,ax=ax)
        # ax.set_aspect('equal','box')
        # ax.set_xlim([0,20000])
        ax.set_xlabel('x in m')
        ax.set_ylabel('depth in m')

def plot_3D_mesh(mesh,
                 T=None
                 ):
    '''
    plot the 3D mesh with pyvista, either with markers or with temperature results

    Parameters
    ----------
    mesh : pygilmi mesh from build_world
    T : Temperatures on mesh
    
    Returns
    -------
    Plot
    '''
    
    if T==None:
        pg.show(mesh,showMesh=True,alpha=0.7)
    else:
        pg.show(mesh,data=T,showMesh=True, label='Temperature in °C', cMap="inferno")
