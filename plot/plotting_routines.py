#----------------------------------------------------------------------------------------------
# Module for plotting routines
#----------------------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, cm, colors, colorbar
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import math

griddir ='../grid/' # must exist
nbfaces = 6
rad2deg = 180/np.pi

#----------------------------------------------------------------------------------------------
# Plot the scalar field q given in a A grid.
#----------------------------------------------------------------------------------------------
def plot_scalarfield(q, map_projection, title, filename, colormap, qmin, qmax,\
    dpi, figformat):

    #Get grid info
    N = np.shape(q)[0]
    lons = np.zeros((N+1,N+1,6))
    lats = np.zeros((N+1,N+1,6))

    for p in range(0, nbfaces):
        input_grid_bin_lon  = griddir+'grid_N'+str(N)+'_face'+str(p+1)+'_lon.dat'    
        input_grid_bin_lat  = griddir+'grid_N'+str(N)+'_face'+str(p+1)+'_lat.dat'

        x = open(input_grid_bin_lon, 'rb')
        x = np.fromfile(x, dtype=np.float64)
        x = np.reshape(x, (N+1,N+1))*rad2deg

        y = open(input_grid_bin_lat, 'rb')
        y = np.fromfile(y, dtype=np.float64)
        y = np.reshape(y, (N+1,N+1))*rad2deg
        lons[:,:,p] = x
        lats[:,:,p] = y


    # map projection
    if map_projection == "mercator":
        plt.figure(figsize=(1832/dpi, 977/dpi), dpi=dpi)
        plateCr = ccrs.PlateCarree()        
    elif map_projection == "sphere":
        plt.figure(figsize=(800/dpi, 800/dpi), dpi=dpi) 
        plateCr = ccrs.Orthographic(central_longitude=0.0, central_latitude=0.0)
        plateCr = ccrs.Orthographic(central_longitude=0.25*180, central_latitude=180/6.0)


    projection=ccrs.PlateCarree(central_longitude=0)
    plateCr._threshold = plateCr._threshold/10.

    ax = plt.axes(projection=plateCr)
    ax.set_global()
    ax.stock_img()
    ax.gridlines(draw_labels=True)


    # Color of each cubed panel
    colors = ('black','black','black','black','black','black')
    colors = ('gray','gray','gray','gray','gray','gray')
    # plot for each tile
    for p in range(0,nbfaces):
        # Get grid
        lon = lons[:,:,p]
        lat = lats[:,:,p]

        # Plot cube edges
        A_lon, A_lat = lon[0,0], lat[0,0]
        B_lon, B_lat = lon[N, 0], lat[N, 0]
        C_lon, C_lat = lon[N, N], lat[N, N]
        D_lon, D_lat = lon[0, N], lat[0, N]
        lw = 0.1
        plt.rcParams["axes.axisbelow"] = True

        ax.plot([A_lon, B_lon], [A_lat, B_lat], '-.', linewidth=lw, color=colors[p], transform=ccrs.Geodetic(), zorder=11)
        ax.plot([B_lon, C_lon], [B_lat, C_lat], '-.', linewidth=lw, color=colors[p], transform=ccrs.Geodetic(), zorder=11)
        ax.plot([C_lon, D_lon], [C_lat, D_lat], '-.', linewidth=lw, color=colors[p], transform=ccrs.Geodetic(), zorder=11)
        ax.plot([D_lon, A_lon], [D_lat, A_lat], '-.', linewidth=lw, color=colors[p], transform=ccrs.Geodetic(), zorder=11)

        # Plot scalar field
        im = ax.pcolormesh(lon, lat, q[:,:,p], alpha=1, transform=ccrs.PlateCarree(), \
        zorder=10, vmin = qmin, vmax=qmax,  cmap=colormap)

    plt.title(title)

    # Plot colorbar
    cax,kw = colorbar.make_axes(ax,orientation='vertical' , fraction=0.046, pad=0.04, shrink=0.8, format='%.1e')
    cb=plt.colorbar(im, cax=cax, extend='both',**kw)
    cb.ax.tick_params(labelsize=8)
    plt.savefig(filename+'.'+figformat, format=figformat)
    print('plotted ', filename)
    plt.close()
