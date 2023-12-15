import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cpf
import cartopy.io.img_tiles as cimgt
from cartopy.io.shapereader import Reader
from cartopy.io.shapereader import natural_earth
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


def mapa_base(llat, llon):
    """
    Mapa base para graficar las variables
    """
    l_lat = llat
    l_lon = np.array(llon) % 360  #Pasamos lon en [-180, 180] a [0, 360]
    # Trabajamos con el SHAPEFILE de IGN para provincias
    shp = Reader(natural_earth(resolution='10m', category='cultural',
                               name='admin_1_states_provinces_lines'))
    countries = shp.records()
    

    # Comenzamos la Figura
    fig = plt.figure(figsize=(6, 8))
    proj_lcc = ccrs.PlateCarree()
    ax = plt.axes(projection=proj_lcc)
    ax.coastlines(resolution='10m', zorder=1)
    ax.add_feature(cpf.BORDERS, linestyle='-', zorder=1)
    for country in countries:
        print(country)
        exit()
        if country.attributes['adm0_name'] == 'Argentina':
            ax.add_geometries( [country.geometry], ccrs.PlateCarree(),
                                edgecolor='black', facecolor='none',
                                linewidth=0.7, zorder=1 )
    # Colocamos reticula personalizada
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=0.4, color='gray', alpha=0.7, linestyle=':', zorder=1)
    gl.top_labels = False
    gl.right_labels = False
    gl.bottom_labels = False
    gl.left_labels = False
    #gl.xlocator = mticker.FixedLocator(np.linspace(llon[0], llon[1], 5))
    #gl.ylocator = mticker.FixedLocator(np.linspace(l_lat[0], l_lat[1], 9))
    #gl.xformatter = LONGITUDE_FORMATTER
    #gl.yformatter = LATITUDE_FORMATTER
    # Extension del mapa
    ax.set_extent([l_lon[0], l_lon[1], l_lat[0], l_lat[1]], crs=proj_lcc)
    # Posicion del eje (desplazamos un poco a la izquierda y más abajo)
    pos1 = ax.get_position() # get the original position
    pos2 = [pos1.x0 - 0.02, pos1.y0 - 0.06,  pos1.width*1.14, pos1.height*1.22]
    ax.set_position(pos2) # set a new position

    return fig, ax


def colores_tmin():
    
    from matplotlib import colors as c
    # Preparamos los intervalos y colores a usar
    cMap = c.ListedColormap(['#A880C1', '#D1D5F0', '#ffffff'])
    bounds = np.array([-100., 0., 3., 100.])
    norm = c.BoundaryNorm(boundaries=bounds, ncolors=3)

    return cMap, bounds, norm

def legend_temp(ax, x, y, clr, txt, proj):
    """
    Define legends as square for specific values.
    x, y, clr  and txt MUST be arrays of same len to work.
    the proj value, comes from Cartopy projection used in ax axis.
    """
    for i in np.arange(0, len(x)):
        ax.add_patch(mpatches.Rectangle(xy=[x[i], y[i]], width=1.8, height=1.5,
                                         facecolor=clr[i],
                                         transform=proj))
        ax.text(x[i] - -2.2, y[i] - -0.2, txt[i], horizontalalignment='left',
                 fontweight='bold', fontsize=14, transform=ccrs.Geodetic())


def mapa_temperatura_minima_superficial(x, y, LSTmin, l_lat, l_lon, nfig):
    fig, ax = mapa_base(l_lat, l_lon)
    cMap, bounds, norm = colores_tmin()
    ax.pcolormesh(x, y, LSTmin, transform=ccrs.PlateCarree(), cmap=cMap, norm=norm, zorder=0)
    legend_temp(ax, [-60, -60], [-45, -47], ['#D1D5F0', '#A880C1'],
                ['0 - 3 ' + u'\u2103', '< 0' + u'\u2103'],
                ccrs.PlateCarree())
    #ax.set_title('Temperatura mínima para: ' + fecha_tiff)
    #ax.add_image(stamen_terrain, 5) #c=df['prcp24'],cmap='brg'
    #nfig = 'fig_test_LST_GOES' + fecha_nfig + '.jpg'
    plt.savefig(nfig, bbox_inches='tight', dpi=100)