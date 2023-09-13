from osgeo import gdal
from osgeo import ogr
from osgeo import osr
from affine import Affine
import numpy as np
import pandas as pd
import datetime as dt

import sys
sys.path.append('./lib')
from lst_class import lst_class


import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cpf
import cartopy.io.img_tiles as cimgt
from cartopy.io.shapereader import Reader
from cartopy.io.shapereader import natural_earth
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches

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
        if country.attributes['adm0_name'] == 'Argentina':
            ax.add_geometries( [country.geometry], ccrs.PlateCarree(),
                                edgecolor='black', facecolor='none',
                                linewidth=0.7, zorder=1 )
    # Colocamos reticula personalizada
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=0.4, color='gray', alpha=0.7, linestyle=':', zorder=1)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlocator = mticker.FixedLocator(np.linspace(llon[0], llon[1], 5))
    gl.ylocator = mticker.FixedLocator(np.linspace(l_lat[0], l_lat[1], 9))
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    # Extension del mapa
    ax.set_extent([l_lon[0], l_lon[1], l_lat[0], l_lat[1]], crs=proj_lcc)
    # Posicion del eje (desplazamos un poco a la izquierda y más abajo)
    pos1 = ax.get_position() # get the original position
    pos2 = [pos1.x0 - 0.02, pos1.y0 - 0.06,  pos1.width*1.14, pos1.height*1.22]
    ax.set_position(pos2) # set a new position

    return fig, ax


#Objeto para leer IMERG
a = lst_class('20230720')
# Lineas para leer OBS
#carpeta = '../../data_test/'
#archivo = 'precip24h_emas+conv_2020-03-10.txt'
#df = pd.read_csv(carpeta+archivo, sep='\\s+', header=0)
#list_est = '../../obs_tools/lista_est/estaciones_referencia.txt'
#df = pd.read_csv(list_est,delimiter=r"\s+")
#print(df.dtypes)

#stamen_terrain = cimgt.Stamen('terrain')

l_lat = [-60., -14.]
l_lon = [-89., -43.]

fig, ax = mapa_base(l_lat, l_lon)
ax.pcolormesh(a.lon, a.lat, a.lst_min, transform=ccrs.PlateCarree())
#ax.add_image(stamen_terrain, 5) #c=df['prcp24'],cmap='brg'
plt.savefig('fig_test_LST_GOES.jpg', dpi=100)

