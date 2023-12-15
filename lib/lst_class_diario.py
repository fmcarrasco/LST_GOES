import os
import numpy as np
import numpy.ma as ma
import pandas as pd
from pyproj import Proj

from netCDF4 import Dataset

from datetime import datetime, timedelta

from warnings import filterwarnings
filterwarnings(action='ignore', category=DeprecationWarning, message='`np.bool` is a deprecated alias')


def getGeoTransform(extent, nlines, ncols):
    resx = (extent[2] - extent[0]) / ncols
    resy = (extent[3] - extent[1]) / nlines
    return [extent[0], resx, 0, extent[3] , 0, -resy]


class lst_class_diario:
    def __init__(self, fecha, out_carpeta, filtrado=1):
        self.fecha0 = datetime.strptime(fecha, '%Y%m%d')
        self.out_carpeta = out_carpeta
        os.makedirs(self.out_carpeta, exist_ok=True)
        self.get_data(filtrado)

    def get_data(self, filtrado):
        from lst_class_horario import lst_class_horario
        # Extraemos los datos entre 0 y 12 UTC (21-9 Hora Local)
        horas = np.arange(0,13)
        arrays = []
        for hora in horas:
            fecha = self.fecha0.strftime('%Y%m%d') + str(hora).zfill(2) + '00'
            a = lst_class_horario(fecha, self.out_carpeta)
            if hora == 0:
                self.lats = a.lats
                self.lons = a.lons
            if filtrado == 1:
                arrays.append(a.LST_filtrado)
            else:
                arrays.append(a.LST)
        F = ma.stack(arrays, axis=0)
        c_val = ma.count(F, axis=0)  # Contamos la cantidad de datos
        l_val = ma.min(F, axis=0)  # Calculamos la Tmin
        #new_mask = (c_val < 7) | ma.getmask(l_val)
        #new_mask = (c_val < 7)
        new_mask = ma.getmask(l_val)
        self.LSTmin = ma.masked_array(l_val.data, new_mask)
        self.LSTcount = c_val
    
    def save_gtiff(self, file_name):
        from osgeo import gdal, osr, ogr
        fillvalue = -9999.9
        extent = [-76.5, -52, -56.5, -20.]
        #extent = [-93.0, -60.00, -25.00, 18.00]
        # Get GDAL driver GeoTiff
        driver = gdal.GetDriverByName('GTiff')
        # Get dimensions
        data = ma.filled(self.LSTmin, fill_value=fillvalue)
        nlines = data.shape[0]
        ncols = data.shape[1]
        data2 = ma.filled(self.LSTcount, fill_value=fillvalue)
        data_type = gdal.GDT_Float32
        ########################
        # Create a temp grid
        grid_data = driver.Create('grid_data', ncols, nlines, 1, data_type)#, options)
        # Write data for LSTmin y LSTcount
        grid_data.GetRasterBand(1).SetNoDataValue(fillvalue)
        grid_data.GetRasterBand(1).WriteArray(data)
        grid_data.GetRasterBand(2).SetNoDataValue(fillvalue)
        grid_data.GetRasterBand(2).WriteArray(data2)
        # Lat/Lon WSG84 Spatial Reference System
        srs = osr.SpatialReference()
        srs.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
        # Setup projection and geo-transform
        grid_data.SetProjection(srs.ExportToWkt())
        grid_data.SetGeoTransform(getGeoTransform(extent, nlines, ncols)) 
        
        grid_data.FlushCache()
        # Save the file
        print(f'Generated GeoTIFF: {file_name}')
        driver.CreateCopy(self.out_carpeta + file_name, grid_data, 0)
        
        # Close the file
        driver = None
        grid_data = None



    def mapa_helada_ORA(self):
        from figure_functions import mapa_temperatura_minima_superficial
        l_lat = [-56.5, -20. ]
        l_lon = [-76.5, -52. ]
        nfig = self.out_carpeta + 'LSTmin_ORA_' + self.fecha0.strftime('%Y%m%d') + '.jpg'
        mapa_temperatura_minima_superficial(self.lons, self.lats, self.LSTmin, l_lat, l_lon, nfig)
    
    def mapa_helada_diario(self):
        import matplotlib.colors as mcolors
        import matplotlib.cm as cm
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        import cartopy
        import cartopy.crs as ccrs

        from matplotlib.axes import Axes
        from cartopy.mpl.geoaxes import GeoAxes
        GeoAxes._pcolormesh_patched = Axes.pcolormesh

        from cartopy.io.shapereader import Reader
        from cartopy.feature import ShapelyFeature

        
        figname = self.out_carpeta + 'LSTmin_' + self.fecha0.strftime('%Y%m%d') + '.jpg'
        
        provincias = cartopy.feature.NaturalEarthFeature(category='cultural',
                                                 name='admin_1_states_provinces_lines',
                                                 scale='10m',
                                                 facecolor='none')
        paises = cartopy.feature.NaturalEarthFeature(category='cultural',
                                             name='admin_0_countries',
                                             scale='10m',
                                             facecolor='none')
        lst_colors = [ [255/255,255/255,255/255], # +12 to +10°C
               [ 74/255,  8/255,143/255], # +10 to  +8°C
               [ 51/255,  4/255,180/255], #  +8 to  +6°C
               [ 45/255, 46/255,255/255], #  +6 to  +4°C
               [ 38/255,105/255,247/255], #  +4 to  +2°C
               [  2/255,128/255,250/255], #  +2 to   0°C
               [ 82/255,176/255,250/255], #   0 to  -2°C
               [132/255,215/255,249/255], #  -2 to  -4°C
               [160/255,227/255,196/255], #  -4 to  -6°C
               [182/255,235/255,199/255], #  -6 to  -8°C
               [207/255,236/255,204/255], #  -8 to -10°C
               [221/255,242/255,209/255], # -10 to -12°C
               [233/255,242/255,209/255], # -12 to -14°C
               [218/255,220/255,217/255]  # -14 to -16°C
             ]
        cmap_lst = mcolors.ListedColormap(lst_colors)
        cmap_lst.set_bad(color='green')
        bounds = np.array([-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12,14,16])
        norm_topo = mcolors.BoundaryNorm(boundaries=bounds, ncolors=14)

        # Definimos proyeccion
        projection = ccrs.PlateCarree()
        extent = [-76.5, -52, -56.5, -20.]

        fig=plt.figure(figsize=[11,10])

        ax = fig.add_subplot(111, projection=projection)
        ax.set_extent(extents=extent, crs=projection)

        img = plt.pcolormesh(self.lons, self.lats, self.LSTmin, vmin=-16, vmax=12, cmap=cmap_lst,
                           transform=ccrs.PlateCarree())
        plt.colorbar(img, pad=0.01, aspect=42, shrink=0.5, ticks=np.arange(-16,12.1,2))

        ax.add_feature(provincias, linewidth=0.25, edgecolor='k', facecolor='None')
        ax.add_feature(paises, linewidth=1.0, edgecolor='k', facecolor='None')

        # Agregamos lineas y marcas lat/lon
        gl = ax.gridlines(crs=ccrs.PlateCarree(), color='gray', alpha=1.0, linestyle='--', linewidth=0.1,
                          xlocs=np.arange(-180, 180, 5), ylocs=np.arange(-90, 90, 5), draw_labels=True)
        gl.top_labels = False
        gl.right_labels = False
        
        plt.title('GOES-16 LST minima', fontsize=12, loc='left')
        plt.title(self.fecha0.strftime('%Y-%m-%d 0-12 UTC'), fontsize=12, loc='right')
        plt.savefig(figname, dpi=150, bbox_inches='tight')
