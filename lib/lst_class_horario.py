import netCDF4
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

# FECHA Y HORA DE INTERES PARA ANALIZAR
# NOTA: este producto se disponibiliza uno por hora
#fecha Formato: yyyymmddHHMM (La hora en UTC)

class lst_class_horario:
    def __init__(self, fecha, out_carpeta):
        self.fecha0 = datetime.strptime(fecha, '%Y%m%d%H%M')
        self.prefijo = 'noaa-goes16/ABI-L2-LSTF/'
        self.carpeta = self.prefijo + self.fecha0.strftime('%Y') +\
              '/' + self.fecha0.strftime('%j') + '/' + self.fecha0.strftime('%H')+'/'
        self.out_carpeta = out_carpeta
        os.makedirs(self.out_carpeta, exist_ok=True)
        self.get_data()

    def get_data(self):
        import s3fs
        # ROI (Argentina)
        x1_lon = 520
        x2_lon = 765
        y1_lat = 735 
        y2_lat = 1025
        #
        fs = s3fs.S3FileSystem(anon=True)
        files = fs.ls(self.carpeta)
        # Nos quedamos con el archivo que contenga el M6: FullDisk
        self.archivo = [s for s in files if 'OR_ABI-L2-LSTF-M6' in s][0]
        if not self.archivo:
            print('No existe el archivo para la fecha indicada')
            print(files)
            print(self.carpeta)
            exit()
        fname = self.archivo
        fname0 = fname.split('/')[-1]
        with fs.open(fname) as f:
            with netCDF4.Dataset(fname0, memory=f.read()) as goes:
                H = goes.variables['goes_imager_projection'].getncattr('perspective_point_height')
                lon_0 = goes.variables['goes_imager_projection'].getncattr('longitude_of_projection_origin')
                sat_sweep = goes.variables['goes_imager_projection'].getncattr('sweep_angle_axis')
                x = goes.variables['x'][x1_lon:x2_lon] * H
                y = goes.variables['y'][y1_lat:y2_lat] * H
                xv, yv = np.meshgrid(np.array(x), np.array(y))
                # Doc: https://proj.org/operations/projections/geos.html
                geo = Proj(proj='geos', h=H, lon_0=lon_0, sweep=sat_sweep)
                lons_lst, lats_lst = geo(xv, yv, inverse=True)
                add_seconds = int(goes.variables['time_bounds'][0])
                date = datetime(2000,1,1,12) + timedelta(seconds=add_seconds)

                LST = goes.variables['LST'][y1_lat:y2_lat,x1_lon:x2_lon][::1,::1]
                DQF = goes.variables['DQF'][y1_lat:y2_lat,x1_lon:x2_lon][::1,::1]
                PQI = goes.variables['PQI'][y1_lat:y2_lat,x1_lon:x2_lon][::1,::1]

                # Guardamos variables a la clase
                self.H = H
                self.lon_0 = lon_0
                self.sat_sweep = sat_sweep
                self.lats = lats_lst
                self.lons = lons_lst
                self.LST = LST - 273.15
                self.DQF = DQF
                self.PQI = PQI
                mask_lst = np.copy(DQF)
                # Asigno 0 donde DQF==0 (calidad buena de la estimacion LST)
                mask_lst[DQF < 0.5] = 0
                # Asigno 1 donde la calidad de la estimacion de LST es media o baja.
                mask_lst[DQF >= 0.5] = 1
                # Generamos la nueva mascara incluyendo lo que ya venia enmascarado.
                self.mascara = mask_lst | ma.getmask(LST)
                self.LST_filtrado = ma.masked_array(LST.data, self.mascara)
    
    def save_gtiff(self, file_name):
        from osgeo import gdal, osr, ogr
        fillvalue = -9999.9
        extent = [-76.5, -52, -56.5, -20.]
        #extent = [-93.0, -60.00, -25.00, 18.00]
        # Get GDAL driver GeoTiff
        driver = gdal.GetDriverByName('GTiff')
        # Get dimensions
        data = ma.filled(self.LST_filtrado, fill_value=fillvalue)
        nlines = data.shape[0]
        ncols = data.shape[1]
        #nbands = len(data.shape)
        data_type = gdal.GDT_Float32
        # Create a temp grid
        #options = ['COMPRESS=JPEG', 'JPEG_QUALITY=80', 'TILED=YES']
        grid_data = driver.Create('grid_data', ncols, nlines, 1, data_type)#, options)
        # Write data for each bands
        grid_data.GetRasterBand(1).SetNoDataValue(fillvalue)
        grid_data.GetRasterBand(1).WriteArray(data)
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

    def save_map_lst(self, opt=1):
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

        if opt == 1:
            LST = self.LST_filtrado
            figname = self.out_carpeta + 'LST_filtrada_' + self.fecha0.strftime('%Y%m%d%H') + '.jpg'
        else:
            LST = self.LST
            figname = self.out_carpeta + 'LST_' + self.fecha0.strftime('%Y%m%d%H') + '.jpg'
        
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

        img = plt.pcolormesh(self.lons,self.lats,LST, vmin=-16, vmax=12, cmap=cmap_lst,
                           transform=ccrs.PlateCarree())
        plt.colorbar(img, pad=0.01, aspect=42, shrink=0.5, ticks=np.arange(-16,12.1,2))

        ax.add_feature(provincias, linewidth=0.25, edgecolor='k', facecolor='None')
        ax.add_feature(paises, linewidth=1.0, edgecolor='k', facecolor='None')

        # Agregamos lineas y marcas lat/lon
        gl = ax.gridlines(crs=ccrs.PlateCarree(), color='gray', alpha=1.0, linestyle='--', linewidth=0.1,
                          xlocs=np.arange(-180, 180, 5), ylocs=np.arange(-90, 90, 5), draw_labels=True)
        gl.top_labels = False
        gl.right_labels = False
        
        plt.title('GOES-16 LST', fontsize=12, loc='left')
        plt.title(self.fecha0.strftime('%Y-%m-%d %H:%M UTC'), fontsize=12, loc='right')
        plt.savefig(figname, dpi=150, bbox_inches='tight')

    def save_map_DQF(self):
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
        
        provincias = cartopy.feature.NaturalEarthFeature(category='cultural',
                                                 name='admin_1_states_provinces_lines',
                                                 scale='10m',
                                                 facecolor='none')
        paises = cartopy.feature.NaturalEarthFeature(category='cultural',
                                             name='admin_0_countries',
                                             scale='10m',
                                             facecolor='none')
        dqf_colors = [
               [ 51/255,204/255,  0/255], # Alta
               [255/255,153/255,  0/255], # Media
               [255/255,  0/255,  0/255], # Baja
               [218/255,220/255,217/255]  # No Data
             ]
        cmap_dqf = mcolors.ListedColormap(dqf_colors)
        # Definimos proyeccion
        projection = ccrs.PlateCarree()
        extent = [-76.5, -52, -56.5, -20.]
        fig=plt.figure(figsize=[11,10])
        ax = fig.add_subplot(111, projection=projection)
        ax.set_extent(extents=extent, crs=projection)

        img=plt.pcolormesh(self.lons,self.lats,self.DQF,vmin=0,vmax=4,cmap=cmap_dqf,
                           transform=ccrs.PlateCarree())
        cbar = plt.colorbar(img, pad=0.01, aspect=18, shrink=0.5)
        cbar.ax.get_yaxis().set_ticks([])
        for j, lab in enumerate(['$ALTA$','$MEDIA$','$BAJA$','$SinDato$']):
            cbar.ax.text(.5, j + 0.6, lab, ha='center', va='center', color='k', rotation=270)
        ax.add_feature(provincias, linewidth=0.25, edgecolor='k', facecolor='None')
        ax.add_feature(paises, linewidth=1.0, edgecolor='k', facecolor='None')
        # Agregamos lineas y marcas lat/lon
        gl = ax.gridlines(crs=ccrs.PlateCarree(), color='gray', alpha=1.0, linestyle='--',
                          linewidth=0.1, xlocs=np.arange(-180, 180, 5), ylocs=np.arange(-90, 90, 5), draw_labels=True)
        gl.top_labels = False
        gl.right_labels = False

        plt.title('GOES-16 LST Calidad [DQF]', fontsize=12, loc='left')
        plt.title(self.fecha0.strftime('%Y-%m-%d %H:%M UTC'), fontsize=12, loc='right')
        figname = self.out_carpeta + 'DQF_' + self.fecha0.strftime('%Y%m%d%H') + '.jpg'
        plt.savefig(figname, dpi=150, bbox_inches='tight')

        