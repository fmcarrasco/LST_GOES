import glob
import os
import pandas as pd
from osgeo import gdal
from affine import Affine
import numpy as np
import datetime as dt
import netCDF4 as nc
import xarray as xr
from pyproj import Transformer, transform

'''
La carpeta de datos se va a guardar en esta computadora en:

D:/LST_GOES/TIF/${year}/${month}

y los archivos tienen el siguiente formato:

nombre: TIF_LST_goes_min_${yyyymmdd}.tif
'''

class lst_class:
    def __init__(self, fechai, fechaf='SinFecha'):
        self.prefijo = 'TIF_LST_goes_min_'
        self.fecha_i = dt.datetime.strptime(fechai, '%Y%m%d')
        self.carpeta = 'D:/LST_GOES/TIF/'
        if fechaf != 'SinFecha':
            self.fecha_f = dt.datetime.strptime(fechaf, '%Y%m%d')
        else:
            self.fecha_f = 'SinFecha'
        self.get_data()

    def get_data(self):
        if self.fecha_f != 'SinFecha':
            fechas = pd.date_range(start=self.fecha_i, end=self.fecha_f)
            archivos = []
            for fecha in fechas:
                year = str(fecha.year)
                mes = str(fecha.month).zfill(2)
                str_fecha = fecha.strftime('%Y%m%d')
                fname = self.carpeta + year + '/' + mes + '/' + self.prefijo + str_fecha +'.tif'
                if os.path.exists(fname):
                    archivos.append(fname)
                else:
                    print('No hay dato para LST-GOES en:', fname)
                    continue
        else:
            d1 = self.fecha_i.strftime('%Y%m%d')
            year = str(self.fecha_i.year)
            mes = str(self.fecha_i.month).zfill(2)
            fname = self.carpeta + year + '/' + mes + '/' + self.prefijo + d1 + '.tif'
            if os.path.exists(fname):
                archivos = [fname]
            else:
                print('No hay dato para LST-GOES en:', fname)
                exit()
            
        FirstTime = True
        for it, archivo in enumerate(archivos):
            print(archivo)
            ds = gdal.Open(archivo)
            data = ds.ReadAsArray()
            data[data>1000.] = np.nan
            data[data<-1000.] = np.nan            
            f1 = archivo.split('/')[-1]
            f2 = f1.split('.')[0].split('_')[-1]
            ffecha = dt.datetime.strptime(f2, '%Y%m%d')
            if FirstTime:
                lst_min = np.zeros((len(archivos), data.shape[0], data.shape[1]))
                indices = np.indices(data.shape)
                T0 = Affine.from_gdal(*ds.GetGeoTransform())
                conv = T0*Affine.translation(0.5, 0.5)  # Relativo al centro
                lon, lat = conv*(indices[1,:,:], indices[0,:,:])  # Row/Column
                tiempos = []
                FirstTime = False
            lst_min[it,:,:] = data
            tiempos.append(ffecha)
        # End of LOOP
        if not archivos:
            self.SinDatos = True
        else:
            self.SinDatos = False
        self.lst_min = np.squeeze(lst_min)
        self.lat = lat
        self.lon = lon
        self.tiempos = tiempos

    def get_nearest_index(self, loni, lati):
        # Transformar de EPSG:4326 a ECEF (EPSG:4978)
        # GtiFF viene en EPSG:4326 y pasamos a coord cartesianas EPSG:4978
        transformer = Transformer.from_crs(4326, 4978)
        X, Y = transformer.transform(self.lon, self.lat)
        x, y = transformer.transform(loni,lati)
        distance = ((x-X)**2 + (y-Y)**2)**0.5
        pos_index = np.unravel_index(distance.argmin(), distance.shape)

        return pos_index

    def get_xarray_from_data(self):
        if len(self.lst_min.shape) > 2:
            lstmin = self.lst_min
        else:
            lstmin = np.expand_dims(self.lst_min, axis=0)
        time_unit = 'days since 2000-01-01 00:00:00'
        calendar_u = 'proleptic_gregorian'
        tiempos_nc = nc.date2num(self.tiempos, units=time_unit, calendar=calendar_u)
        da = xr.DataArray(data=lstmin, dims=['time','y', 'x'],
                          coords={'time': tiempos_nc,
                                  'lat': (['y', 'x'], self.lat),\
                                  'lon': (['y', 'x'], self.lon)})
        
        d1 = 'Temperatura superficial minima del suelo estimada por GOES entre 21 y 9 AM por SMN'
        da.attrs = dict(description=d1, units='mm/day')
        da.lon.attrs = {"unit":"degrees_east","standard_name": "longitude", "long_name" : "longitude"}
        da.lat.attrs = {"unit":"degrees_north","standard_name": "latitude", "long_name" : "latitude"}
        da.time.attrs['units'] = time_unit
        return da
