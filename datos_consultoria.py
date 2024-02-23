import time
import numpy as np
import numpy.ma as ma
import pandas as pd
import datetime as dt

import sys
sys.path.append('./lib')
from lst_class_horario_consultoria import lst_horario

start = time.time()

mes = 10
nmes = 'Octubre'

if mes in [10,11,12]:
    years = np.arange(2017,2024)
elif mes in [1,2,3,4]:
    years = np.arange(2019,2025)

for year in years:
    print('Trabajando el año', year)
    ini = str(year)+'1001'
    fin = str(year)+'1101'
    fechas = pd.date_range(start=ini, end=fin, freq='H')
    print(fechas[0:-1])
    ntime = len(fechas[0:-1])
    ny = 290
    nx = 245
    LST = np.empty((ntime, ny, nx))
    LST[:] = np.nan
    DQF = np.empty((ntime, ny, nx))
    DQF[:] = np.nan
    for i,fecha in enumerate(fechas[0:-1]):
        fechad = fecha.strftime('%Y%m%d%H%M')
        print('Extrayendo datos de', fechad)
        b = lst_horario(fechad, './salidas/' + fecha.strftime('%Y%m%d') + '/')
        LST[i,:,:] = b.LST
        DQF[i,:,:] = b.DQF
        '''
        if (i == 0) & (year == years[0]):
            lat = b.lats
            lon = b.lons
            nfile2 = './salidas/lats.npy'
            nfile3 = './salidas/lons.npy'
            np.save(nfile2, lat)
            np.save(nfile3, lon)
        '''
    nfile0 = './salidas/output_2017/LST_' + nmes + '_' + str(year) + '.npy'
    nfile1 = './salidas/output_2017/DQF_' + nmes + '_' + str(year) + '.npy'
    np.save(nfile0, LST)
    np.save(nfile1, DQF)
    



end = time.time()

secs = np.round(end - start, 2)
print('Se demoro ', secs, ' segundos.')