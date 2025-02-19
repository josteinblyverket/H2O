import sys, os, glob
import os.path
from sys import argv
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
from shapely.geometry import Polygon, MultiPolygon, shape, Point
import geopandas as gp
import xesmf as xe
import pyproj
import cartopy.crs as ccrs
from datetime import datetime, timedelta

#usage python Malestasjon_production_runoff.py 202411201100 /lustre/storeA

#read date from argument format YYYYMMDDHH00
#using exising weights runs about 1 min

#In [1]: %time run Malestasjon_production_runoff.py 202411210300 /lustre/storeA
#CPU times: user 55.5 s, sys: 8.57 s, total: 1min 4s
#Wall time: 1min 5s
def datespan(startDate, endDate, delta=timedelta(days=1)):
     currentDate = startDate
     while currentDate < endDate:
         yield currentDate
         currentDate += delta

def listOfDates(dtg_start, dtg_end, dtg_step):
    
    year_start = int(dtg_start[0:4])
    year_end = int(dtg_end[0:4])
    mm_start = int(dtg_start[4:6].strip("0"))
    mm_end = int(dtg_end[4:6].strip("0"))

    day_start = int(dtg_start[6:8].strip("0"))
    day_end = int(dtg_end[6:8].strip("0"))

    day_start = 12
    day_end = 1
    
    hour = int(dtg_start[8:10])
    
    time = []
    print(day_start)
    print(day_end)
    
    for timestamp in datespan(datetime(year_start, mm_start, day_start, hour, 00),
                               datetime(year_end, mm_end, day_end, hour, 00),
                               delta=timedelta(hours=dtg_step)):
          
          time.append(str(timestamp))
                    
    return time

def convert_3D_2D(geometry):
    '''
    Takes a GeoSeries of 3D Multi/Polygons (has_z) and returns a list of 2D Multi/Polygons
    '''
    new_geo = []
    for p in geometry:
        #print(p)
        if p.has_z:
            if p.geom_type == 'Polygon':
                #print('polygon')
                lines = [xy[:2] for xy in list(p.exterior.coords)]
                new_p = Polygon(lines)
                new_geo.append(new_p)
            elif p.geom_type == 'MultiPolygon':
                new_multi_p = []
                for ap in p.geoms:
                    lines = [xy[:2] for xy in list(ap.exterior.coords)]
                    new_p = Polygon(lines)
                    new_multi_p.append(new_p)
                new_geo.append(MultiPolygon(new_multi_p))
    return new_geo


dtg_start = "2018030100"
dtg_end = "2018080100"
dtg_step = 24

time = listOfDates(dtg_start, dtg_end, dtg_step)    


#date="201803200600" #argv[1]
#lustrep=argv[2] #/lustre/storeA
#sys.exit()

#to make the weights eg
#qlogin -l h_rss=32G,mem_free=32G,h_data=32G -pe shmem-1 1 -q bigmem-r8.q 

#redhat 
#source /modules/rhel8/user-apps/suv-modules/miniconda3/24.7.1/etc/profile.d/conda.sh
#conda activate /home/helenebe/.conda/env/squeegee

'''
 env recipie:
 2069  conda create -y -p .conda/env/squeegee
 2070  conda activate /home/helenebe/.conda/env/squeegee
 2071  conda install -c conda-forge xesmf
 2072  conda install -c conda-forge cartopy matplotlib dask netCDF4
conda deactivate #(ESMpy 8.4 bug https://xesmf.readthedocs.io/en/stable/installation.html)
conda activate /home/helenebe/.conda/env/squeegee
'''

#starting offline runs to feed into future SURFEX DIFF ES runs and products
#one product is estimates of river runoff DRAINC and RUNOFFC
#gathered 1st by squeegee xesmf

#here I test using squeeege on full domain
# Do Norway first

#first data with the variable is 20th nov 2024 stamped 10
#/2024/11/20/09/SURFOUT.20241120_10h00.nc
#fix SURFOUT.
pathf= "/ec/res4/scratch/nor3005/sfx_data/CARRA_Land_Pv1_v2/archive/"

for i in range(len(time)):

    dated=pd.to_datetime(time[i])         #pd.to_datetime(date) #('202411200900')
    print(dated)
    datef=dated+pd.Timedelta('24h')
    filesp=pathf+dated.strftime('%Y/%m/%d/06/')+datef.strftime('SURFOUT.%Y%m%d_06h00.nc')

    print(filesp)
    geoin=xr.open_mfdataset(filesp,cache=False)

    Rarome=6371229 
    #grid on sphere (unlike real globe)

    proj_string = "+proj=lcc +lat_0="+str(geoin.LAT0.data.item())+" +lon_0="+str(geoin.LON0.data.item())+" +lat_1="+str(geoin.LAT0.data.item())+" +lat_2="+str(geoin.LAT0.data.item())+" +units=m +no_defs +R=" + str(Rarome)

    XX=np.arange(geoin.DX[0,0],geoin.DX[0,0]*(geoin['xx'].shape[0]+1), geoin.DX[0,0])
    YY=np.arange(geoin.DY[0,0],geoin.DY[0,0]*(geoin['yy'].shape[0]+1), geoin.DY[0,0])

    myP=pyproj.Proj(proj_string)

    false_easting, false_northing = myP(geoin.LONORI.data, geoin.LATORI.data,inverse=False)

    #works with LLC grid at least
    tMNx,tMNy = XX-XX[0]/2.+false_easting,YY-YY[0]/2.+false_northing

    geoin = geoin.rename({'xx': 'x', 'yy': 'y'})
    geoin=geoin.assign_coords(x=tMNx)
    geoin=geoin.assign_coords(y=tMNy)
    lons, lats = myP(*np.meshgrid(tMNx,tMNy),inverse=True)


    geoin.coords['lat'] = (('y','x'),lats)
    geoin.coords['lon'] = (('y','x'),lons)
    geoin.lon.attrs={'units':'degrees_east'}#)
    geoin.lat.attrs={'units':'degrees_north'}
    geoin.set_coords(['lat','lon'])

    geoin.attrs['pyproj_srs']=myP.srs

    #need lon_b lat_b for conservative regridding
    dx, dy = geoin.DX.data[0,0], geoin.DY.data[0,0]
    Xcorners=np.arange(geoin['x'].data[0]-dx/2., geoin['x'].data[-1]+3*dx/2., dx)
    Ycorners=np.arange(geoin['y'].data[0]-dy/2., geoin['y'].data[-1]+3*dy/2., dy)

    Lon2b, Lat2b = myP(*np.meshgrid(Xcorners,Ycorners),inverse=True) #
    geoin.coords['xb'] = (Xcorners)
    geoin.coords['yb'] = (Ycorners)
    geoin.coords['lat_b'] = (('yb','xb'),Lat2b)
    geoin.coords['lon_b'] = (('yb','xb'),Lon2b)
    geoin.set_coords(['lat_b','lon_b'])
    geoin['mask']=geoin.RUNOFFC_ISBA.isnull()
    geoin['mask'].values=np.where(~geoin.RUNOFFC_ISBA.isnull(),1,0)

    #read shape files 
    #Load some polygons 
    #ccm2_basins = gp.read_file('shapefiles/ccm21/WGS84_W2008.gdb',layer='SEAOUTLETS')
    #regs2 = gp.read_file('shapefiles/Nedborfelt/Nedborfelt_Vassdragsomr.shp')
    regs2 = gp.read_file('/ec/res4/scratch/sbjb/Projects/CERISE/Hydrology/Hydrologi_TotalNedborfeltMalestasjon.shp')

    print(regs2.columns)

    #subset to MN domain
    print(geoin.lon.min().values)
    print(geoin.lon.max().values)
    print(geoin.lat.min().values)
    print(geoin.lat.max().values)
    
    
#    regs2sub = regs2.cx[geoin.lon.min().values:geoin.lon.max().values,geoin.lat.min().values:geoin.lat.max().values]
    regs2sub = regs2.cx[10.0:32.0,geoin.lat.min().values:72.0]
    

    #Is needed bc shape file format has Z dim
    print(regs2sub.columns)
    regs2sub.loc[:,'geometry'] = convert_3D_2D(regs2sub.geometry)

    #print(regs2sub)
    print(regs2sub.columns)

    #remove any catchments / basins below 1000 km2
    resultn=regs2sub[regs2sub.areal_km2 > 50]

    #Do the Spatial averaging, use exsisting file with weights.
    #If first time access more memory and set reuse_weights to False.
    #This will take a very long time, many hours for this case..
    savg = xe.SpatialAverager(geoin, resultn.geometry, geom_dim_name="stID",filename='/ec/res4/scratch/sbjb/Projects/CERISE/Hydrology/CARRA_NE_sub_1000km2_Malestasjon_xesmf_weights.nc',reuse_weights=True)


    #Add runoff and discharge accumulated last hour
    tmp1=geoin['RUNOFFC_ISBA'] + geoin['DRAINC_ISBA']
    tmp1=tmp1.to_dataset(name='mrunoff')

    #Make sure the data has a mask
    tmp1['mask']=geoin['mask']

    out=savg(tmp1)
    out = out.assign_coords(stID=xr.DataArray(resultn.stID.values, dims=("stID",)))

    #write output to file.
    #The output is further processed in Merge_csvs_Malestasjon.py
    outpd=out.to_pandas()
    outpd.to_csv('/ec/res4/scratch/sbjb/Projects/CERISE/Hydrology/RiverDischarge/CTL/SpatialAverager_ldas_runoff_Malestasjon_'+datef.strftime('%Y%m%d%H00')+'.csv')
    print('written SpatialAverager_ldas_runoff_Malestasjon_'+datef.strftime('%Y%m%d%H00')+'.csv')

    writeRunNC=False
    if writeRunNC:
        #can add coastlines in ncview with the meta on the file
        tmp1.attrs['pyproj_srs']=myP.srs
        tmp1.to_netcdf('/ec/res4/scratch/sbjb/Projects/CERISE/Hydrology/test.nc')

        #close input file to be tidy
    geoin.close()


