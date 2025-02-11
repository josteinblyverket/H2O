import sys, os, glob
import os
import numpy as np
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
#matplotlib.use('Agg')
import pandas as pd
from shapely.geometry import Polygon, MultiPolygon, shape, Point
import geopandas as gp
import xesmf as xe
import pyproj
import requests
import cartopy.crs as ccrs
from cartopy.io.img_tiles import Stamen
import cartopy.feature as cfeature
tiler = Stamen('terrain-background')


def convert_3D_2D(geometry):
    '''
    Takes a GeoSeries of 3D Multi/Polygons (has_z) and returns a list of 2D Multi/Polygons
    '''
    new_geo = []
    for p in geometry:
        if p.has_z:
            if p.geom_type == 'Polygon':
                lines = [xy[:2] for xy in list(p.exterior.coords)]
                new_p = Polygon(lines)
                new_geo.append(new_p)
            elif p.geom_type == 'MultiPolygon':
                new_multi_p = []
                for ap in p:
                    lines = [xy[:2] for xy in list(ap.exterior.coords)]
                    new_p = Polygon(lines)
                    new_multi_p.append(new_p)
                new_geo.append(MultiPolygon(new_multi_p))
    return new_geo


def getDomainData(expdirs):

    # Assuming same forcing
    #FORCING=xr.open_dataset(expdirs+'forcing/2023050100/FORCING.nc', cache=False,use_cftime=False)
    FORCING=xr.open_dataset('/lustre/storeB/users/josteinbl/sfx_data/LDAS_NOR/archive/2022/12/12/06/raw.nc', cache=False,use_cftime=False)
    print(FORCING.latitude.min(), FORCING.latitude.max())
    print(FORCING.longitude.min(), FORCING.longitude.max())

    fPREP=xr.open_dataset('/lustre/storeB/users/josteinbl/sfx_data/LDAS_NOR/'+'archive/2023/04/01/00/SURFOUT.nc', cache=False,use_cftime=False)
    fPGD=xr.open_dataset('/lustre/storeB/users/josteinbl/sfx_data/LDAS_NOR/'+'climate/PGD.nc', cache=False,use_cftime=False)

    # Add coords etc so that surfex files may be processed
    R_pysurfex=6371000 # and met nordic
    R_arome   =6371229
    Rpy=6.371229e+06
    proj_string = "+proj=lcc +lat_0=63 +lon_0=15 +lat_1=63 +lat_2=63 +units=m +no_defs +R=" + str(Rpy)
    myP=pyproj.Proj(proj_string)

    lcc = ccrs.LambertConformal(#globe=globe, # no datumshift in GCMs
        central_longitude=15, central_latitude=63,
        standard_parallels=(63,63))#,
    geodetic=ccrs.Geodetic() #default WGS84

    return myP, proj_string, fPREP, fPGD, FORCING, lcc


def readStationList(stationfile):

    #read a list of stations
    stationst=stationfile
    stations=pd.read_csv(stationst, header=0,index_col=0,sep=None,
                         dtype={'GauID': 'str'})
    stations=stations.drop_duplicates(subset=('GauID'), keep='last')

    # Nr.	GauID	Lat	Long	Ara_km2	Group
    daily = 1440
    station='2.25.0'# stations.GauID[2]

    stations['GauID']=stations.GauID+'.0'

    return stations

    
def readShapeFiles(shapefiles, stations):

    # or Outlets to sea defined for Europe CCM2
#    ccm2_basins = gp.read_file(shapefiles,layer='SEAOUTLETS')
    #ccm2sub = ccm2_basins.cx[FORCING.LON.min().values:FORCING.LON.max().values,FORCING.LAT.min().values:FORCING.LAT.max().values]    

    # or cathements where NVE have measuring stations
    #regs = gp.read_file('/lustre/storeB/project/nwp/H2O/wp4/RRdata/shapefiles/utm33shp/NVEData/Hydrologi/Hydrologi_TotalNedborfeltMalestasjon.shp')
    regs2 = gp.read_file(shapefiles)

    print("prior")
    print(regs2.columns)
    regs2=regs2.loc[regs2['stID'].isin(stations.GauID)]

    print("posterior")
    #inspect the shape file
    print(regs2.crs)
    print(regs2.columns)

    #remove extra dim in nve shape-files
    regs2.geometry = convert_3D_2D(regs2.geometry)

    plot_polygons=False
    #make a plot of the polygons:
    if plot_polygons==True:
    
        fig = plt.figure(figsize=(4,8))
        ax1 = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(), frameon=False)
 
        ax1.coastlines(resolution='10m')
        toplot= regs2.cx[FORCING.longitude.min().values:FORCING.longitude.max().values,FORCING.latitude.min().values:FORCING.latitude.max().values]

        toplot.plot(column='stID',kind='geo',edgecolor="g",linewidth=0.2,
                   ax=ax1,legend=True,cmap='Reds',transform=ccrs.Geodetic())
        #,legend_kwds={"orientation": "horizontal", "pad": 0.01, "label": "VassOmrNr"})
        #ccm2sub.plot(column='AREA_KM2',kind='geo',edgecolor="g",linewidth=0.2,ax=ax1,legend=True,cmap='Reds',legend_kwds={"orientation": "horizontal", "pad": 0.01})
        ax1.axis('off')
        plt.title("")
        plt.savefig("polygons.png")
        plt.show()

    return regs2

def readModelData(expdirs, expdirfs, fPREP, fPGD, ncstart, expnam, variables):

    numruns=len(expnam)

    #-----------------------------------------
    # read data in this way because a large amount of data
    def preprocess(ds):
        return ds[variables]

    dictds={}

    for i in np.arange(0,numruns):
        dictds[expnam[i]] = xr.open_mfdataset(expdirs+expdirfs[i]+ncstart[i]+
                                          'sel2023??.nc',
                                          cache=False,
                                          preprocess=preprocess,
                                          concat_dim="time",
                                          combine="nested",
                                          chunks={'time': 1})
        print(expnam[i] + ' given to: ' +expdirs+expdirfs[i]+ncstart[i]+
              'sel2023??.nc')

    print(dictds)

    return dictds

def generateWeights(dictds, fPREP, fPGD, FORCING, regs2, weightfile):

    # xesmf needs dataset to generate weights
    TSr=dictds[expnam[-1]]['RUNOFFC_ISBA'].isel(time=0).to_dataset()
    print('TSr')
    print(TSr)

    xfalseaR, yfalseaR = myP(fPREP.LONORI.data, fPREP.LATORI.data,inverse=False)
    #Surfex makes all coords positive by setting the origin to the lower left corner

    TSr=TSr.assign_coords(x= fPREP.XX[0,:].data+xfalseaR)
    TSr=TSr.assign_coords(y= fPREP.YY[:,0].data+yfalseaR)
    TSr = TSr.rename({'xx': 'x', 'yy': 'y'})

    dx, dy = fPGD.DX.data[0,0], fPGD.DY.data[0,0]
    Xcorners=np.arange(TSr['x'].data[0]-dx/2., TSr['x'].data[-1]+3*dx/2., dx)
    Ycorners=np.arange(TSr['y'].data[0]-dy/2., TSr['y'].data[-1]+3*dy/2., dy)
    #Lon2, Lat2 = myP(fPREP.XX.data+xfalseaR,fPREP.YY.data+yfalseaR,inverse=True)
    Lon2b, Lat2b = myP(*np.meshgrid(Xcorners,Ycorners),inverse=True) 

    TSr.coords['xb'] = (Xcorners)
    TSr.coords['yb'] = (Ycorners)
    TSr.coords['lat_b'] = (('yb','xb'),Lat2b)
    TSr.coords['lon_b'] = (('yb','xb'),Lon2b)
    TSr.set_coords(['lat_b','lon_b'])

    print("Forcing lat data")
    print(np.shape(FORCING.latitude.data))
    #print(TSr.coords['lat_b'])

    TSr.coords['lat'] = (('y','x'),FORCING.latitude.data)
    TSr.coords['lon'] = (('y','x'),FORCING.longitude.data)
    TSr.lon.attrs=FORCING.longitude.attrs
    TSr.lat.attrs=FORCING.longitude.attrs
    TSr.set_coords(['lat','lon'])

    TSr['mask']=TSr.RUNOFFC_ISBA.isnull()
    TSr['mask'].values=np.where(~TSr.RUNOFFC_ISBA.isnull(),1,0)

    TSr.attrs['pyproj_srs']=proj_string

    #make an outer domain polygon to crop the input shape file

    dompoly=Polygon(zip([TSr['lon_b'][0,0].data,TSr['lon_b'][0,-1].data,TSr['lon_b'][-1,-1].data,
                     TSr['lon_b'][-1,0].data],[TSr['lat_b'][0,0].data,TSr['lat_b'][0,-1].data,
                                               TSr['lat_b'][-1,-1].data,TSr['lat_b'][-1,0].data]))

    gdom = gp.GeoSeries([dompoly])

    # assing defined polygon to a new dataframe
    pol_gpd= gp.GeoDataFrame()
    pol_gpd['geometry'] = None
    pol_gpd.loc[0,'geometry'] = dompoly
    pol_gpd.crs=regs2.crs #hope okay

    # crop shape file to domain, remove basins > 2 km2 
    result = gp.sjoin(regs2, pol_gpd, how='inner')#, op='within')
    resultn=result[result.areal_km2 > 2]
    resultn=resultn[resultn.stID!='2.11.0'] #only missing data 
    resultn.to_csv('basins_used.csv')

    # This generates the weights, takes time the first time only
    # ...
    savg = xe.SpatialAverager(TSr, resultn.geometry, geom_dim_name="stID",
                          filename=weightfile,reuse_weights=True)

    plotweights=False

    # this takes a lot of memory for all these basins, not recommended
    if plotweights:
        w = xr.DataArray(
            savg.weights.toarray().reshape(resultn.geometry.size, *TSr.lat.shape),
            dims=("stID", *TSr.lat.dims),
            coords=dict(stID=out.stID, **TSr.lon.coords),
        )

        plt.subplots_adjust(top=0.9)
        facets = w.plot(col="stID", col_wrap=6, aspect=2, vmin=0, vmax=0.05)
        facets.cbar.set_label("Averaging weights")
        plt.savefig('weights.png')

    return TSr, savg, resultn

#store the dfs in dict

# here the sparse matrix with grid cell weights for each ploygon is generated
# see xesmf doc: 
# https://pangeo-xesmf.readthedocs.io/en/latest/notebooks/Spatial_Averaging.html
# the matrix can be stored as a nc-file and re-used, so that the program runs 
# very fast once it is generated. 
# The weights are specific to the domain, projection, resolution, and to 
# options like 
#&nam_pgd_arrange_cover 
#  lwater_to_nature = .true. 
#  ltown_to_rock = .true.
#  which alters where RUNOFFC etc is defined 


def run(dictds, expname, TSr, savg, resultn):

    numruns=len(expnam)
    cachrunoff={}

    for i in np.arange(0,numruns):

        tmp1=dictds[expnam[i]]['RUNOFFC_ISBA'] +dictds[expnam[i]]['DRAINC_ISBA']
    
        #some experiments have accumlated runoff (i.e. do not use option
        # LRESETCUMUL = .true. in NAM_WRITE_DIAG_SURFn 
        # see https://www.umr-cnrm.fr/surfex/spip.php?article406 )
        # so need to deaccumulate
        #if expnam[i] in ['LDAS_FC']:
        #    tmp1=tmp1.diff(dim='time')
        tmp1=tmp1.to_dataset(name=expnam[i])
        tmp1['mask']=TSr.RUNOFFC_ISBA.isnull()
        tmp1['mask'].values=np.where(~TSr.RUNOFFC_ISBA.isnull(),1,0)
        out=savg(tmp1[expnam[i]])
        #out = out.assign_coords(roms_id=xr.DataArray(resultn.index.values, dims=("roms_id",)))
        out = out.assign_coords(stID=xr.DataArray(resultn["stID"], dims=("stID",)))
        print("out")
        print(out)
        cachrunoff[expnam[i]]=out.to_pandas()
        rainfRun= cachrunoff[expnam[i]]=out.to_pandas()
        print("rainfRun")
        print(rainfRun.keys())

        rainfRun[:].plot()
        #rainfRun.plot()
        plt.title(expnam[i], fontsize=12)
        plt.show()
        cachrunoff[expnam[i]].to_csv(expnam[i]+'.csv')


if __name__ == "__main__":

    expdirs='/lustre/storeB/users/josteinbl/sfx_data/'
    #forcingfile = '/lustre/storeB/project/nwp/H2O/wp4/FORCING/527_450_2020/FORCING_527_450_202010.nc'
    stationfile = '/lustre/storeB/project/nwp/H2O/wp4/RRdata/stations_Huang20.txt'
#    shapefiles = '/lustre/storeB/project/nwp/H2O/wp4/RRdata/shapefiles/CCM2/ccm21/WGS84_W2008.gdb'
    shapefiles = '/lustre/storeB/project/nwp/H2O/wp4/RRdata/shapefiles/latlonHyd/latlonhydorder/NVEData/Hydrologi/Hydrologi_TotalNedborfeltMalestasjon.shp'
    #shapefiles = '/lustre/storeB/project/nwp/H2O/wp4/RRdata/shapefiles/latlonREG/Nedborfelt/Nedborfelt_Vassdragsomr.shp'

    weightfile = '/lustre/storeB/users/josteinbl/sfx_data/LDAS_NOR/spatial_avg_4catchments_527nature2lim.nc'

    # folder where runs are, this could have been a long list of different exps.
    #expdirfs=['LDAS_NOR/', 'LDAS_NOR_eps_05/']
    expdirfs=['LDAS_FC/']
    # start of nc-file-output from runs
    #openoflsel201903.nc
    
    # Analysis mode:
    #ncstart=  ['openofl', 'openofl']
    ncstart=  ['openofl']
    # choose a name for legend in plots
    #expnam=['LDAS_NOR', 'LDAS_NOR_eps_05']
    expnam=['LDAS_FC']

    variables = ['RUNOFFC_ISBA', 'DRAINC_ISBA']

    myP, proj_string, fPREP, fPGD, FORCING, lcc = getDomainData(expdirs)
    stations = readStationList(stationfile)

    regs2 = readShapeFiles(shapefiles, stations)

    dictds = readModelData(expdirs, expdirfs, fPREP, fPGD, ncstart, expnam, variables)
    
    TSr, savg, resultn = generateWeights(dictds, fPREP, fPGD, FORCING, regs2, weightfile)

    run(dictds, expnam, TSr, savg, resultn)

    plot_catchments = False

    if plot_catchments == True:

        #for sid in regs2.stID:
        #for sid in stations["GauID"]:

        fig = plt.figure(figsize=(12,7))
        #gs = fig.add_gridspec(1, 1)
        #ax0=fig.add_subplot(gs[0,0], projection=ccrs.PlateCarree())
        #ax0.set_extent([FORCING.longitude.data.min()-15, FORCING.longitude.data.max()+15, FORCING.latitude.data.min()-15, FORCING.latitude.data.max()+15], crs=ccrs.PlateCarree())
        #ax0.add_feature(cfeature.LAND)
        #ax0.add_feature(cfeature.COASTLINE)
        #ax0.add_image(tiler, 7)
        #gl = ax0.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
        #              linewidth=2, color='gray', alpha=0.5, linestyle='--')
        #ax0.plot(FORCING.longitude[[0,0,-1,-1,0],[0,-1,-1,0,0]].values.diagonal(),FORCING.latitude[[0,0,-1,-1,0],[0,-1,-1,0,0]].values.diagonal(),color='blue', linewidth=2,transform=ccrs.Geodetic())
        #gl.top_labels = gl.right_labels = False
        ax0b = fig.add_subplot(projection=lcc)
        ax0b.set_extent([TSr.x.min(), TSr.x.max(), TSr.y.min(), TSr.y.max()], crs=lcc)
        gl = ax0b.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=0.8, color='gray', alpha=0.5, linestyle='--')
        gl.top_labels = gl.right_labels = False
        ax0b.add_image(tiler, 7)
        #regs2[regs2.stID==sid].boundary.plot(ax=ax0b,transform=ccrs.Geodetic(),color='r')            
        regs2.boundary.plot(ax=ax0b,transform=ccrs.Geodetic(),color='r')            
        #plt.title(sid+' '+str(resultn.areal_km2[resultn.stID==sid])+' km2')
        #regs2.boundary.plot(ax=ax0b,transform=ccrs.Geodetic(),color='r')
        plt.tight_layout()
        #plt.savefig("/lustre/storeB/users/josteinbl/TOPD/Figures/"+sid+'_snowDA_test.png')
        #plt.savefig('/lustre/storeB/users/josteinbl/TOPD/Figures/All_catchments_snowDA_test.png')
        #plt.close()
        #plt.show()
