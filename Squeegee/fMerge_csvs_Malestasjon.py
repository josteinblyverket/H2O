import pandas as pd
import glob, os
import geopandas as gp
import cartopy.crs as ccrs
import geopandas as gp
import matplotlib.pyplot as plt
#import contextily as cx #import folium
#from contextily import Place
plt.style.use('ggplot')
files_CTL = sorted(glob.glob('/ec/res4/scratch/sbjb/Projects/CERISE/Hydrology/RiverDischarge/CTL/SpatialAverager_ldas_runoff_Malestasjon_*.csv'))
files_mbr000 = sorted(glob.glob('/ec/res4/scratch/sbjb/Projects/CERISE/Hydrology/RiverDischarge/mbr000/SpatialAverager_ldas_runoff_Malestasjon_*.csv'))
#first date with new measurement stations is 202412121300
#print (files)

#read the shapefile
regs2 = gp.read_file('/ec/res4/scratch/sbjb/Projects/CERISE/Hydrology/Hydrologi_TotalNedborfeltMalestasjon.shp')
#read catchement metadata
mml=pd.read_csv('/ec/res4/scratch/sbjb/Projects/CERISE/Hydrology/stationmetasbig.csv',index_col=0)
#riverName                                       Lakselva
#percentLake
#qStartYear                                       2003.0
#qEndYear                                         2022.0
#culQm                                            25.5967 gul middelflor
#culQ5                                            33.0561 oransje 5 års
#culQ10                                           39.8306  
#3culQ20                                          46.636  
#culQ50                                           55.9146 rødt 50 års


#Automatisk beregnet statistikk for årene 2003 - 2022

#Gjentaksintervall: Gjennomsnittlig antall år mellom hver gang en hendelse (flom) vil inntreffe.

#Gult nivå (Middelflom): Gjennomsnittet av den største vannføringen hvert år.

#Oransje nivå (5-årsflom): Flom med gjentaksintervall 5 år. Det er 20% sannsynlighet, hvert år, for at en flom av denne størrelse vil overskrides.

#Rødt nivå (50-årsflom): Flom med gjentaksintervall 50 år. Det er 2% sannsynlighet, hvert år, for at en flom av denne størrelse vil overskrides.


df = pd.concat([pd.read_csv(fp,usecols=['stID','mrunoff'],index_col=0).
                T.assign(time=[pd.to_datetime(fp[-16:-4])]).
                set_index(['time']) for fp in files_CTL])

#add missing leading 0s
df.columns=df.columns.astype(str).str.zfill(3)

df_mbr = pd.concat([pd.read_csv(fp,usecols=['stID','mrunoff'],index_col=0).
                T.assign(time=[pd.to_datetime(fp[-16:-4])]).
                set_index(['time']) for fp in files_mbr000])

#add missing leading 0s
df_mbr.columns=df_mbr.columns.astype(str).str.zfill(3)


#first date with new measurement stations is 202412121300
#print (df)
#print(len(df.iloc[0,:]))
testplot=False

if testplot:
    # option mul by area_land and hr/sec
    # make m3/s:
    # 1E-3 m /hr |*1E6 m2 |*1/(60*60) hr/sec


    for i in range(0,len(df.iloc[0,:])):
      meta=regs2[regs2.stID==df.columns[i]]
      st = str(meta.stID.item())
      print(st)
#      if st == "189.3.0":
      if df.iloc[0,i] > 0.0:
          
        fig, axs = plt.subplot_mosaic("A;B", gridspec_kw=
                                      dict(height_ratios=[0.5, 0.5] ),
                                      figsize=(12, 8))

        meta=regs2[regs2.stID==df.columns[i]]
        if meta.stNavn.item()!=None:
            navn=meta.stNavn.item().replace(' ','').replace('.','')
        else:
            navn=meta.stID.item().replace(' ','').replace('.','_')

        meta.plot(ax=axs['A'])
        (1E-3*df.iloc[:,i]*meta.areal_km2.item()
                               *1E6/3600.).plot(ax=axs['B'])
        meta.plot(ax=axs['A'])
        (1E-3*df_mbr.iloc[:,i]*meta.areal_km2.item()
                               *1E6/3600.).plot(ax=axs['B'])
        if meta.stNavn.item()!=None:
            plt.title(meta.stID.item()+' '+navn+' '+
                      str(meta.areal_km2.item())+' km2')
        else:
            plt.title(meta.stID.item()+' '+str(meta.areal_km2.item())+' km2')
        plt.ylabel('m3/s')
        plt.legend(["CTL","mbr000"])
        #sys.exit()
        plt.savefig('/ec/res4/scratch/sbjb/Projects/CERISE/Hydrology/Figs/Malestasjon_mbr000_'+navn+'_'+
                    df.index[0].strftime('%Y%m%d%H_')+
                    df.index[-1].strftime('%Y%m%d%H.png'))
        #plt.show()
        plt.close()
                                  
plotwTSobs=True
#sys.exit()
if plotwTSobs:
    both=pd.DataFrame(index=df.columns.values)
    metadf=pd.DataFrame(index=both.index, columns=['name', 'area','glacier'])

    kgedf=pd.DataFrame(index=both.index, columns=['squegee'])
    biasdf=pd.DataFrame(index=both.index, columns=['squegee'])
    ccdf=pd.DataFrame(index=both.index, columns=['squegee'])
    mrdf=pd.DataFrame(index=both.index, columns=['squegee'])
    sdrdf=pd.DataFrame(index=both.index, columns=['squegee'])

    #no overlapping dates, obs startes later than runs
    #j=0
    print(both.index)
    
    for sid in both.index:
        #sid=resultn.index[k]
        meta=regs2[regs2.stID==sid]
        if meta.stNavn.item()!=None:
            navn=meta.stNavn.item().replace(' ','').replace('.','')
        else:
            navn=meta.stID.item().replace(' ','').replace('.','_')
        sim=df.loc[:,sid].resample('1D', closed='right',label='right').sum().rename('ldas')
        sim_mbr=df_mbr.loc[:,sid].resample('1D', closed='right',label='right').sum().rename('mbr000')

        simu=(1E-3*sim*meta.areal_km2.item()
                               *1E6/86400.)

        simu_mbr = (1E-3*sim_mbr*meta.areal_km2.item()
                               *1E6/86400.)

        print(glob.glob('/ec/res4/scratch/sbjb/Projects/CERISE/Hydrology/RR/'+sid+'_1001_*.csv'))
        fil=sorted(glob.glob('/ec/res4/scratch/sbjb/Projects/CERISE/Hydrology/RR/'+sid+'_1001_*.csv'))
        print(fil)
        if len(fil)==0:
            print('no obs ', str(sid))
           
        
        else:
            obs=pd.read_csv(fil[-1],index_col=0)
            obs.index=pd.to_datetime(obs.index)
            obs=obs.dropna()
            obspd=pd.DataFrame(index=pd.to_datetime(obs.index.strftime('%Y-%m-%d')), data=obs.values,columns=['Obs']) #*1000*86400/(resultn.areal_km2[resultn.stID==sid].iloc[0]*1E6),columns=['Obs'])  #
            obspd=obspd[ (obspd.index > pd.to_datetime('2018010123', format='%Y%m%d%H'))]#, errors='ignore'))]
            obspd=obspd[ (obspd.index < pd.to_datetime('2018080123', format='%Y%m%d%H'))]#, errors='ignore'))]
            both = obspd.join([simu,simu_mbr], how='outer')
            #both = both2.join(both2, how='outer')                                    

            
            
            ''' 
            if len(both>10):
                allboth[sid]=both
                kgestring=' kge: '
                biasstring=' bias: '
                ccstring=r' $ \rho $: '
                mrstring=r' $ \mu_s / \mu_o $: '
                sdrstring=r' $ \sigma_s / \sigma_o $: '
                for sim in ['CTRIP', 'squegee']:
                    KGE, cc, sdr, mr = KGEf(both[sim],both['Obs'])
                    kgedf.loc[sid,sim]=round(KGE,2)
                    mrdf.loc[sid,sim]=round(mr,2)
                    sdrdf.loc[sid,sim]=round(sdr,2)
                    ccdf.loc[sid,sim]=round(cc,2)          
                    biasdf.loc[sid,sim]=round(both.mean()[sim]-both.mean()['Obs'],1)
                    kgestring+=sim+' :'+str(kgedf.loc[sid,sim])+' '
                    ccstring+=sim+' :'+str(ccdf.loc[sid,sim])+' '
                    mrstring+=sim+' :'+str(mrdf.loc[sid,sim])+' '
                    sdrstring+=sim+' :'+str(sdrdf.loc[sid,sim])+' '
                    biasstring+=sim+' :'+str(biasdf.loc[sid,sim])+' '
                #print(sid,str(resultn.areal_km2[resultn.index==sid].iloc[0])+ ' km') 
            '''

            if len(both)>10:
                fig = plt.figure(figsize=(12,5))
                gs = fig.add_gridspec(1, 2,width_ratios=[1, 4])      
                #map = folium.Map(location=[13.406, 80.110], tiles="OpenStreetMap", zoom_start=9)
                ax2 = fig.add_subplot(gs[:,0], projection=ccrs.PlateCarree())
                ax2.set_extent([meta.bounds.minx.item()-2, meta.bounds.maxx.item()+3,meta.bounds.miny.item()-2, meta.bounds.maxy.item()+2], crs=ccrs.PlateCarree())
                gl = ax2.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                              linewidth=0.7, color='gray', alpha=0.5, linestyle='--')
                ax2.coastlines(resolution='10m',linewidth=0.5)
                
                #regs2.boundary.plot(ax=ax2,color='r')
                mm=meta.boundary.plot(ax=ax2, color='b')
#                cx.add_basemap(mm,zoom=6,crs=ccrs.PlateCarree())

                ax1 = fig.add_subplot(gs[:,1])
                both[['Obs','ldas', 'mbr000']].plot(style=['-' ,'--','--'],ax=ax1)
                plt.ylabel('m3/s')
                if (sid in mml.index):
                    plt.axhline(mml.loc[sid,'culQm'],color='yellow',ls='--')
                    plt.axhline(mml.loc[sid,'culQ5'],color='orange',ls='--')
                    plt.axhline(mml.loc[sid,'culQ50'],color='red',ls='--')
              
                if (meta.stNavn.item()!=None) & (sid in mml.index) :
                     plt.title(meta.stID.item()+' '+navn+', '+ mml.loc[sid,'riverName'] + ', '+
                                  str(int(round(meta.areal_km2.item(),0)))+' km2, lakes:' + str(mml.loc[sid,'percentLake'].item())+ '%, glaciers: '+str(mml.loc[sid,'percentGlacier'].item()) +'%')
                elif meta.stNavn.item()!=None:
                    plt.title(meta.stID.item()+' '+navn+' '+
                          str(int(round(meta.areal_km2.item(),0)))+' km2')
                else:
                    plt.title(meta.stID.item()+' '+str(meta.areal_km2.item())+' km2')

                #plt.grid(True)
                plt.tight_layout()
                plt.savefig('/ec/res4/scratch/sbjb/Projects/CERISE/Hydrology/Figs/Maalest_NVEid_'+str(sid)+'_testplotn.png')
                #plt.show()
                plt.close()
                print('obs for ', sid)
            
           
        #sys.exit()
