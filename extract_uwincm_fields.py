import numpy as np
from netCDF4 import Dataset
import sys
import datetime as dt
from pyhycom import pyhycom

## Import local helpers module (helpers.py).
import helpers


## Using compression, the data size of the Alaska awo_2014110700_fnl_4km_nudging run was
## reduced from 24 GB to ~7 GB. However, it takes LONG to write the file in the last step.
do_compression = True
#do_compression = False

## Specify the compression level here (1-9).
## The default (4) yielded a file size of 6.8 GB.
## However, it was very slow to write and to view in ncview.
## complevel = 1 yielded a file size of 6.9 GB and is much faster.
complevel = 1

fmt_wrf='%Y-%m-%d_%H:%M:%S'

def time_stamp_to_datetime(time_stamp):
    return dt.datetime.strptime(time_stamp, '%Y-%m-%d_%H:%M:%S')


if len(sys.argv) < 2:
    
    print('Usage: extract_uwincm_fields /path/to/uwincm/run')

else:

    model_run_dir = sys.argv[1]
    
    time_stamp_list = helpers.get_wrfout_timestamp_list(model_run_dir, domain=1)

    print(time_stamp_list)

    
    """
    WRF Data
    """
    WRF = {}
    WRF['dt'] = []
    WRF2 = {}
    WRF2['dt'] = []

    
    ## D01
    
    tt = -1    
    for this_time_stamp in time_stamp_list:
        
        this_dt = time_stamp_to_datetime(this_time_stamp)
        print(this_dt)

        try:
            fn = (model_run_dir + '/wrfout_d01_' + this_time_stamp)
            print(fn)
            F = helpers.read_sfc_wind(fn)
            DS=Dataset(fn)
            F['land_mask'] = DS['LANDMASK'][:]
            F['taux'] = DS['TAUX_ESMF'][:]
            F['tauy'] = DS['TAUY_ESMF'][:]
            F['psfc'] = DS['PSFC'][:]
            DS.close()
            
            if not 'lon' in WRF:
                WRF['lon'] = F['lon']
                WRF['lat'] = F['lat']
                WRF['land_mask'] = F['land_mask']
                wrf_shape = F['lon'].shape
                WRF['u10'] = np.zeros((len(time_stamp_list), wrf_shape[0], wrf_shape[1]))
                WRF['v10'] = np.zeros((len(time_stamp_list), wrf_shape[0], wrf_shape[1]))
                WRF['wspd10'] = np.zeros((len(time_stamp_list), wrf_shape[0], wrf_shape[1]))
                WRF['u_wind_stress'] = np.zeros((len(time_stamp_list), wrf_shape[0], wrf_shape[1]))
                WRF['v_wind_stress'] = np.zeros((len(time_stamp_list), wrf_shape[0], wrf_shape[1]))
                WRF['psfc'] = np.zeros((len(time_stamp_list), wrf_shape[0], wrf_shape[1]))

            tt += 1
            WRF['dt'].append(this_dt)
            WRF['u10'][tt,:,:] = F['u10'][:,:]
            WRF['v10'][tt,:,:] = F['v10'][:,:]
            WRF['wspd10'][tt,:,:] = F['wspd10'][:,:]
            WRF['u_wind_stress'][tt,:,:] = F['taux'][:,:]
            WRF['v_wind_stress'][tt,:,:] = F['tauy'][:,:]
            WRF['psfc'][tt,:,:] = F['psfc'][:,:]

        except:
            print('WRF d01: There was an issue with this time.')


    ## D02

    tt = -1    
    for this_time_stamp in time_stamp_list:
        
        this_dt = time_stamp_to_datetime(this_time_stamp)
        print(this_dt)
        
        try:
            fn2 = (model_run_dir + '/wrfout_d02_' + this_time_stamp)
            print(fn2)
            F2 = helpers.read_sfc_wind(fn2)
            DS2=Dataset(fn2)
            F2['land_mask'] = DS2['LANDMASK'][:]
            F2['taux'] = DS2['TAUX_ESMF'][:]
            F2['tauy'] = DS2['TAUY_ESMF'][:]
            F2['psfc'] = DS2['PSFC'][:]
            DS2.close()
    
            if not 'lon' in WRF2:
                WRF2['lon'] = F2['lon']
                WRF2['lat'] = F2['lat']
                WRF2['land_mask'] = F2['land_mask']
                wrf_shape2 = F2['lon'].shape
                WRF2['u10'] = np.zeros((len(time_stamp_list), wrf_shape2[0], wrf_shape2[1]))
                WRF2['v10'] = np.zeros((len(time_stamp_list), wrf_shape2[0], wrf_shape2[1]))
                WRF2['wspd10'] = np.zeros((len(time_stamp_list), wrf_shape2[0], wrf_shape2[1]))
                WRF2['u_wind_stress'] = np.zeros((len(time_stamp_list), wrf_shape2[0], wrf_shape2[1]))
                WRF2['v_wind_stress'] = np.zeros((len(time_stamp_list), wrf_shape2[0], wrf_shape2[1]))
                WRF2['psfc'] = np.zeros((len(time_stamp_list), wrf_shape2[0], wrf_shape2[1]))

            tt += 1
            WRF2['dt'].append(this_dt)
            WRF2['u10'][tt,:,:] = F2['u10'][:,:]
            WRF2['v10'][tt,:,:] = F2['v10'][:,:]
            WRF2['wspd10'][tt,:,:] = F2['wspd10'][:,:]
            WRF2['u_wind_stress'][tt,:,:] = F2['taux'][:,:]
            WRF2['v_wind_stress'][tt,:,:] = F2['tauy'][:,:]
            WRF2['psfc'][tt,:,:] = F2['psfc'][:,:]

        except:
            print('WRF d02: There was an issue with this time.')


    ## HYCOM

    HY = {}
    HY['dt'] = []

    tt = -1    
    for this_time_stamp in time_stamp_list:
        
        this_dt = time_stamp_to_datetime(this_time_stamp)
        print(this_dt)
        hycom_format_time_stamp = this_dt.strftime('%Y_%j_%H')
        
        try:
            
            fnhy = (model_run_dir + '/archv.' + hycom_format_time_stamp + '.a')
            print(fnhy)
            u = pyhycom.getField('u-vel', fnhy, np.nan, [0])
            v = pyhycom.getField('v-vel', fnhy, np.nan, [0])

            if not 'lon' in HY:
                grid_file = model_run_dir + '/regional.grid.a'
                HY['lon'] = pyhycom.getField('plon', grid_file, np.nan)
                HY['lat'] = pyhycom.getField('plat', grid_file, np.nan)
                dims = pyhycom.getDims(fnhy)
                HY['bathy'] = pyhycom.getBathymetry((model_run_dir + '/regional.depth.a'), [dims[1], dims[2]], np.NaN)
                HY['u'] = np.zeros((len(time_stamp_list), dims[1], dims[2]))
                HY['v'] = np.zeros((len(time_stamp_list), dims[1], dims[2]))

            tt += 1
            HY['dt'].append(this_dt)
            HY['u'][tt,:,:] = u[0,:,:]
            HY['v'][tt,:,:] = v[0,:,:]
                
        except:
            print('HYCOM: There was an issue with this time.')

    ## UMWM

    UMWM = {}
    UMWM['dt'] = []

    tt = -1    
    for this_time_stamp in time_stamp_list:
        
        this_dt = time_stamp_to_datetime(this_time_stamp)
        print(this_dt)
        
        try:
            
            fnumwm = (model_run_dir + '/umwmout_' + this_time_stamp + '.nc')
            print(fnumwm)
            DS = Dataset(fnumwm)
            
            if not 'lon' in UMWM:
                UMWM['lon'] = DS['lon'][:]
                UMWM['lat'] = DS['lat'][:]
                umwm_shape = UMWM['lon'].shape
                UMWM['bathy'] = DS['depth'][:]
                UMWM['seamask'] = DS['seamask'][:]
                UMWM['swh'] = np.zeros((len(time_stamp_list), umwm_shape[1], umwm_shape[2]))
                UMWM['u_sfc_curr'] = np.zeros((len(time_stamp_list), umwm_shape[1], umwm_shape[2]))
                UMWM['v_sfc_curr'] = np.zeros((len(time_stamp_list), umwm_shape[1], umwm_shape[2]))

            tt += 1
            UMWM['dt'].append(this_dt)
            UMWM['swh'][tt,:,:] = DS['swh'][:]
            UMWM['u_sfc_curr'][tt,:,:] = DS['uc'][:]
            UMWM['v_sfc_curr'][tt,:,:] = DS['vc'][:]            
            DS.close()
            
        except:
            print('UMWM: There was an issue with this time.')
       
    ## Set WRF lon to be 0 - 360
    WRF['lon'][WRF['lon'] < 0.0] += 360.0
    WRF2['lon'][WRF2['lon'] < 0.0] += 360.0

    
    """
    NetCDF Output
    """

    fn_out = 'alaska_awo_2014110700_fnl_4km_nudging.nc4'
    print('Writing NetCDF4 output to: ' + fn_out)
    if do_compression:
        print('You specified compression, so this may take awhile!')
    
    FOUT = Dataset(fn_out, 'w', format='NETCDF4')

    ## Dimensions
    print('Dimensions', flush=True)
    FOUT.createDimension('wrf_time_d01',size=len(WRF['dt']))
    FOUT.createDimension('wrf_we_d01',size=wrf_shape[1])
    FOUT.createDimension('wrf_sn_d01',size=wrf_shape[0])

    FOUT.createDimension('wrf_time_d02',size=len(WRF2['dt']))
    FOUT.createDimension('wrf_we_d02',size=wrf_shape2[1])
    FOUT.createDimension('wrf_sn_d02',size=wrf_shape2[0])

    FOUT.createDimension('hycom_time',size=len(HY['dt']))
    FOUT.createDimension('hycom_x',size=dims[2])
    FOUT.createDimension('hycom_y',size=dims[1])

    FOUT.createDimension('umwm_time',size=len(UMWM['dt']))
    FOUT.createDimension('umwm_x',size=umwm_shape[2])
    FOUT.createDimension('umwm_y',size=umwm_shape[1])
    
    ## Variables
    print('Variables', flush=True)
    FOUT.createVariable('time_wrf_d01', 'd', ('wrf_time_d01',),zlib=do_compression,complevel=complevel)
    FOUT.createVariable('lon_wrf_d01', 'd', ('wrf_sn_d01','wrf_we_d01',),zlib=do_compression,complevel=complevel)
    FOUT.createVariable('lat_wrf_d01', 'd', ('wrf_sn_d01','wrf_we_d01',),zlib=do_compression,complevel=complevel)
    FOUT.createVariable('land_mask_wrf_d01', 'i', ('wrf_sn_d01','wrf_we_d01',),zlib=do_compression,complevel=complevel)
    FOUT.createVariable('u10_wrf_d01', 'd', ('wrf_time_d01','wrf_sn_d01','wrf_we_d01',),zlib=do_compression,complevel=complevel)
    FOUT.createVariable('v10_wrf_d01', 'd', ('wrf_time_d01','wrf_sn_d01','wrf_we_d01',),zlib=do_compression,complevel=complevel)
    FOUT.createVariable('wspd10_wrf_d01', 'd', ('wrf_time_d01','wrf_sn_d01','wrf_we_d01',),zlib=do_compression,complevel=complevel)
    FOUT.createVariable('u_wind_stress_wrf_d01', 'd', ('wrf_time_d01','wrf_sn_d01','wrf_we_d01',),zlib=do_compression,complevel=complevel)
    FOUT.createVariable('v_wind_stress_wrf_d01', 'd', ('wrf_time_d01','wrf_sn_d01','wrf_we_d01',),zlib=do_compression,complevel=complevel)
    FOUT.createVariable('psfc_wrf_d01', 'd', ('wrf_time_d01','wrf_sn_d01','wrf_we_d01',),zlib=do_compression,complevel=complevel)
    
    FOUT.createVariable('time_wrf_d02', 'd', ('wrf_time_d02',),zlib=do_compression,complevel=complevel)
    FOUT.createVariable('lon_wrf_d02', 'd', ('wrf_sn_d02','wrf_we_d02',),zlib=do_compression,complevel=complevel)
    FOUT.createVariable('lat_wrf_d02', 'd', ('wrf_sn_d02','wrf_we_d02',),zlib=do_compression,complevel=complevel)
    FOUT.createVariable('land_mask_wrf_d02', 'i', ('wrf_sn_d02','wrf_we_d02',),zlib=do_compression,complevel=complevel)
    FOUT.createVariable('u10_wrf_d02', 'd', ('wrf_time_d02','wrf_sn_d02','wrf_we_d02',),zlib=do_compression,complevel=complevel)
    FOUT.createVariable('v10_wrf_d02', 'd', ('wrf_time_d02','wrf_sn_d02','wrf_we_d02',),zlib=do_compression,complevel=complevel)
    FOUT.createVariable('wspd10_wrf_d02', 'd', ('wrf_time_d02','wrf_sn_d02','wrf_we_d02',),zlib=do_compression,complevel=complevel)
    FOUT.createVariable('u_wind_stress_wrf_d02', 'd', ('wrf_time_d02','wrf_sn_d02','wrf_we_d02',),zlib=do_compression,complevel=complevel)
    FOUT.createVariable('v_wind_stress_wrf_d02', 'd', ('wrf_time_d02','wrf_sn_d02','wrf_we_d02',),zlib=do_compression,complevel=complevel)
    FOUT.createVariable('psfc_wrf_d02', 'd', ('wrf_time_d02','wrf_sn_d02','wrf_we_d02',),zlib=do_compression,complevel=complevel)
    
    FOUT.createVariable('time_hycom', 'd', ('hycom_time',),zlib=do_compression,complevel=complevel)
    FOUT.createVariable('lon_hycom', 'd', ('hycom_y','hycom_x',),zlib=do_compression,complevel=complevel)
    FOUT.createVariable('lat_hycom', 'd', ('hycom_y','hycom_x',),zlib=do_compression,complevel=complevel)
    FOUT.createVariable('bathy_hycom', 'd', ('hycom_y','hycom_x',),zlib=do_compression,complevel=complevel)
    FOUT.createVariable('u_top_hycom', 'd', ('hycom_time','hycom_y','hycom_x',),zlib=do_compression,complevel=complevel)
    FOUT.createVariable('v_top_hycom', 'd', ('hycom_time','hycom_y','hycom_x',),zlib=do_compression,complevel=complevel)

    FOUT.createVariable('time_umwm', 'd', ('umwm_time',),zlib=do_compression,complevel=complevel)
    FOUT.createVariable('lon_umwm', 'd', ('umwm_y','umwm_x',),zlib=do_compression,complevel=complevel)
    FOUT.createVariable('lat_umwm', 'd', ('umwm_y','umwm_x',),zlib=do_compression,complevel=complevel)
    FOUT.createVariable('bathy_umwm', 'd', ('umwm_y','umwm_x',),zlib=do_compression,complevel=complevel)
    FOUT.createVariable('seamask_umwm', 'd', ('umwm_y','umwm_x',),zlib=do_compression,complevel=complevel)
    FOUT.createVariable('swh_umwm', 'd', ('umwm_time','umwm_y','umwm_x',),zlib=do_compression,complevel=complevel)
    FOUT.createVariable('u_sfc_current_umwm', 'd', ('umwm_time','umwm_y','umwm_x',),zlib=do_compression,complevel=complevel)
    FOUT.createVariable('v_sfc_current_umwm', 'd', ('umwm_time','umwm_y','umwm_x',),zlib=do_compression,complevel=complevel)

    # Values
    print('Values: WRF d01', flush=True)
    FOUT['time_wrf_d01'][:] = [(x - dt.datetime(1970,1,1,0,0,0)).total_seconds() / 3600.0 for x in WRF['dt']]
    FOUT['lon_wrf_d01'][:] = WRF['lon']
    FOUT['lat_wrf_d01'][:] = WRF['lat']
    FOUT['land_mask_wrf_d01'][:] = WRF['land_mask']
    FOUT['u10_wrf_d01'][:] = WRF['u10']
    FOUT['v10_wrf_d01'][:] = WRF['v10']
    FOUT['wspd10_wrf_d01'][:] = WRF['wspd10']
    FOUT['u_wind_stress_wrf_d01'][:] = WRF['u_wind_stress']
    FOUT['v_wind_stress_wrf_d01'][:] = WRF['v_wind_stress']
    FOUT['psfc_wrf_d01'][:] = WRF['psfc']

    print('Values: WRF d02', flush=True)
    FOUT['time_wrf_d02'][:] = [(x - dt.datetime(1970,1,1,0,0,0)).total_seconds() / 3600.0 for x in WRF2['dt']]
    FOUT['lon_wrf_d02'][:] = WRF2['lon']
    FOUT['lat_wrf_d02'][:] = WRF2['lat']
    FOUT['land_mask_wrf_d02'][:] = WRF2['land_mask']
    FOUT['u10_wrf_d02'][:] = WRF2['u10']
    FOUT['v10_wrf_d02'][:] = WRF2['v10']
    FOUT['wspd10_wrf_d02'][:] = WRF2['wspd10']
    FOUT['u_wind_stress_wrf_d02'][:] = WRF2['u_wind_stress']
    FOUT['v_wind_stress_wrf_d02'][:] = WRF2['v_wind_stress']
    FOUT['psfc_wrf_d02'][:] = WRF2['psfc']

    print('Values: HYCOM', flush=True)
    FOUT['time_hycom'][:] = [(x - dt.datetime(1970,1,1,0,0,0)).total_seconds() / 3600.0 for x in HY['dt']]
    FOUT['lon_hycom'][:] = HY['lon']
    FOUT['lat_hycom'][:] = HY['lat']
    FOUT['bathy_hycom'][:] = HY['bathy']
    FOUT['u_top_hycom'][:] = HY['u']
    FOUT['v_top_hycom'][:] = HY['v']

    print('Values: UMWM', flush=True)
    FOUT['time_umwm'][:] = [(x - dt.datetime(1970,1,1,0,0,0)).total_seconds() / 3600.0 for x in UMWM['dt']]
    FOUT['lon_umwm'][:] = UMWM['lon']
    FOUT['lat_umwm'][:] = UMWM['lat']
    FOUT['bathy_umwm'][:] = UMWM['bathy']
    FOUT['seamask_umwm'][:] = UMWM['seamask']
    FOUT['swh_umwm'][:] = UMWM['swh']
    FOUT['u_sfc_current_umwm'][:] = UMWM['u_sfc_curr']
    FOUT['v_sfc_current_umwm'][:] = UMWM['v_sfc_curr']

    
    ## Attributes
    for field in ['time_wrf_d01','time_wrf_d02','time_hycom','time_umwm']:
        FOUT[field].setncattr('units','hours since 1970-1-1')

    for field in ['u10_wrf_d01','v10_wrf_d01','wspd10_wrf_d01'
                  ,'u10_wrf_d02','v10_wrf_d02','wspd10_wrf_d02'
                  ,'u_top_hycom','v_top_hycom'
                  ,'u_sfc_current_umwm','v_sfc_current_umwm']:
        FOUT[field].setncattr('units','m/s')

    for field in ['u_wind_stress_wrf_d01','v_wind_stress_wrf_d01','u_wind_stress_wrf_d02','v_wind_stress_wrf_d02']:
        FOUT[field].setncattr('units','N/m2')

    for field in ['swh_umwm','bathy_hycom','bathy_umwm']: 
        FOUT[field].setncattr('units','m')

    FOUT['psfc_wrf_d01'].setncattr('units','Pa')
    FOUT['psfc_wrf_d02'].setncattr('units','Pa')
        
    FOUT.close()
    print('Done.)', flush=True)
    
