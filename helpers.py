import numpy as np
from netCDF4 import Dataset
import os

def get_wrfout_timestamp_list(model_output_dir, domain = 1):

    DOM = str(domain).zfill(2)
    file_list = sorted(os.listdir(model_output_dir))
    time_stamp_list = [x[11:] for x in file_list if x[0:10] == ('wrfout_d' + DOM)]

    return time_stamp_list


def read_sfc_wind(fn):
    
    DS = Dataset(fn)
    WRF = {}
    WRF['lon'] = DS['XLONG'][:][0,:,:]
    WRF['lat'] = DS['XLAT'][:][0,:,:]
    WRF['u10'] = DS['U10'][:][0,:,:]
    WRF['v10'] = DS['V10'][:][0,:,:]
    WRF['wspd10'] = np.sqrt(WRF['u10'] * WRF['u10'] + WRF['v10'] * WRF['v10'])

    return WRF



def read_hycom_archv_3d(data_dir, time_stamp):
    
    ##
    ## Get grid.
    ##
    
    grid_file = data_dir + '/regional.grid.a'
    lon = awovispy.pyhycom.getField('plon', grid_file, np.nan)[0,:]
    lat = awovispy.pyhycom.getField('plat', grid_file, np.nan)[:,0]

    ##
    ## Get 3D data
    ##
    
    this_timetuple = dt.datetime.strptime(time_stamp, '%Y-%m-%d_%H:%M:%S').timetuple()
    this_doy = this_timetuple.tm_yday
    this_year = this_timetuple.tm_year
    this_hour = this_timetuple.tm_hour
    this_hycom_grid_file = (data_dir + "/regional.grid.a")
    #this_hycom_archv_file = (data_dir + "/archv." + str(this_year)
    #                         + "_" + str(this_doy).zfill(3) + "_" + str(this_hour).zfill(2) + ".a")
    this_hycom_time_step = time_stamp_to_hycom_time_step(time_stamp, "2017-07-17_00:00:00", 300, 12259296)
    this_hycom_archv_file = (data_dir + "/archv."  + str(int(this_hycom_time_step)).zfill(11) + ".a")
       
    u = awovispy.pyhycom.getField('u-vel', this_hycom_archv_file, np.nan)
    v = awovispy.pyhycom.getField('v-vel', this_hycom_archv_file, np.nan)
    t = awovispy.pyhycom.getField('temp', this_hycom_archv_file, np.nan)
    s = awovispy.pyhycom.getField('salin', this_hycom_archv_file, np.nan)
    
    dz = awovispy.pyhycom.getField('thknss', this_hycom_archv_file, np.nan) / 9806.0
    z = np.zeros(dz.shape)
    kdm, jdm, idm = z.shape

    for k in range(1, kdm):
        z[k,:,:] = z[k-1,:,:] + dz[k-1,:,:]

    mld = awovispy.pyhycom.getField('mix_dpth', this_hycom_archv_file, np.nan) / 9806.0

    return (lon, lat, z, t, s, u, v, mld)


def read_uwincm_waves(data_dir, time_stamp):
    pass
