#!/usr/bin/env python
import numpy as np
from scipy.io import netcdf

def load_data(file_grid, file_ice, file_ocean):
    with netcdf.netcdf_file(file_grid, "r") as filegrid, netcdf.netcdf_file(file_ice, "r") as fileice:
        II = filegrid.variables["I_1to2"].data.copy()
        JJ = filegrid.variables["J_1to2"].data.copy()
        isocean = filegrid.variables["isocean"].data.copy()
        time_ice = fileice.variables["time"].data.copy()
        time_ice0 = np.size(time_ice) - 1
        bmb = fileice.variables["bmelt"].data.copy()
        fmb = fileice.variables["thk"].data.copy()
        fmb = np.where(fmb>0., 1., 0.)
        return II, JJ, isocean, time_ice, time_ice0, bmb, fmb




def main():
    pass

if __name__ == '__main__':
    # Program Flow Switches
    do_fortran_boundIJ = True
    write_latlon = False
    do_write_time = False
    # Filenames # TODO: Later as command line switches
    file_grid = "gridpoint_map_ice_to_ocean.nc"
    file_ice = "test_data/pism_ini.nc"
    file_ocean = "./ice2oce.nc"

    # Physical Parameters
    latent_heat = 333700.       # J/kg, PG: Of what?
    rho_ice = 910.               # kg/m3
    rho_oce = 1026.0
    salinity_seaice = 5.

    # Load the data
    II, JJ, isocean, time_ice, time_ice0, bmb, fmb = load_data(file_grid, file_ice, file_ocean)
    ilen1, jlen1 = np.shape(bmb[0, :, :])
    ilen2, jlen2 = np.shape(isocean)
    tlen = np.size(time_ice)

    if tlen>=1:
        sizeIce2Oce = (tlen, ilen2, jlen2)
    else:
        sizeIce2Oce = (ilen2, jlen2)
    
    if do_fortran_boundIJ:
        II -= 1
        JJ -= 1
    bmb_oce = np.zeros(sizeIce2Oce)
    fmb_oce = np.zeros(sizeIce2Oce)

    bhb_oce = np.zeros(sizeIce2Oce)
    fhb_oce = np.zeros(sizeIce2Oce)

    if tlen >= 1:
        for it in range(tlen):
            for i in range(ilen1):
                for j in range(jlen1):
                    io = II[i, j]
                    jo = JJ[i, j]
                    bmb_oce[it, io, jo] += bmb[it, i, j]
                    fmb_oce[it, io, jo] += fmb[it, i, j]
    else:
        for i in range(ilen1):
            for j in range(jlen1):
                io = II[i, j]
                jo = JJ[i, j]
                bmb_oce[io, jo] += bmb[i, j]
                fmb_oce[io, jo] += fmb[i, j]
    # Convert latent heat to a heat flux by assuming everything must melt completely:
    conversion_factor = latent_heat * rho_ice
    bhb_oce = bmb_oce * conversion_factor
    fhb_oce = fmb_oce * conversion_factor

    ### Save the data
    print("Saving output...")
    fout = netcdf.netcdf_file(file_ocean, mode="w")

    #
    # Global NetCDF attributes
    #
    print('    Global attributes')
    fout.authors = 'Paul Gierz'
    fout.instiution = "Alfred Wegener Institute for Polar and Marine Research"
    fout.address = "Bussesstrasse 24, 27570 Bremerhaven, Germany"
    fout.web = "http://www.awi.de"
    import os   as os
    fout.working_directory = os.getcwd()
    fout.username = os.environ["USER"]
    #?fout.hostname = os.environ["HOST"]
    #fout.hosttype = os.environ["HOSTTYPE"]
    fout.uname = os.popen('uname -a').read()
    #fout.gridfile1=fileice
    #fout.gridfile2=fileoce


    #fout.gridfile1=fileice
    #fout.gridfile2=fileoce

    #
    # Definition of the dimensions
    #
    print('    Create dimension')
    if ( tlen >= 1 ):
        fout.createDimension('time',tlen)
    else :
        fout.createDimension('time',tlen)
    fout.createDimension('x1',ilen1)
    fout.createDimension('y1',jlen1)
    fout.createDimension('x2',ilen2)
    fout.createDimension('y2',jlen2)
    #fout.createDimension('z',klen)
    fout.createDimension('one',(1))


    #
    # Define output variables and specify its attributes
    #
    print('    Define variables')

    ## x-, y-arrays
    x1_var       = fout.createVariable('x1', 'd', ('x1', ))
    x1_var.standard_name = 'projection_x_coordinate'
    x1_var.long_name     = 'X-coordinate in Cartesian system' ;
    #x1_var.unit          = 'm'
    y1_var       = fout.createVariable('y1', 'd', ('y1', ))
    y1_var.standard_name = 'projection_y_coordinate'
    y1_var.long_name     = 'Y-coordinate in Cartesian system' ;
    #y1_var.unit          = 'm'

    x2_var       = fout.createVariable('x2', 'd', ('x2', ))
    x2_var.standard_name = 'projection_x_coordinate'
    x2_var.long_name     = 'X-coordinate in Cartesian system' ;
    #x2_var.unit          = 'm'
    y2_var       = fout.createVariable('y2', 'd', ('y2', ))
    y2_var.standard_name = 'projection_y_coordinate'
    y2_var.long_name     = 'Y-coordinate in Cartesian system' ;
    #y2_var.unit          = 'm'


    if ( write_latlon ):
        lat1_var       = fout.createVariable('lat1', 'd', ('x1', 'y1', ))
        lat1_var.standard_name = 'latitude'
        lat1_var.long_name     = 'latitude,  grid 1'
        lat1_var.unit          = 'degree'
        lon1_var       = fout.createVariable('lon1', 'd', ('x1', 'y1', ))
        lon1_var.standard_name = 'longitude'
        lon1_var.long_name     = 'longitude, grid 1' ;
        lon1_var.unit          = 'degree'


        lat2_var       = fout.createVariable('lat2', 'd', ('x2', 'y2', ))
        lat2_var.standard_name = 'latitude'
        lat2_var.long_name     = 'latitude,  grid 2'
        lat2_var.unit          = 'degree'
        lon2_var       = fout.createVariable('lon2', 'd', ('x2', 'y2', ))
        lon2_var.standard_name = 'longitude'
        lon2_var.long_name     = 'longitude, grid 2' ;
        lon2_var.unit          = 'degree'


    isocean_var       = fout.createVariable('isocean', 'd', ('x2', 'y2', ))
    isocean_var.long_name   = 'ocean mask (1:ocean, 0:else)' ;
    isocean_var.coordinates ='lat2 lon2'

    #
    if ( tlen >= 1 ):
        bmr_var       = fout.createVariable('bmr', 'd', ('time', 'x2', 'y2', ))
    else :
        bmr_var       = fout.createVariable('bmr', 'd', ('x2', 'y2', ))
    bmr_var.long_name   ='basal_melting_rate'
    bmr_var.units       ='m3/s'
    bmr_var.coordinates ='lat2 lon2'
    #bmr_var.reference_density=dens_ice

    if ( tlen >= 1 ):
        fmr_var       = fout.createVariable('fmr', 'd', ('time', 'x2', 'y2', ))
    else :
        fmr_var       = fout.createVariable('fmr', 'd', ('x2', 'y2', ))    
    fmr_var.long_name   ='frontal_melting_rate'
    fmr_var.units       ='m3/s'
    fmr_var.coordinates ='lat2 lon2'

    #
    #
    if ( tlen >= 1 ):
        bhr_var       = fout.createVariable('bhr', 'd', ('time', 'x2', 'y2', ))
    else :
        bhr_var       = fout.createVariable('bhr', 'd', ('x2', 'y2', ))
    bhr_var.long_name   ='basal_heating_rate'
    bhr_var.units       ='Jm/s'
    bhr_var.coordinates ='lat2 lon2'

    if ( tlen >= 1 ):
        fhr_var       = fout.createVariable('fhr', 'd', ('time', 'x2', 'y2', ))
    else :
        fhr_var       = fout.createVariable('fhr', 'd', ('x2', 'y2', ))    
    fhr_var.long_name   ='frontal_heating_rate'
    fhr_var.units       ='Jm/s'
    fhr_var.coordinates ='lat2 lon2'
    #fhr_var.reference_density=dens_ice

    #
    #
    #
    print('    Populate variables')

    # x-, y-arrarys
    x1_var[:]   = range(ilen1); y1_var[:]   = range(jlen1)
    x2_var[:]   = range(ilen2); y2_var[:]   = range(jlen2)

    if ( write_latlon ):
        lat1_var[:] = lat1
        lon1_var[:] = lon1

        lat2_var[:] = lat2
        lon2_var[:] = lon2

    isocean_var[:]=isocean

    bmr_var[:] = bmb_oce
    fmr_var[:] = fmb_oce

    bhr_var[:] = bhb_oce
    fhr_var[:] = fhb_oce
    ################
    #
    # Close the files
    #
    print "Close all open files"


    fout.close()

    # -----------
    #
    # Some information
    #
    print("")
    if ( do_write_time ) :
        print(' Done at '+time.ctime(time.time()))
    else :
        print(' Done')
