# Compute from a given (smaller ice sheet) Lat- and Lon-fields, fields #1, for
# a (larger) other lat-, lon-fields, field #2, of the corresponding grid point
# indexes.
#
# This routine uses Grid1to2-File to PERFORM THE TRANSFORMATION
#
#
# This stuff is needed to determine for given ice sheet model grids (PISM) the
# corresponding ocean grid points.

file_grid    = 'Grid1to2.nc'
file_ice     = 'pism_ini.nc'

file_ocean   = 'ice2oce.nc'




#
# Switches that control the behavior of the code
#
do_plot         = False # Switch off: False, Default: False
#do_plot         = True
do_write_time   = True  # Switch off: False, Default: True
do_fortran_boundIJ = True       # Index start at one/1 (=Fortran): True; zero/0: False
write_latlon    = True



#
# Important physical parameters
#
latent_heat= 333700.#J/kg
rho_ice    =  910.0 #kg/m3
rho_oce    = 1026.0 #kg/m3
salinie_seaice = 5  #g/kg = psu

print("")
print("Import packages / modules")
import numpy as N
import Scientific.IO.NetCDF as NC
if ( do_plot ):         import matplotlib.pyplot    as PLT
if ( do_write_time ) :  import time as time
    
print("       Done import packages / modules")

if ( do_write_time ) :
    print('Start at '+time.ctime(time.time()))



#=============================================    
#
# -------------- main program ----------------
#


# Write logging information
#
print("What will I do")
print("  do_plot                : "+str(do_plot))
print("  do_write_time          : "+str(do_write_time))
print("  do_fortran_boundIJ     : "+str(do_fortran_boundIJ))
if ( do_fortran_boundIJ ) :
    print(" \t Index start at one/1 (Fortran-bound)")
else :
    print(" \t Index start at zero/0 (C-bound, python-bound)")

print("  write_latlon           : "+str(write_latlon))

#Future: do_parallel     = False #default !!
#Future: if ( do_parallel ) :
#Future:     import pp as pp
#
#Future: print("  do_parallel            : "+str(do_parallel))
#Future: #ppservers = ("*", )
#Future: ppservers = ()
#Future: if (do_parallel):
#Future:     print " Do have started outside of the python script 'ppserver -s lala007'"
#Future:     job_server = pp.Server(ppservers=ppservers, secret="lala007")
#Future:     #job_server = pp.Server(ncpus=2)
#Future:     print(" \t Starting pp with", job_server.get_ncpus(), "workers")

#
#

print("")

# -------------- 
#
# Open the files and access the data
#
print("Load data")

# -----------
# -----------
#
# IJ grid point 
#
# -----------
#
# ocean data
#
print(" Grid data file: "+file_grid)
filegrid  = NC.NetCDFFile(file_grid,mode="r")
II      = filegrid.variables['I_1to2'].getValue()
JJ      = filegrid.variables['J_1to2'].getValue()
#WRONG II      = filegrid.variables['TI_1to2'].getValue()
#WRONG JJ      = filegrid.variables['TJ_1to2'].getValue()
isocean = filegrid.variables['isocean'].getValue()
if ( write_latlon ):
    lat1     = filegrid.variables['lat1'].getValue()
    lon1     = filegrid.variables['lon1'].getValue()
    lat2     = filegrid.variables['lat2'].getValue()
    lon2     = filegrid.variables['lon2'].getValue()



#
# Ice melt fields
#
print(" Ice model data file  : "+file_ice)
fileice  = NC.NetCDFFile(file_ice,mode="r")
time_ice = fileice.variables['time'].getValue()
time_ice0= N.size(time_ice) - 1
#lat1     = fileice.variables['lat'].getValue()
#lon1     = fileice.variables['lon'].getValue()
#xpism    = fileice.variables['x'].getValue()
#ypism    = fileice.variables['y'].getValue()

#fmb     = fileice.variables['fmb'].getValue()
#bmb     = fileice.variables['bmb'].getValue()
bmb     = fileice.variables['bmelt'].getValue()
fmb     = fileice.variables['thk'].getValue()
fmb = N.where(fmb>0., 1., 0.)
print "********** Manipulated for tests fmb ***************"


#=============================================
#
#
# Deduced variables
#

[ilen1, jlen1]    = N.shape(bmb[0, :, :])

[ilen2, jlen2]   = N.shape(isocean) #I_1to2
#WRONG [jlen2, ilen2]   = N.shape(isocean) #TI_1to2

tlen            = N.size(time_ice)

# FIXME: Assumption about the NEMO/ocean data file
#   1. time is first index in the 
#   2. lat/lon are the second and third indexes
if ( tlen >= 1 ):
    sizeIce2Oce = (tlen, ilen2, jlen2)
else :
    sizeIce2Oce = (      ilen2, jlen2)


#
# Switch from zero-bound index scheme (C, python) to one-bound index scheme (fortran) ?
#
if ( do_fortran_boundIJ ) :
    print(" Index start at one/1 (Fortran-bound)")
    II = II - 1
    JJ = JJ - 1
else :
    print(" Index start at zero/0 (C-bound, python-bound)")

print("")



bmb_oce = N.zeros(sizeIce2Oce)
fmb_oce = N.zeros(sizeIce2Oce)

bhb_oce = N.zeros(sizeIce2Oce)
fhb_oce = N.zeros(sizeIce2Oce)

if ( tlen >= 1 ):
    for it in range(tlen) :
        for i in range(ilen1) :
            for j in range(jlen1) :
                io = II[i, j]
                jo = JJ[i, j]
                #bmb_oce[it, io, jo] = bmb_oce[it, io, jo] + bmb[it, i, j]
                bmb_oce[it, io, jo] += bmb[it, i, j]
                fmb_oce[it, io, jo] += fmb[it, i, j]
                
else :                
  for i in range(ilen1) :
    for j in range(jlen2) :
        io = II[i, j]
        jo = JJ[i, j]
        bmb_oce[io, jo] += bmb[i, j]
        fmb_oce[io, jo] += fmb[i, j]


#FIXME: Melting rate into heat rate
print ("Dummy calculation of heat rate !!!!!!!!!!!!!!!!!!!!")


conversion_factor = latent_heat * rho_ice
bhb_oce = bmb_oce * conversion_factor
fhb_oce = fmb_oce * conversion_factor

# -------------- Output file ----------------
#
print("Save data")

# -----------
#
# ice sheet/shelf data
#
print(" Ice-ocean interaction data file  : "+file_ocean)
fileout = NC.NetCDFFile(file_ocean,mode='w')  # Open the output file

#
# Global NetCDF attributes
#
print('    Global attributes')
fileout.authors = 'Christian Rodehacke'
fileout.insitution = 'Danish Meteorological Institute (DMI)'
fileout.address = 'Lyngbyvej 100,  DK-2100 Copenhagen O,  Denmark'
fileout.web = 'http://www.dmi.dk'

if ( do_write_time ) : fileout.date = 'Created on '+time.ctime(time.time())

import os   as os
fileout.working_directory = os.getcwd()
fileout.username = os.environ["USER"]
#?fileout.hostname = os.environ["HOST"]
#fileout.hosttype = os.environ["HOSTTYPE"]
fileout.uname = os.popen('uname -a').read()
#fileout.gridfile1=fileice
#fileout.gridfile2=fileoce

#
# Definition of the dimensions
#
print('    Create dimension')
if ( tlen >= 1 ):
    fileout.createDimension('time',tlen)
else :
    fileout.createDimension('time',tlen)
fileout.createDimension('x1',ilen1)
fileout.createDimension('y1',jlen1)
fileout.createDimension('x2',ilen2)
fileout.createDimension('y2',jlen2)
#fileout.createDimension('z',klen)
fileout.createDimension('one',(1))


#
# Define output variables and specify its attributes
#
print('    Define variables')

## x-, y-arrays
x1_var       = fileout.createVariable('x1', 'd', ('x1', ))
x1_var.standard_name = 'projection_x_coordinate'
x1_var.long_name     = 'X-coordinate in Cartesian system' ;
#x1_var.unit          = 'm'
y1_var       = fileout.createVariable('y1', 'd', ('y1', ))
y1_var.standard_name = 'projection_y_coordinate'
y1_var.long_name     = 'Y-coordinate in Cartesian system' ;
#y1_var.unit          = 'm'

x2_var       = fileout.createVariable('x2', 'd', ('x2', ))
x2_var.standard_name = 'projection_x_coordinate'
x2_var.long_name     = 'X-coordinate in Cartesian system' ;
#x2_var.unit          = 'm'
y2_var       = fileout.createVariable('y2', 'd', ('y2', ))
y2_var.standard_name = 'projection_y_coordinate'
y2_var.long_name     = 'Y-coordinate in Cartesian system' ;
#y2_var.unit          = 'm'


if ( write_latlon ):
    lat1_var       = fileout.createVariable('lat1', 'd', ('x1', 'y1', ))
    lat1_var.standard_name = 'latitude'
    lat1_var.long_name     = 'latitude,  grid 1'
    lat1_var.unit          = 'degree'
    lon1_var       = fileout.createVariable('lon1', 'd', ('x1', 'y1', ))
    lon1_var.standard_name = 'longitude'
    lon1_var.long_name     = 'longitude, grid 1' ;
    lon1_var.unit          = 'degree'
    
    
    lat2_var       = fileout.createVariable('lat2', 'd', ('x2', 'y2', ))
    lat2_var.standard_name = 'latitude'
    lat2_var.long_name     = 'latitude,  grid 2'
    lat2_var.unit          = 'degree'
    lon2_var       = fileout.createVariable('lon2', 'd', ('x2', 'y2', ))
    lon2_var.standard_name = 'longitude'
    lon2_var.long_name     = 'longitude, grid 2' ;
    lon2_var.unit          = 'degree'


isocean_var       = fileout.createVariable('isocean', 'd', ('x2', 'y2', ))
isocean_var.long_name   = 'ocean mask (1:ocean, 0:else)' ;
isocean_var.coordinates ='lat2 lon2'

#
if ( tlen >= 1 ):
    bmr_var       = fileout.createVariable('bmr', 'd', ('time', 'x2', 'y2', ))
else :
    bmr_var       = fileout.createVariable('bmr', 'd', ('x2', 'y2', ))
bmr_var.long_name   ='basal_melting_rate'
bmr_var.units       ='m3/s'
bmr_var.coordinates ='lat2 lon2'
#bmr_var.reference_density=dens_ice

if ( tlen >= 1 ):
    fmr_var       = fileout.createVariable('fmr', 'd', ('time', 'x2', 'y2', ))
else :
    fmr_var       = fileout.createVariable('fmr', 'd', ('x2', 'y2', ))    
fmr_var.long_name   ='frontal_melting_rate'
fmr_var.units       ='m3/s'
fmr_var.coordinates ='lat2 lon2'

#
#
if ( tlen >= 1 ):
    bhr_var       = fileout.createVariable('bhr', 'd', ('time', 'x2', 'y2', ))
else :
    bhr_var       = fileout.createVariable('bhr', 'd', ('x2', 'y2', ))
bhr_var.long_name   ='basal_heating_rate'
bhr_var.units       ='Jm/s'
bhr_var.coordinates ='lat2 lon2'

if ( tlen >= 1 ):
    fhr_var       = fileout.createVariable('fhr', 'd', ('time', 'x2', 'y2', ))
else :
    fhr_var       = fileout.createVariable('fhr', 'd', ('x2', 'y2', ))    
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


filegrid.close()
fileice.close()
fileout.close()

# -----------
#
# Some information
#
print("")
if ( do_write_time ) :
    print(' Done at '+time.ctime(time.time()))
else :
    print(' Done')
