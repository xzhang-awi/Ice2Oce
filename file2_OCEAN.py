# Compute from a given (smaller ice sheet) Lat- and Lon-fields, fields #1, for
# a (larger) other lat-, lon-fields, field #2, of the corresponding grid point
# indexes.
#
# This routine computes the Grid1to2-File which represents the transformation
#
# This stuff is needed to determine for given ice sheet model grids (PISM) the
# corresponding ocean grid points.


## 4 Tests
##
#fileocean = 'mask.select.nc'
#fileice = 'pism_ini.selected.nc'
#
#fileocean = 'mask.pism_ini.nc'
#fileice = 'pism_ini.nc'

# Greenland -> EC-Earth
#fileocean   = 'ORCA1_MM_1850_grid_T.nc'
#Nemo ORCA1: fileocean   = 'mask.nc'
#IFS T159   fileocean   = 'IFS_grid.nc' # ncrename -v atmo.msk,landmask masks_orig.nc IFS_grid.nc # ncks -A -v atmo.lon,atmo.lat atmo_grid.nc IFS_grid.nc 
fileocean   = 'IFS_grid.nc'
fileice     = 'pism_ini.nc' #file1
fileout_name= 'Grid1to2.nc' #file2

# Antarctic -> EC-Earth
fileocean   = 'IFS_grid.nc'
fileice     = 'AntSynne.nc' #file1
fileout_name= 'Grid1to2.Ant2ATM.nc' #file2


## Antarctic -> EC-Earth
fileocean   = 'mask.nc'
fileice     = 'AntSynne.nc' #file1
fileout_name= 'Grid1to2.AntPISM20km.allOcean.nc' #file2

#
# Switches that control the behavior of the code
#
do_plot         = False # Switch off: False, Default: False
#do_plot         = True
do_write_time   = True  # Switch off: False, Default: True

do_fortran_boundIJ = True       # Index start at one/1 (=Fortran): True; zero/0: False
#Nemo ORCA1: do_reorder_bad_points = True   # Reorder some points (to avoid problems: eg. move points out of fjord to the next coast)
#IFS T159: do_reorder_bad_points = False   # Reorder some points (to avoid problems: eg. move points out of fjord to the next coast)
do_reorder_bad_points = False
#Nemo ORCA1: do_grid1_2_same_orientation = False  # True: I_1to2=i , J_1to2=j; False: I_1to2=j , J_1to2=i
#IFS T159: do_grid1_2_same_orientation = False  # True: I_1to2=i , J_1to2=j; False: I_1to2=j , J_1to2=i
do_grid1_2_same_orientation = False  # false , when PISM grid is tilted..

do_parallel     = False #False=default !!

#Nemo ORCA1: do_restrict2coast = True
#IFS T159: do_restrict2coast = False
do_restrict2coast = True
do_restrict2coast = False

#
#
#
print("")
print("Import packages / modules")
import numpy as N
import Scientific.IO.NetCDF as NC
if ( do_plot ):         import matplotlib.pyplot    as PLT # import matplotlib.pylab    as PLT 
if ( do_write_time ) :  import time as time
if ( do_parallel ) :    import pp as pp
    
print("       Done import packages / modules")

if ( do_write_time ) :  print('Start at '+time.ctime(time.time()))


# Reading the flags from the command line

def latlon_distance(lat1,  lon1,  lat2,  lon2):
    # calculate the distance between points given as pair of lat/lon values
    # Input
    # lat1  : latitude point 1 (degree)
    # lon1  : longitude point 1 (degree , -180..180)
    # lat2  : latitude point 2 (degree)
    # lon2  : longitude point 2 (degree , -180..180)
    # Output
    # dist  : distance between the points (meter)


    #pi = 3.1415926535897932384626433832795028841971693993751058 # https://www.wolframalpha.com 2013-01-21
    #earthradius = 6.36751E6 #earth radius meters # https://www.wolframalpha.com 2013-01-21
    #deg2meter = 2.*pi*earthradius/360.
    deg2meter = 111134.014654 #meter (1 degree latitude equals this distance)
    
    #rad2deg = 57.2957795131 # rad2deg = 360./(2.*pi)
    deg2rad = 0.0174532925199 # rad2deg = (2.*pi)/360
    
    dlon0 = lon1 - lon2
    #dlon0 = N.where( N.abs(dlon0)> 180. ,  N.sign(dlon0) * (360.- N.abs(dlon0)) , dlon0 )
    if ( N.abs(dlon0) > 180.) :
        if ( dlon0 < -180. ):
            dlon0 = -1. * (360. + dlon0)
        else :
            dlon0 = -1. * (360. + dlon0)                
     
    lat0 = (lat1+lat2) * 0.5
    dlon = N.cos( lat0 * deg2rad ) * dlon0 # CAUTION: cos expects rad Not degree
    dlat = (lat1 - lat2)
    # distance
    dist = N.sqrt(dlat*dlat + dlon*dlon) * deg2meter
    # phase angle
    #angle = 
    #% CALCUALTE ANGLE TO X AXIS % angle buildin function of octave/matlab
    #phaseangle  = angle(dep+dlat*sqrt(-1))*rad2deg
    return dist

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
print("  do_restrict2coast      : "+str(do_restrict2coast))

print("  do_reorder_bad_points  : "+str(do_reorder_bad_points))
print("  do_grid1_2_same_orientation: "+str(do_grid1_2_same_orientation))
#if ( do_grid1_2_same_orientation ) :
#    print(" \t Same orientation I=i, J=j")
#else :
#    print(" \t Transposed orientation I=j, J=i")
    
print("  do_parallel            : "+str(do_parallel))
#ppservers = ("*", )
ppservers = ()
if (do_parallel):
    print " Do have started outside of the python script 'ppserver -s lala007'"
    job_server = pp.Server(ppservers=ppservers, secret="lala007")
    #job_server = pp.Server(ncpus=2)
    print(" \t Starting pp with", job_server.get_ncpus(), "workers")



print("weiter")
#
#

print("")

# -------------- 
#
# Open the files and access the data
#
print("Load data")

# -----------
#
# ice sheet/shelf data
#
# FIXME: Assumption about the PISM/ice sheet data file: time is the first index
print(" Ice model data file  : "+fileice)
fileice  = NC.NetCDFFile(fileice,mode="r")
time_ice = fileice.variables['time'].getValue()
time_ice0= N.size(time_ice) - 1
lat1     = fileice.variables['lat'].getValue()
lon1     = fileice.variables['lon'].getValue()
xpism    = fileice.variables['x'].getValue()
ypism    = fileice.variables['y'].getValue()

# -----------
#
# ocean data
#
print(" Ocean model data file: "+fileocean)
fileoce  = NC.NetCDFFile(fileocean,mode="r")
#time_oce = fileoce.variables['time_counter'].getValue()
# FIXME: Assumption about the NEMO/ocean data file
#   1. time is first index in the 
#   2. lat/lon are the second and third indexes
#Nemo ORCA1
#Nemo ORCA1: #sss      = fileoce.variables['sosaline'].getValue()[0, :, :]
#Nemo ORCA1: tmask    = fileoce.variables['tmask'].getValue()[0,0, :, :]
#Nemo ORCA1: lat2     = fileoce.variables['nav_lat'].getValue()
#Nemo ORCA1: lon2     = fileoce.variables['nav_lon'].getValue()

tmask    = fileoce.variables['tmask'].getValue()[0,0, :, :]
lat2     = fileoce.variables['nav_lat'].getValue()
lon2     = fileoce.variables['nav_lon'].getValue()

#IFS T159
#IFS T159: tmask    = fileoce.variables['landmask'].getValue() # 1=land
#IFS T159: lat2     = fileoce.variables['atmo.lat'].getValue()
#IFS T159: lon2     = fileoce.variables['atmo.lon'].getValue()


#=============================================
#
#
# Deduced variables
#
ilen1 = N.size(xpism)   # ilen = N.size(temp[0,0,:,0])   
jlen1 = N.size(ypism)   # jlen = N.size(temp[0,0,0,:])

[ilen2, jlen2] = N.shape(lat2)

#tlen = N.size(time_oce) # N.size(temp[:, 1, 1, 1])
#
size2Dice   = (ilen1,  jlen1)
size2Doce   = (ilen2,  jlen2)
# --
# From file #2 (ocean file)
#
#isocean     = N.where(sss > 0 , True, False)
isocean     = N.where(tmask > 0 , True, False)

    #print " For test: isocean = True" ; isocean = N.ones(N.shape(isocean), dtype='bool')

#
# Determine coast point in the ocean
#
iscoast = N.zeros(N.shape(isocean), dtype='bool')
if ( do_restrict2coast ) :
    #
    # Corner points
    i2=0 ; j2=0
    if ( isocean[i2, j2] and not (isocen[i2+1, j2] and isocen[i2, j2+1]) ) : iscoast[i2, j2]=True       
    i2=ilen2-1 ; j2=0
    if ( isocean[i2, j2] and not (isocen[i2-1, j2] and isocen[i2, j2+1]) ) : iscoast[i2, j2]=True       
    i2=0 ; j2=jlen2-1
    if ( isocean[i2, j2] and not (isocen[i2+1, j2] and isocen[i2, j2-1]) ) : iscoast[i2, j2]=True       
    i2=ilen2-1 ; j2=jlen2-1
    if ( isocean[i2, j2] and not (isocen[i2-1, j2] and isocen[i2-1, j2]) ) : iscoast[i2, j2]=True       
    
    #
    # Edges along the i2=0 and i2=ilen2-1
    for i2 in range(1, ilen2-1):
        j2=0
        if ( isocean[i2, j2] and not 
            (isocean[i2-1, j2]   and isocean[i2+1, j2]   and  isocean[i2, j2+1] and 
             isocean[i2-1, j2+1] and isocean[i2+1, j2+1]  )) :
            iscoast[i2, j2] = True
        j2=jlen2-1
        if ( isocean[i2, j2] and not 
            (isocean[i2-1, j2]   and isocean[i2+1, j2]   and isocean[i2, j2-1] and  
             isocean[i2-1, j2-1] and isocean[i2+1, j2-1]  )) :
            iscoast[i2, j2] = True
                
                
    #
    # Edges along the j2=0 and j2=jlen2-1           
    for j2 in range(1, jlen2-1):                
        i2=0
        if ( isocean[i2, j2] and not 
            ( isocean[i2+1, j2]   and isocean[i2, j2-1]   and isocean[i2, j2+1] and
              isocean[i2+1, j2-1] and isocean[i2+1, j2+1]  )) :
            iscoast[i2, j2] = True
        i2=ilen2-1
        if ( isocean[i2, j2] and not 
            (isocean[i2-1, j2]   and isocean[i2, j2-1]   and isocean[i2, j2+1]  and  
             isocean[i2-1, j2-1] and isocean[i2-1, j2+1]  )) :
            iscoast[i2, j2] = True
    
    #
    # All inner points
    for i2 in range(1, ilen2-1):    
        for j2 in range(1, jlen2-1):                
            if ( isocean[i2, j2] and not 
                (isocean[i2-1, j2]   and isocean[i2+1, j2]   and isocean[i2, j2-1]   and isocean[i2, j2+1]  and  
                 isocean[i2-1, j2-1] and isocean[i2-1, j2+1] and isocean[i2+1, j2-1] and isocean[i2+1, j2+1]  )) :
                    iscoast[i2, j2] = True
else :
    iscoast = isocean

#
# Longitudes in the range [-180 to 180 deg]
# Latitudes in the range [-90 to 90 deg]
#
#pi = 3.14159 ; Rad2Dec = 180./pi
#lat1 = lat1*Rad2Dec
#lon1 = lon1*Rad2Dec
lon1 = N.where(lon1>180.0,  lon1-360.0, lon1 )
lon2 = N.where(lon2>180.0,  lon2-360.0, lon2 )


I_1to2 = N.ones(size2Dice,  dtype='i') * -9 #I_1to2 = N.zeros(size2Dice,  dtype='i')
J_1to2 = N.ones(size2Dice,  dtype='i') * -9 #J_1to2 = N.zeros(size2Dice,  dtype='i')
Distance= N.ones(size2Dice) * -9.99999
#Angle   = N.zeros(size2Dice)

ii  = 0;    jj  = 0
iall= ilen1*jlen1;  fall= iall*1.0
print(" Main loop")
ic  = 0
# Loop i1,j1 over all ice grid points (grid #1)
for i1 in range(ilen1):
    for j1 in range(jlen1): 
        print("   Grid #1 : point \t"+str(ic)+" of "+str(iall)+"\t\t("+str(ic*100./fall)+'%)')
        distance = 9.9E24
        Flag    = False
        # For a given point (i1,j2) Loop i2,j2 over all ocean grid points (grid #2)
        for i2 in range(ilen2):
            for j2 in range(jlen2):                
                #if ( isocean[i2, j2] ): # Only ocean point should be considered
                if ( iscoast[i2, j2] ): # Only costal ocean point should be considered
                    dist = latlon_distance(lat1[i1, j1],  lon1[i1, j1],  lat2[i2, j2],  lon2[i2, j2])
                    if ( dist < distance ) :
                        ii=i2
                        jj=j2
                        distance = dist
                        #angle = angl
                        Flag = True

        if ( Flag ) :
            Distance[i1, j1] = distance        
            #Angle[i1, j1]   = angle
            I_1to2[i1, j1] = ii
            J_1to2[i1, j1] = jj
        ic += 1

#
# Reorder points?
#
if ( do_reorder_bad_points ) :
    print(" I might reorder certain points?")
    # All these vectors have to have the same length!
    # Enter the indexces in grid #2
    #
    # Caution
    # 1) C-style, zero-bound: index start at zero
    # 2) I/J index might be switched depending on the orientation of
    #    both input fields. Check results.

    Ibad    = [244, 244, 244]
    Jbad    = [255, 256, 257]
    Igood   = [242, 243, 242]
    Jgood   = [255, 256, 257]
    
    Jbad    = [244, 244]
    Ibad    = [255, 257]
    Jgood   = [242, 243]
    Igood   = [255, 257]    
    if ( N.size(Ibad) == N.size(Jbad) and N.size(Igood) == N.size(Jgood) and N.size(Igood) == N.size(Ibad) ):
        ijlen = N.size(Igood)
        print(" \t I will move "+str(ijlen)+" data points")
        print(" \t\t Size I_1to2"+str(N.size(I_1to2)))
        for ic in range(ijlen):
            print(" \t\t (\t"+str(Ibad[ic])+" /\t"+str(Jbad[ic])+" )\t -> (\t"+str(Igood[ic])+" /\t"+str(Jgood[ic])+" )")
            for i1 in range(ilen1):
                for j1 in range(jlen1):
                    if ( I_1to2[i1, j1] == Ibad[ic] and J_1to2[i1, j1] == Jbad[ic] ) :
                        I_1to2[i1, j1] = Igood[ic]
                        J_1to2[i1, j1] = Jgood[ic]
    else :
        print("** Size of the input data arrays do not fit. Please correct")
        print("\t Size Ibad  "+str(N.size(Ibad)))
        print("\t Size Jbad  "+str(N.size(Jbad)))
        print("\t Size Igood "+str(N.size(Igood)))
        print("\t Size Jgood "+str(N.size(Jgood)))
        print("S T O P")
        exit()
else :
    print(" Skip any reorder, since 'do_reorder_bad_points'="+str(do_reorder_bad_points))

#
# Switch from zero-bound index scheme (C, python) to one-bound index scheme (fortran) ?
#
if ( do_fortran_boundIJ ) :
    print(" Index start at one/1 (Fortran-bound)")
    I_1to2 = I_1to2 + 1
    J_1to2 = J_1to2 + 1
else :
    print(" Index start at zero/0 (C-bound, python-bound)")

print("")

# -------------- Output file ----------------
#
print("Save data")

# -----------
#
# ice sheet/shelf data
#
print(" Ice-ocean interaction data file  : "+fileout_name)
fileout = NC.NetCDFFile(fileout_name,mode='w')  # Open the output file

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
#fileout.createDimension('time',tlen)
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

# x-, y-arrays
x1_var       = fileout.createVariable('x1', 'd', ('x1', ))
x1_var.standard_name = 'projection_x_coordinate'
x1_var.long_name     = 'X-coordinate in Cartesian system' ;
x1_var.unit          = 'm'
y1_var       = fileout.createVariable('y1', 'd', ('y1', ))
y1_var.standard_name = 'projection_y_coordinate'
y1_var.long_name     = 'Y-coordinate in Cartesian system' ;
y1_var.unit          = 'm'

x2_var       = fileout.createVariable('x2', 'd', ('x2', ))
x2_var.standard_name = 'projection_x_coordinate'
x2_var.long_name     = 'X-coordinate in Cartesian system' ;
#x2_var.unit          = 'm'
y2_var       = fileout.createVariable('y2', 'd', ('y2', ))
y2_var.standard_name = 'projection_y_coordinate'
y2_var.long_name     = 'Y-coordinate in Cartesian system' ;
#y2_var.unit          = 'm'



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

iscoast_var       = fileout.createVariable('iscoast', 'd', ('x2', 'y2', ))
iscoast_var.long_name   = 'coast mask (1:ocean coast, 0:else)' ;

#

I_1to2_var       = fileout.createVariable('I_1to2', 'd', ('x1', 'y1', ))
if ( do_fortran_boundIJ ) :
    I_1to2_var.long_name     = 'i-coordinate of grid #1 in grid #2, I=i (start at one/1, Fortran-bound)'
else :
    I_1to2_var.long_name     = 'i-coordinate of grid #1 in grid #2, I=i (start at zero/0, C-bound)'    

J_1to2_var       = fileout.createVariable('J_1to2', 'd', ('x1', 'y1', ))
if ( do_fortran_boundIJ ) :
    J_1to2_var.long_name     = 'j-coordinate of grid #1 in grid #2, J=j (start at one/1, Fortran-bound)'
else :
    J_1to2_var.long_name     = 'j-coordinate of grid #1 in grid #2, J=j (start at zero/0, C-bound)'    

TI_1to2_var      = fileout.createVariable('TI_1to2', 'd', ('x1', 'y1', ))
if ( do_fortran_boundIJ ) :
    TI_1to2_var.long_name    = 'Switched i-coordinate of grid #1 in grid #2, I=j (start at one/1, Fortran-bound)'
else :
    TI_1to2_var.long_name    = 'Switched i-coordinate of grid #1 in grid #2, I=j (start at zero/0, C-bound)'    

TJ_1to2_var      = fileout.createVariable('TJ_1to2', 'd', ('x1', 'y1', ))
if ( do_fortran_boundIJ ) :
    TJ_1to2_var.long_name    = 'Switched j-coordinate of grid #1 in grid #2, J=i (start at one/1, Fortran-bound)'
else :
    TJ_1to2_var.long_name    = 'Switched j-coordinate of grid #1 in grid #2, J=i (start at zero/0, C-bound)'    


dist_var       = fileout.createVariable('distance', 'd', ('x1', 'y1', ))
dist_var.long_name     = 'distance to next vaild point (grid#1 -> grid#2)' ;
dist_var.unit          = 'm'

#angle_var       = fileout.createVariable('angle', 'd', ('x1', 'y1', ))
#angle_var.long_name     = 'angle / direction to next vaild point' ;
#angle_var.unit          = 'degree'


#
#
#
print('    Populate variables')

# x-, y-arrarys
x1_var[:]   = xpism;        y1_var[:]   = ypism;
x2_var[:]   = range(ilen2); y2_var[:]   = range(jlen2)

lat1_var[:] = lat1
lon1_var[:] = lon1

lat2_var[:] = lat2
lon2_var[:] = lon2

isocean_var[:]=isocean
iscoast_var[:]=iscoast

#I_1to2_var[:] = I_1to2
#J_1to2_var[:] = J_1to2
#TI_1to2_var[:]= J_1to2
#TJ_1to2_var[:]= I_1to2

if ( do_grid1_2_same_orientation ) :
    I_1to2_var[:] = I_1to2
    J_1to2_var[:] = J_1to2
    TI_1to2_var[:]= J_1to2
    TJ_1to2_var[:]= I_1to2
else :
    I_1to2_var[:] = J_1to2
    J_1to2_var[:] = I_1to2
    TI_1to2_var[:]= I_1to2
    TJ_1to2_var[:]= J_1to2

dist_var[:] = Distance



################
#
# Close the files
#
print("Close all open files")

fileoce.close()
fileice.close()
fileout.close()


################
#
# -------------- ploting ----------------
#
if ( do_plot ):
    print("Figures / plots")
    # ---    
    PLT.figure()    
    mypcolor=PLT.pcolor(isocean)
    PLT.colorbar(mypcolor, orientation='vertical')

    mycontour1 = PLT.contour(lon2, range(-180,180,15),colors='k')
    PLT.clabel(mycontour1,fontsize=12,fmt="%3.0f")
    mycontour2 = PLT.contour(lat2, range( -90, 90,15),colors='w')
    PLT.clabel(mycontour2,fontsize=12,fmt="%3.0f")
    PLT.title("ocean (ocean: 1, else:0), black lines: lon, white lines: lat")
    # ---
    
    # ---    
    PLT.figure()    
    mypcolor=PLT.pcolor(iscoast)
    PLT.colorbar(mypcolor, orientation='vertical')

    mycontour1 = PLT.contour(lon2, range(-180,180,15),colors='k')
    PLT.clabel(mycontour1,fontsize=12,fmt="%3.0f")
    mycontour2 = PLT.contour(lat2, range( -90, 90,15),colors='w')
    PLT.clabel(mycontour2,fontsize=12,fmt="%3.0f")
    PLT.title("oceanic coast (coast: 1, else:0), black lines: lon, white lines: lat")
    # ---    
    # ---    
    PLT.figure()    
    
    PLT.subplot(121)
    mypcolor=PLT.pcolor(I_1to2)
    PLT.colorbar(mypcolor, orientation='vertical')
    PLT.title("I_1to2")

    PLT.subplot(122)
    mypcolor=PLT.pcolor(J_1to2)
    PLT.colorbar(mypcolor, orientation='vertical')
    PLT.title("J_1to2")
    
    #PLT.title("ocean")
    # ---
    
    PLT.draw()
    

    print("  continue")
else :
    print(" NO Figures / plots")
    print(" ")
    
    
################

# -----------
#
# Some information
#
print("")
if ( do_write_time ) :
    print(' Done at '+time.ctime(time.time()))
else :
    print(' Done')

if ( do_plot ):
    print("Close figure window(s) to complete the program... ")
    PLT.show()
    
print("Bye bye ..")
