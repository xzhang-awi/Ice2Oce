#!/usr/bin/env python
"""This program takes an ice grid (file1) and generates a map of (i, j) coordinate
pairs of the ice grid onto the corresponding points of the ocean grid (I, J).

ice_file
ocean_file
output_file

optional things:
ice_grid_orientation_as_ocean -- True  # PG: This will be false for Greenland, and ?? for nhem as well as shem
parallelize_calculation -- True
restrict_coast -- True


Paul Gierz, June 2017

"""

import numpy as np
from scipy.io import netcdf
import scipy.ndimage
import logging

# TODO: optional imports


def latlon_distance(lat1, lon1, lat2, lon2):
    """
    calculate the distance between points given as pair of lat/lon values

    Input
    lat1  : latitude point 1 (degree)
    lon1  : longitude point 1 (degree , -180..180)
    lat2  : latitude point 2 (degree)
    lon2  : longitude point 2 (degree , -180..180)

    Output
    dist  : distance between the points (meter)
    """

    # pi = 3.1415926535897932384626433832795028841971693993751058 # https://www.wolframalpha.com 2013-01-21
    # earthradius = 6.36751E6 #earth radius meters # https://www.wolframalpha.com 2013-01-21
    # deg2meter = 2.*pi*earthradius/360.
    deg2meter = 111134.014654  #meter (1 degree latitude equals this distance)
    # rad2deg = 57.2957795131 # rad2deg = 360./(2.*pi)
    deg2rad = 0.0174532925199  # rad2deg = (2.*pi)/360
    dlon0 = lon1 - lon2
    #dlon0 = N.where( N.abs(dlon0)> 180. ,  N.sign(dlon0) * (360.- N.abs(dlon0)) , dlon0 )
    if (np.abs(dlon0) > 180.):
        if (dlon0 < -180.):
            dlon0 = -1. * (360. + dlon0)
        else:
            dlon0 = -1. * (360. + dlon0)
    lat0 = (lat1 + lat2) * 0.5
    dlon = np.cos(
        lat0 * deg2rad) * dlon0  # CAUTION: cos expects rad Not degree
    dlat = (lat1 - lat2)
    # distance
    dist = np.sqrt(dlat * dlat + dlon * dlon) * deg2meter
    return dist


def read_ice(filename):
    with netcdf.netcdf_file(filename, mode="r") as fin_ice:
        time_ice = fin_ice.variables["time"][:].copy().squeeze()
        ntime_ice = np.size(time_ice) - 1
        lat_ice = fin_ice.variables["lat"][:].copy().squeeze()
        lon_ice = fin_ice.variables["lon"][:].copy().squeeze()
        x_ice = fin_ice.variables["x"][:].copy().squeeze()
        y_ice = fin_ice.variables["y"][:].copy().squeeze()
    return time_ice, ntime_ice, lat_ice, lon_ice, x_ice, y_ice


def read_ocean(filename):
    with netcdf.netcdf_file(filename, mode="r") as fin_oce:
        tmask = fin_oce.variables["oces.msk"][:].copy()
        lat2 = fin_oce.variables["oces.lat"][:].copy()
        lon2 = fin_oce.variables["oces.lon"][:].copy()
    return tmask, lat2, lon2


def iscoast_explicit(ilen2, jlen2, isocean):
    iscoast = np.zeros(np.shape(isocean), dtype="bool")
    do_restrict2coast = True
    if do_restrict2coast:
        #
        # Corner points
        i2 = 0
        j2 = 0
        if (isocean[i2, j2] and
                not (isocean[i2 + 1, j2] and isocean[i2, j2 + 1])):
            iscoast[i2, j2] = True
        i2 = ilen2 - 1
        j2 = 0
        if (isocean[i2, j2] and
                not (isocean[i2 - 1, j2] and isocean[i2, j2 + 1])):
            iscoast[i2, j2] = True
        i2 = 0
        j2 = jlen2 - 1
        if (isocean[i2, j2] and
                not (isocean[i2 + 1, j2] and isocean[i2, j2 - 1])):
            iscoast[i2, j2] = True
        i2 = ilen2 - 1
        j2 = jlen2 - 1
        if (isocean[i2, j2] and
                not (isocean[i2 - 1, j2] and isocean[i2 - 1, j2])):
            iscoast[i2, j2] = True
        #
        # Edges along the i2=0 and i2=ilen2-1
        for i2 in range(1, ilen2 - 1):
            j2 = 0
            if (isocean[i2, j2] and
                    not (isocean[i2 - 1, j2] and isocean[i2 + 1, j2] and
                         isocean[i2, j2 + 1] and isocean[i2 - 1, j2 + 1] and
                         isocean[i2 + 1, j2 + 1])):
                iscoast[i2, j2] = True
            j2 = jlen2 - 1
            if (isocean[i2, j2] and
                    not (isocean[i2 - 1, j2] and isocean[i2 + 1, j2] and
                         isocean[i2, j2 - 1] and isocean[i2 - 1, j2 - 1] and
                         isocean[i2 + 1, j2 - 1])):
                iscoast[i2, j2] = True
        #
        # Edges along the j2=0 and j2=jlen2-1
        for j2 in range(1, jlen2 - 1):
            i2 = 0
            if (isocean[i2, j2] and
                    not (isocean[i2 + 1, j2] and isocean[i2, j2 - 1] and
                         isocean[i2, j2 + 1] and isocean[i2 + 1, j2 - 1] and
                         isocean[i2 + 1, j2 + 1])):
                iscoast[i2, j2] = True
                i2 = ilen2 - 1
            if (isocean[i2, j2] and
                    not (isocean[i2 - 1, j2] and isocean[i2, j2 - 1] and
                         isocean[i2, j2 + 1] and isocean[i2 - 1, j2 - 1] and
                         isocean[i2 - 1, j2 + 1])):
                iscoast[i2, j2] = True
        #
        # All inner points
        for i2 in range(1, ilen2 - 1):
            for j2 in range(1, jlen2 - 1):
                if (isocean[i2, j2] and not (
                        isocean[i2 - 1, j2] and isocean[i2 + 1, j2] and
                        isocean[i2, j2 - 1] and isocean[i2, j2 + 1] and
                        isocean[i2 - 1, j2 - 1] and isocean[i2 - 1, j2 + 1] and
                        isocean[i2 + 1, j2 - 1] and isocean[i2 + 1, j2 + 1])):
                    iscoast[i2, j2] = True
    else:
        iscoast = isocean
    return iscoast


def iscoast_binary_erosion(isocean):
    struct = scipy.ndimage.generate_binary_structure(2, 2)
    erode = scipy.ndimage.binary_erosion(isocean, struct, border_value=1)
    iscoast_quick = isocean ^ erode
    return iscoast_quick


def calculate_distances(ilen1, jlen1, lat_ice, lon_ice, lat2, lon2, iscoast):
    ii = jj = 0
    iall = ilen1*jlen1
    fall = iall*1.0
    ic = 0
    import progressbar
    bar = progressbar.ProgressBar(widgets=[
        ' [', progressbar.Timer(), '] ',
        progressbar.Bar(),
        ' (', progressbar.ETA(), ') ',])
    for i1 in bar(range(ilen1)):
        for j1 in range(jlen1):
            #print(" Grid #1: point \t %s of %s (%i\%)" % (ic, iall, ic*100./iall))
            distance = 9.9E24
            Flag = False
            for i2 in range(ilen2):
                for j2 in range(jlen2):
                    if iscoast[i2, j2]:
                        dist = latlon_distance(lat_ice[i1, j1], lon_ice[i1, j1], lat2[i2, j2], lon2[i2, j2])
                        if dist < distance:
                            ii, jj, distance, Flag = i2, j2, dist, True
            if Flag:
                Distance[i1, j1] = distance
                I_1to2[i1, j1] = ii
                J_1to2[i1, j1] = jj
            ic += 1
    return Distance, I_1to2, J_1to2




#def main():
if __name__ == '__main__':
    ice_file = "test_data/pism_ini.nc"
    # TODO: Turn these prints into logging messages
    print("The following steps will be performed...")
    print("Load Data files")
    print("...ice data")
    time_ice, ntime_ice, lat_ice, lon_ice, x_ice, y_ice = read_ice(
        "test_data/pism_ini.nc")
    print("...ocean_data")
    tmask, lat2, lon2 = read_ocean(
        "/home/ollie/pgierz/reference_stuff/GR15_lsm_lon_lat.nc")
    print("...finished!")
    print("Deriving some values...")
    lon_ice, lon2 = np.where(lon_ice>180.0, lon_ice-360.0, lon_ice), np.where(lon2>180.0, lon2-360.0, lon2)
    ilen1, jlen1 = np.size(x_ice), np.size(y_ice)
    ilen2, jlen2 = np.shape(lat2)
    size2D_ice = (ilen1, jlen1)
    size2D_ocean = (ilen2, jlen2)
    isocean = np.where(tmask > 0, True, False)
    iscoast = iscoast_explicit(ilen2, jlen2, isocean)
    iscoast_quick = iscoast_binary_erosion(isocean)
    I_1to2 = J_1to2 = np.ones(size2D_ice, dtype="i") * -9
    Distance = np.ones(size2D_ice, dtype="f") * -9.99999
    print("...finished!")
    print("Main loop...")
    Distance, I_1to2, J_1to2 = calculate_distances(ilen1, jlen1, lat_ice, lon_ice, lat2, lon2, iscoast)
    print("...finished!")
    import matplotlib.pyplot as plt
    f, (ax1, ax2, ax3) = plt.subplots(1, 3)
    ax1.pcolormesh(Distance, cmap="jet")
    ax2.pcolormesh(I_1to2, cmap="jet")
    ax3.pcolormesh(J_1to2, cmap="jet")
    plt.show()

# if __name__ == '__main__':
#     main()
