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
import logging
# TODO: optional imports


def latlon_distance(lat1,  lon1,  lat2,  lon2):
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
    deg2meter = 111134.014654 #meter (1 degree latitude equals this distance)
    # rad2deg = 57.2957795131 # rad2deg = 360./(2.*pi)
    deg2rad = 0.0174532925199 # rad2deg = (2.*pi)/360
    dlon0 = lon1 - lon2
    #dlon0 = N.where( N.abs(dlon0)> 180. ,  N.sign(dlon0) * (360.- N.abs(dlon0)) , dlon0 )
    if ( np.abs(dlon0) > 180.):
        if ( dlon0 < -180.):
            dlon0 = -1. * (360. + dlon0)
        else :
            dlon0 = -1. * (360. + dlon0)
    lat0 = (lat1+lat2) * 0.5
    dlon = np.cos( lat0 * deg2rad ) * dlon0 # CAUTION: cos expects rad Not degree
    dlat = (lat1 - lat2)
    # distance
    dist = np.sqrt(dlat*dlat + dlon*dlon) * deg2meter
    return dist





def main():
    ice_file = "pismr_greenland_mask_5km_for_downscaling.nc"
    # TODO: Turn these prints into logging messages
    print("The following steps will be performed...")
    print("0. Load Data files")
    print("...ice data")
    fin_ice = netcdf.netcdf_file(ice_file)
    print fin_ice.variables
    print fin_ice.dimensions
    time_ice = fin_ice.variables["time"].data.squeeze()
    ntime_ice = np.size(time_ice) - 1
    lat_ice = fin_ice.variables["lat"].data.squeeze()
    lon_ice = fin_ice.variables["lon"].data.squeeze()
    x_ice = fin_ice.dimensions["x"]
    y_ice = fin_ice.dimensions["y"]
    print time_ice, ntime_ice, lat_ice, lon_ice, x_ice, y_ice




if __name__ == '__main__':
    main()
