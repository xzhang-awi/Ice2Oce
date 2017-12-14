clc;clear;
cd /Users/xzhang/Dataset/DataOutside/Model_Output/Icesheet_Tarasov/Cal_drainage_point/input_for_MATLAB_to_cal_freshwater_to_ocean_from_point_field/Ice2Oce

%readin files
sfile=ncread('21ka_ice_volume_mask_DHDt_afterPointer_adjLon_GR15_Positive.nc','fwf');

%% convert the row from 254 to 256
tfile=zeros(256,220);
tfile(2:255,:)=sfile(:,:);

%% write to a netcdf file
% make a copy of file
system('cp Refer_256x220_grid.nc 21ka_positive_DhDt_256x220grid.nc');

ncid=netcdf.open('21ka_positive_DhDt_256x220grid.nc','NC_WRITE');

%      i=55
%     [varname,xtype,dimids,natts]=netcdf.inqVar(ncid,i); 
%     i, varname

netcdf.putVar(ncid,7,tfile);
netcdf.close(ncid);

%% plot the freshwater flux with land-sea mask
system('/usr/local/bin/cdo pardes 21ka_positive_DhDt_256x220grid.nc')
system('/usr/local/bin/cdo mulc,-1 -subc,1 LGM_GR15_256x220_OceanMask.nc slm')
system('/usr/local/bin/cdo sub -div 21ka_positive_DhDt_256x220grid.nc 21ka_positive_DhDt_256x220grid.nc slm fwf-slm.nc')





