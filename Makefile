# Compiler Settings
FC=mpif90
FCFLAGS=-O3 -xHost -DTREAT_OVERLAY -Duse_netCDF -Duse_comm_MPI1 -ip -align -fno-alias \
		-L/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-intel14/lib \
		-lnetcdff \
		-I/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-intel14/include \
    -Wl,-rpath,/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-intel14/lib

FCFLAGS_MPIOM=-O3 -fpe0 -i4 -fp-model source -fast-transcendentals -no-prec-sqrt -xHost -heap-arrays -convert big_endian -fpp

# Path to test files
testfile_raw=/pf/a/a270077/freshwater_tests/hosing_mask_IRD_1Sv_std.nc
testfile_post=FBFLOMPIOM.ext8
make:
	$(FC) $(FCFLAGS) $(FCFLAGS_MPIOM) -o read_glacflux read_glacflux.f90

clean:
	rm read_glacflux FBFLOMPIOM.ext8 FBFLOMPIOM.txt

test:
	./read_glacflux

generate_test:
	cdo -f ext -b 64B \
		-selvar,area $(testfile_raw) \
		$(testfile_post)

	cdo -outputf,%12.8e,1 \
		$(testfile_post) > FBFLOMPIOM.txt

all: clean make generate_test test
