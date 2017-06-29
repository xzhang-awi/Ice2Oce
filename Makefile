FC=mpiifort
FCFLAGS=-O3 -xHost -DTREAT_OVERLAY -Duse_netCDF -Duse_comm_MPI1 -ip -align -fno-alias
FCFLAGS_MPIOM=-O3 -fpe0 -i4 -fp-model source -fast-transcendentals -no-prec-sqrt -xHost -heap-arrays -convert big_endian -fpp
testfile_raw=/work/ollie/pgierz/mpiesm-pism-sc/experiments/basic/couple/pismr_greenland_mass_and_heat_flux_to_mpiom.nc
testfile_post=FBFLOMPIOM.ext8
make:
	$(FC) $(FCFLAGS) $(FCFLAGS_MPIOM) -o read_glacflux read_glacflux.f90

clean:
	rm read_glacflux FBFLOMPIOM.ext8

test:
	./read_glacflux

generate_test:
	cdo -f ext -b 64B \
		-add \
		-selvar,bmr $(testfile_raw) \
		-selvar,fmr $(testfile_raw) \
		$(testfile_post)

all: clean make generate_test test
