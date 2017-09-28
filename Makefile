FC      = /home/plaplant/local/hdf5-1.10.1/bin/h5fc
FFLAGS  = -O3 -openmp -g -traceback -mcmodel=medium -shared-intel
DBFLAGS = -check -warn -debug

FGSL_INC = `pkg-config --cflags fgsl`
FGSL_LIB = `pkg-config --libs fgsl`

HEALPIX     = $(HOME)/local/Healpix_3.30
HEALPIX_INC = -I$(HEALPIX)/includef90
HEALPIX_LIB = -L$(HEALPIX)/libf90 -lhealpix

CFITSIO     = $(HOME)/local/cfitsio-3.370
CFITSIO_INC = -I$(CFITSIO)/include
CFITSIO_LIB = -L$(CFITSIO)/lib -lcfitsio

HDF5 = $(HOME)/local/hdf5-1.10.1
HDF5_INC = -I$(HDF5)/include
HDF5_LIB = -L$(HDF5)/lib -lhdf5 -lhdf5_hl

INC = $(HEALPIX_INC) $(CFITSIO_INC) $(HDF5_INC) $(FGSL_INC)
LIB = $(HEALPIX_LIB) $(CFITSIO_LIB) $(HDF5_LIB) $(FGSL_LIB)

OBJ = xi_global.o io_tools.o xi_tools.o

ifeq ($(DEBUG),1)
FFLAGS += $(DBFLAGS)
else
FFLAGS += -w
endif

# executable
calculate_xi.x: $(OBJ)
	$(FC) $(FFLAGS) $(INC) $(OBJ) calculate_xi.f90 -o $@ $(LIB)

tests.x: $(OBJ)
	$(FC) $(FFLAGS) $(INC) $(OBJ) tests.f90 -o $@ $(LIB)

%.o: %.f90
	$(FC) $(FFLAGS) $(INC) -c $*.f90

.PHONY: all
all: calculate_xi.x

.PHONY: test
test: tests.x

.PHONY: clean
clean:
	rm -rf *.o *.mod *.x
