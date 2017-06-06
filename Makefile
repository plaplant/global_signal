FC      = ifort
FFLAGS  = -O3 -openmp -g -traceback -mcmodel=medium -shared-intel
DBFLAGS = -check -warn -debug

FGSL_INC = `pkg-config --cflags fgsl`
FGSL_LIB = `pkg-config --libs fgsl`

HEALPIX     = $(HOME)/local/Healpix_3.30
HEALPIX_INC = -I$(HEALPIX)/include
HEALPIX_LIB = -L$(HEALPIX)/lib -lhealpix

CFITSIO     = $(HOME)/local/cfitsio-3.390
CFITSIO_INC = -I$(CFITSIO)/include
CFITSIO_LIB = -L$(CFITSIO)/lib -lcfitsio

INC = $(HEALPIX_INC) $(CFITSIO_INC) $(FGSL_INC)
LIB = $(HEALPIX_LIB) $(CFITSIO_LIB) $(FGSL_LIB)

OBJ = xi_global.o io_tools.o xi_tools.o

ifeq ($(DEBUG),1)
FFLAGS += $(DBFLAGS)
else
FFLAGS += -w
endif

# executable
calculate_xi.x: $(OBJ)
	$(FC) $(FFLAGS) $(INC) $(OBJ) calculate_xi.f90 $(LIB) -o $@

%.o: %.f90
	$(FC) $(FFLAGS) $(INC) -c $*.f90

.PHONY: all
all: calculate_xi.x

.PHONY: clean
clean:
	rm -rf *.o *.mod *.x
