# Makefile for Surface
# Last Modification on 2014.09.22

# F90 Compiler
FC = ifort

# Flags for Maximum Performance
#ifeq ($(mode),debug)
FFLAGS = -debug full -ftrapuv -warn all -traceback -fpe0 -fpp -openmp-stubs -check all
#else
#FFLAGS = -debug full -ftrapuv -warn all -traceback -fpe0 -fpp -op
#FFLAGS = -O2 -openmp-stubs -mkl -fpe0 -fpp
#endif
NFLAGS = -ffixed-line-length-none -ffree-line-length-none

# What archiving to use
AR = ar rcs

# OBJ/SRC List
OBJS =	./SVD/SVD.o	\
	./SURFACE_MODULE/SURFACE_MODULE_2D.o	\
	./SURFACE_MODULE/SURFACE_MODULE_3D.o	\
	./DATA_TRANSFER/DATA_TRANSFER_2D.o	\
	./DATA_TRANSFER/DATA_TRANSFER_3D.o	\
	./HASHING/HASHING_2D.o		\
	./HASHING/HASHING_3D.o		\
	./SURFACES/SURFACES_2D.o	\
	./SURFACES/SURFACES_3D.o	\
	./OPERATORS/OPERATORS_2D.o	\
	./OPERATORS/OPERATORS_3D.o	\
	./HIGHORDER/HIGHORDER_3D.o	\
	./EXTRAPOLATION/EXTRAPOLATION_2D.o	\
	./EXTRAPOLATION/EXTRAPOLATION_3D.o	\
	./PROPA_RECONST_REINITIAL/PROPA_RECONST_REINITIAL_2D.o	\
	./PROPA_RECONST_REINITIAL/PROPA_RECONST_REINITIAL_3D.o	\
	./REMESHING/REMESHING_2D.o	\
	./REMESHING/REMESHING_3D.o	\
	./TOPCHANGE/TOPCHANGE_3D.o	\
	./EXTERNAL_SUBROUTINES/EXTERNAL_SUBROUTINES_2D.o	\
	./EXTERNAL_SUBROUTINES/EXTERNAL_SUBROUTINES_3D.o	\
	./SURFACE_MAIN/SURFACE_MAIN.o	

%.o: %.f90
	$(FC) -c $< -o $@ $(FFLAGS)  

all: $(OBJS)
	$(AR) ./libsurface.a $(OBJS)

clean:
	find ./ -name '*~' -exec rm {} \;
	find ./ -name '*.o' -exec rm {} \;
	find ./ -name '*.mod' -exec rm {} \;
	find ./ -name '*__genmod.f90' -exec rm {} \;

realclean: clean
	rm -f ./libsurface.a
