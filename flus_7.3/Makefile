# Makefile for FLUS
# by Changsoo Lee

# execution file name
EXE = flus

# src / obj / ext / lib / bin dir
INCLUDE_DIR =   include
SRC_DIR     =   src
OBJ_DIR     =   obj
EXT_DIR     =   ext
LIB_DIR     =   lib
BIN_DIR     =   bin

# f90 compiler for intel MPI
FC = mpif90

# flags for maximum performance
ifeq ($(mode),debug)
FFLAGS = -O0 -debug full -ftrapuv -warn all -traceback -check all -fpe0 -fpp -fp-model precise -module $(OBJ_DIR)
else
FFLAGS = -O3 -fpe0 -fpp -module $(OBJ_DIR)
endif

# include needed for compile
IFLAGS = -I $(INCLUDE_DIR)

# library needed for linking
LFLAGS = -L $(LIB_DIR) -ltri -ltet -lmetis -lstdc++

# define needed for preprocessor
DFLAGS = -D $(def)

# what archiving to use
AR = xiar rcs

# object list
OBJS  = $(OBJ_DIR)/config_module.o\
        $(OBJ_DIR)/grid_module.o\
        $(OBJ_DIR)/numerics_module.o\
        $(OBJ_DIR)/quadtree_module.o\
        $(OBJ_DIR)/octree_module.o\
        $(OBJ_DIR)/data_module.o\
        $(OBJ_DIR)/mixture_module.o\
        $(OBJ_DIR)/time_module.o\
        $(OBJ_DIR)/post_module.o\
        $(OBJ_DIR)/interface_module.o\
        $(OBJ_DIR)/process_module.o\
        $(OBJ_DIR)/flus_module.o
EXEC  = $(OBJ_DIR)/flus_exec.o

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) -c $< -o $@ $(FFLAGS) $(IFLAGS) $(DFLAGS)

.PHONY: all build externals clean distclean

all: externals build

build:	$(OBJS) $(EXEC)
	$(FC) -o $(BIN_DIR)/$(EXE) $(OBJS) $(EXEC) $(FFLAGS) $(LFLAGS)
	$(AR) $(LIB_DIR)/lib$(EXE).a $(OBJS)
	cp $(OBJ_DIR)/flus_module.mod $(INCLUDE_DIR)/

externals:
	@$(MAKE) -C $(EXT_DIR)/trigen
	@$(MAKE) -C $(EXT_DIR)/tetgen tetlib
	@$(MAKE) -C $(EXT_DIR)/converter
	cp $(EXT_DIR)/trigen/libtri.a $(LIB_DIR)
	cp $(EXT_DIR)/tetgen/libtet.a $(LIB_DIR)
	cp $(EXT_DIR)/converter/converter $(BIN_DIR)

clean:
	find $(OBJ_DIR) -name '*~' -exec rm {} \;
	find $(OBJ_DIR) -name '*.o' -exec rm {} \;
	find $(OBJ_DIR) -name '*.mod' -exec rm {} \;
	find $(OBJ_DIR) -name '*__genmod.f90' -exec rm {} \;

distclean: clean
	@$(MAKE) clean -C $(EXT_DIR)/trigen
	@$(MAKE) clean -C $(EXT_DIR)/tetgen
	@$(MAKE) clean -C $(EXT_DIR)/converter
	rm -f $(INCLUDE_DIR)/flus_module.mod
	rm -f $(LIB_DIR)/libtri.a
	rm -f $(LIB_DIR)/libtet.a
	rm -f $(BIN_DIR)/converter
	rm -f $(BIN_DIR)/$(EXE)
	rm -f $(LIB_DIR)/lib$(EXE).a
