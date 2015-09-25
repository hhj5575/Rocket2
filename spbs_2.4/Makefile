# Makefile for spbs
# by Changsoo Lee

# execution file name
EXE = spbs

# src / obj / lib / bin dir
SRC_DIR     =   src
OBJ_DIR     =   obj
LIB_DIR     =   lib
BIN_DIR     =   bin

# f90 compiler
FC = ifort

# flags for maximum performance
ifeq ($(mode),debug)
FFLAGS = -O0 -debug full -ftrapuv -warn all -traceback -check all -fpe0 -fpp -module $(OBJ_DIR)
else
FFLAGS = -O3 -fpe0 -fpp -module $(OBJ_DIR)
endif

# include needed for compile
IFLAGS =

# library needed for linking
LFLAGS =

# whar archiving to use
AR = xiar rcs

# target obj list
OBJS    =   $(OBJ_DIR)/TYPE.o\
            $(OBJ_DIR)/REGION.o\
            $(OBJ_DIR)/INTERFACE.o\
            $(OBJ_DIR)/BURN.o\
            $(OBJ_DIR)/INITIAL.o\
            $(OBJ_DIR)/SOLVER.o\
            $(OBJ_DIR)/POST.o\
            $(OBJ_DIR)/REMESH.o\
            $(OBJ_DIR)/MAIN.o
EXEC    =   $(OBJ_DIR)/EXEC.o

.PHONY: all clean distclean

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) -c $< -o $@ $(FFLAGS) $(IFLAGS)

all: $(OBJS) $(EXEC)
	$(FC) -o $(BIN_DIR)/$(EXE) $(OBJS) $(EXEC) $(FFLAGS) $(LFLAGS)
	$(AR) $(LIB_DIR)/lib$(EXE).a $(OBJS)

clean:
	find $(OBJ_DIR) -name '*~' -exec rm {} \;
	find $(OBJ_DIR) -name '*.o' -exec rm {} \;
	find $(OBJ_DIR) -name '*.mod' -exec rm {} \;
	find $(OBJ_DIR) -name '*__genmod.f90' -exec rm {} \;

distclean: clean
	rm -f $(BIN_DIR)/$(EXE)
	rm -f $(LIB_DIR)/lib$(EXE).a
