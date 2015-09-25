# Makefile for TRIGEN
# by Changsoo Lee
# Last Modification on 2014.04.23

# F90 Compiler
FC = ifort

# Flags for Maximum Performance
FFLAGS = -O3

# What archiving to use
AR = ar rcs

# OBJ/SRC List
OBJS =	./src/header.o		\
	./src/a_trigen.o	\
	./src/b_allocvars.o	\
	./src/c_input.o		\
	./src/d_initbd.o	\
	./src/e_connct.o	\
	./src/f_findpt.o	\
	./src/g_degen2.o	\
	./src/g_degen.o		\
	./src/h_compar.o	\
	./src/i_initqd.o	\
	./src/j_quad.o		\
	./src/k_rebay.o		\
	./src/l_insert.o	\
	./src/m_putin.o		\
	./src/n_pulout.o	\
	./src/o_nearpt.o	\
	./src/p_qualty.o	\
	./src/q_grid.o

%.o: %.f90
	$(FC) -c $< -o $@ $(FFLAGS)

all: $(OBJS)
	$(AR) libtri.a $(OBJS)

clean:
	find ./ -name '*~' -exec rm {} \;
	find ./ -name '*.o' -exec rm {} \;
	find ./ -name '*.mod' -exec rm {} \;
	find ./ -name '*__genmod.f90' -exec rm {} \;
	rm -f libtri.a
