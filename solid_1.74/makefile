# Makefile for FEMULA
#    by Jeeho Lee
# Last Modification on April 15, 2011

.SUFFIXES: .f .f90 .for .o

ifeq ($(mode),debug)
FFLAGS = -g -debug full -ftrapuv -warn all -check all -traceback -fpe0 -free -openmp -real-size 64 -mkl
else
FFLAGS = -O3 -fpe0 -free -openmp -real-size 64 -mkl
endif

# Include Needed for Compile
IFLAGS = -I ./Lib

# Library Needed for Linking
LFLAGS = ./Lib/libtecio.a

METIS=./Lib/libmetis.a

LIB = ./main/fem_gxl_lib.o ./main/fem_gxl_lib_geo.o

MODULES0 = ./main/file_info_house.o ./main/system_info_house.o ./main/fem_gxl_math.o

MODULES =  ./main/pdomain.o ./main/tdomain.o ./main/sdomain.o ./main/elements.o  ./main/element_spec.o ./main/materials.o ./main/coarseproblem.o ./postprocess/output_domain.o ./3D/surface/quick_sort.o ./3D/remesh/remesh_domain.o

SINGLE = ./main/buildmodel.o ./main/contact_control.o ./main/fem_gxl_foro.o ./main/fem_gxl_initialize.o ./main/fem_gxl_prima.o ./main/fem_gxl.o

ELEMENTS =  ./elements/elmt_lib.o ./elements/elmt_fn_lib.o ./elements/elmt_2d.o ./elements/elmt_2d_mixed.o ./elements/elmt_2d_fn.o ./elements/elmt_2d_fn_mx.o ./elements/elmt_3d.o ./elements/elmt_3d_fn.o ./elements/elmt_3d_fn_mx.o ./elements/str_tang_infn.o ./elements/str_tang_fn.o ./elements/spring_elmt.o

MATERIALS = ./materials/linearelastic.o ./materials/viscoelastic.o ./materials/gmnrv.o ./materials/ozupek.o ./materials/umat_ozupek.o ./materials/nonlinearuniaxial.o

INTERFACE = ./main/element_interface.o ./main/material_initialize.o ./main/build_interface.o ./main/mesh_interface.o ./main/recover_step.o ./main/fsi_interface.o 

POSTPROCESS = ./postprocess/ext_quad_data.o ./postprocess/nodal_interpolation.o ./postprocess/output_main.o ./postprocess/utility.o

ALE = ./ale/ale_domain.o ./ale/ale_control.o ./ale/ale_lib.o ./ale/fvm_routine.o ./ale/restore_states.o ./ale/remesh_interpolation.o

FSI = ./fsi/fsi_domain.o ./fsi/fsi_iter_set.o ./fsi/pload_interpolation.o

SURFACE = ./surface/surface_domain.o ./surface/mesh_smooth.o ./surface/move_surface.o ./surface/post_build_interface.o ./surface/remesh_routine.o ./surface/remesh_routine_3D.o ./surface/separation_routine.o ./surface/update_surface.o ./surface/surface_control.o ./surface/surface_lib.o ./fsi/set_fsi_inte.o

MESH = ./mesh/mesh_lib.o ./mesh/angle_check.o ./mesh/angle_recheck.o ./mesh/bound_soft.o ./mesh/calc_d.o ./mesh/check_ele_angle.o ./mesh/clock_check.o ./mesh/combine_two_ele.o ./mesh/con_node_cross.o ./mesh/continue_end_smooth.o ./mesh/cross_check.o ./mesh/cross_near.o ./mesh/divide_ele.o ./mesh/divide_point.o ./mesh/end_check.o ./mesh/inner_soft.o ./mesh/inout.o ./mesh/len_smooth.o ./mesh/line_revise.o ./mesh/make_ele.o ./mesh/mesh_generation.o ./mesh/mesh_gen.o ./mesh/one_ele.o ./mesh/partition.o ./mesh/paving.o ./mesh/point_summation.o ./mesh/remake_bound.o ./mesh/reset_bound.o ./mesh/spline.o ./mesh/thin.o ./mesh/sep_mesh_gen.o

EXTERNALS = ./externals/cohesive1.o ./externals/cohesive1_2d.o 

3D_ALE = ./3D/ale/ale_lib_3D.o ./3D/ale/set_ale_region.o

3D_REMESH = ./3D/remesh/b_rep_lib.o ./3D/remesh/find_edge.o ./3D/remesh/make_hexa_mesh.o ./3D/remesh/remesh_3D_main.o  ./3D/remesh/remesh_ale_region.o ./3D/remesh/remesh_lib.o ./3D/remesh/remesh_smoothing.o ./3D/remesh/set_b_rep.o

3D_SURFACE = ./3D/surface/calc_elem_jacobian.o ./3D/surface/mesh_smooth_3D.o ./3D/surface/set_bound_node_lv.o ./3D/surface/surface_lib_3D.o 

.f.o:
	ifort -c $(FFLAGS) $*.f -o $@  $(IFLAGS)

.f90.o:
	ifort -c $(FFLAGS) $*.f90 -o $@  $(IFLAGS)

.for.o:
	ifort -c -O3 -fpe0 -openmp -real-size 64 -mkl $*.for -o $@ $(IFLAGS)


./libsolid.a: $(MODULES0) $(MODULES) $(LIB) $(FSI) $(SURFACE) $(ALE) $(SINGLE) $(MATERIALS) $(ELEMENTS) $(INTERFACE) $(POSTPROCESS) $(MESH) $(EXTERNALS) $(3D_ALE) $(3D_REMESH) $(3D_SURFACE) $(LFLAGS)
	$(AR) rcv $@ $^
	ranlib $@
clean:
	rm -f *.a
	rm -f *.o
	rm -f *.mod
	rm -f *.f90
	rm -f ./main/*.o
	rm -f ./main/*.mod
	rm -f ./elements/*.o
	rm -f ./elements/*.mod
	rm -f ./materials/*.o
	rm -f ./materials/*.mod
	rm -f ./mesh/*.o
	rm -f ./mesh/*.mod
	rm -f ./postprocess/*.o
	rm -f ./postprocess/*.mod
	rm -f ./ale/*.o
	rm -f ./ale/*.mod
	rm -f ./fsi/*.o
	rm -f ./fsi/*.mod
	rm -f ./surface/*.o
	rm -f ./surface/*.mod
	rm -f ./externals/*.o
	rm -f ./externals/*.mod
	rm -f ./3D/ale/*.o
	rm -f ./3D/ale/*.mod
	rm -f ./3D/remesh/*.o
	rm -f ./3D/remesh/*.mod
	rm -f ./3D/surface/*.o
	rm -f ./3D/surface/*.mod
	rm -f *__genmod.f90
