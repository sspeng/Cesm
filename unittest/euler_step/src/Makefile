OBJS = kinds.o hybrid_mod.o physical_constants.o control_mod.o dimensions_mod.o element_mod.o hybvcoord_mod.o  \
    derivative_mod.o edgetype_mod.o slave_euler_step.o slave_euler_v.o slave_biharmonic_wk_scalar.o\
	  euler_step.o
all: euler_step
euler_step : $(OBJS)
	sw5f90 -hybrid kinds.o hybrid_mod.o physical_constants.o control_mod.o dimensions_mod.o element_mod.o hybvcoord_mod.o \
	   derivative_mod.o edgetype_mod.o slave_euler_step.o slave_euler_v.o         \
		 slave_biharmonic_wk_scalar.o euler_step.o -o ../bin/euler_step
element_mod.o: element_mod.F90
	sw5f90 -host -c element_mod.F90
kinds.o: kinds.F90
	sw5f90 -host -c kinds.F90
dimensions_mod.o: dimensions_mod.F90
	sw5f90 -host -c dimensions_mod.F90
control_mod.o: control_mod.F90
	sw5f90 -host -c control_mod.F90
physical_constants.o: physical_constants.F90
	sw5f90 -host -c physical_constants.F90
derivative_mod.o : derivative_mod.F90
	sw5f90 -host -c derivative_mod.F90
edgetype_mod.o : edgetype_mod.F90
	sw5f90 -host -c edgetype_mod.F90
hybrid_mod.o: hybrid_mod.F90
	sw5f90 -host -c hybrid_mod.F90
hybvcoord_mod.o: hybvcoord_mod.F90
	sw5f90 -host -c hybvcoord_mod.F90
slave_euler_step.o: slave_euler_step.c dma_macros.h ldm_alloc.h
	sw5cc -slave -c slave_euler_step.c -msimd
slave_euler_v.o: slave_euler_v.c dma_macros.h ldm_alloc.h
	sw5cc -slave -c slave_euler_v.c -msimd
slave_biharmonic_wk_scalar.o : slave_biharmonic_wk_scalar.c
	sw5cc -slave -c slave_biharmonic_wk_scalar.c -msimd
euler_step.o: euler_step.F90
	sw5f90 -host -c euler_step.F90
clean:
	rm *.o  *.mod ../bin/euler_step
