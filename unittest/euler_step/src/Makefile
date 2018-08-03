OBJS = kinds.o hybrid_mod.o dimensions_mod.o element_mod.o hybvcoord_mod.o derivative_mod.o slave_euler_step.o euler_step.o
all: euler_step
euler_step : $(OBJS)
	sw5f90 -hybrid kinds.o hybrid_mod.o dimensions_mod.o element_mod.o hybvcoord_mod.o derivative_mod.o slave_euler_step.o euler_step.o -o ../bin/euler_step
element_mod.o: element_mod.F90
	sw5f90 -host -c element_mod.F90
kinds.o: kinds.F90
	sw5f90 -host -c kinds.F90
dimensions_mod.o: dimensions_mod.F90
	sw5f90 -host -c dimensions_mod.F90
derivative_mod.o : derivative_mod.F90
	sw5f90 -host -c derivative_mod.F90
hybrid_mod.o: hybrid_mod.F90
	sw5f90 -host -c hybrid_mod.F90
hybvcoord_mod.o: hybvcoord_mod.F90
	sw5f90 -host -c hybvcoord_mod.F90
slave_euler_step.o: slave_euler_step.c
	sw5cc -slave -c slave_euler_step.c -msimd
euler_step.o: euler_step.F90
	sw5f90 -host -c euler_step.F90
clean:
	rm *.o  *.mod ../bin/euler_step
