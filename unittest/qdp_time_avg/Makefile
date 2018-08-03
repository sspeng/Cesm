OBJS = element_mod.o qdp_time_avg_slave.o qdp_time_avg.o
all: qdp_time_avg
qdp_time_avg : $(OBJS)
	sw5f90 -hybrid qdp_time_avg_slave.o element_mod.o qdp_time_avg.o -g -o ./bin/qdp_time_avg
qdp_time_avg_slave.o: qdp_time_avg_slave.c
	sw5cc -slave -c qdp_time_avg_slave.c -msimd
element_mod.o: element_mod.F90
	sw5f90 -host -c element_mod.F90
qdp_time_avg.o: qdp_time_avg.F90
	sw5f90 -host -c  qdp_time_avg.F90
clean:
	rm *.o  *.mod ./bin/qdp_time_avg
