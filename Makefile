FC=mpif90
CC=gcc

FFLAGS=-O3 -c -ffree-line-length-none -DUSE_MPI
CFLAGS=-O3 -c -fPIC
#CFLAGS= -c -g -O0 -Wall -fPIC

METIS_DIR=/usr/local/
LDFLAGS=$(METIS_DIR)/lib/libmetis.a

SOURCES= my_kinddefs.F90 \
         connection_module.F90 \
         mpi_schedule_module.F90 \
         mesh_module.F90 \
         line_module.F90 \
         partition_module.F90 \
         load_mesh.F90 \
         setup_mesh.F90 \
         partitioner.F90 


INCLUDES=
OBJECTS=$(addsuffix .o,$(basename $(SOURCES)))

EXECUTABLE=partitioner

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(FC) $(OBJECTS) $(LDFLAGS) -o $@.mpi
	rm *.o *.mod

.SUFFIXES: .F90 .c .o

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) $< -o $@

.cpp.o:
	$(CXX) $(CFLAGS) $(INCLUDES) $< -o $@

.F90.o:
	$(FC) $(FFLAGS) $< -o $@
		
clean: 
	$(RM) partitioner.mpi  *.o *.mod
	
print-%: ; @echo '$(subst ','\'',$*=$($*))'
