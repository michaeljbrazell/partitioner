rm *.mod partitioner.serial partitioner.mpi


# use this for gnu compilers
# note the libmetis.a file linked at the end change it to your library path
gfortran -O2  -ffree-line-length-none  my_kinddefs.F90 connection_module.F90 mpi_schedule_module.F90 mesh_module.F90 line_module.F90  partition_module.F90  load_mesh.F90 setup_mesh.F90 partitioner.F90 /usr/local/lib/libmetis.a -o partitioner.serial
rm *.mod
mpif90 -O2  -DUSE_MPI -ffree-line-length-none  my_kinddefs.F90 connection_module.F90 mpi_schedule_module.F90 mesh_module.F90 line_module.F90  partition_module.F90  load_mesh.F90 setup_mesh.F90 partitioner.F90 /usr/local/lib/libmetis.a -o partitioner.mpi

# use this for intel compilers
##ifort -O2  my_kinddefs.F90 connection_module.F90 mpi_schedule_module.F90 mesh_module.F90 line_module.F90  partition_module.F90  load_mesh.F90 setup_mesh.F90 partitioner.F90 -o partitioner.serial -L/glade/u/home/mbrazell/Codes/metis-5.1.0/build/Linux-x86_64/libmetis -lmetis
##rm *.mod
##mpif90 -O2  -DUSE_MPI  my_kinddefs.F90 connection_module.F90 mpi_schedule_module.F90 mesh_module.F90 line_module.F90  partition_module.F90  load_mesh.F90 setup_mesh.F90 partitioner.F90  -o partitioner.mpi -L/glade/u/home/mbrazell/Codes/metis-5.1.0/build/Linux-x86_64/libmetis -lmetis

./partitioner.serial -p 16 gmsh/arcp2.msh
./partitioner.serial -p 4 gmsh/arc.msh
./partitioner.serial -p 5 gmsh/pripyrtet.msh
./partitioner.serial -p 6 gmsh/transfinite.msh

mpirun -np 4 ./partitioner.mpi -p 8 bump2d.mcell.unf


