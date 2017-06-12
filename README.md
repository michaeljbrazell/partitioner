# partitioner


Thanks for using the mixed element line based mesh partitioner!

Some features include:

-Input an unstructured mixed element mcell.unf, gmsh .msh file, or a ugrid file

-Detects cell based lines (prisms and hexes) and collapse the graph

-Partitions the collapsed graph using Metis (in serial) with a simple weighting for element type and line length

-Outputs the partitioned mesh into separate files with mpi schedules

-Designed for cell-based discretizations

-Some support for high order elements (p=2 tetrahedra) use gmsh manual to add new high order element types


How to use:

1) install metis first http://glaros.dtc.umn.edu/gkhome/metis/metis/download

2) modify Makefile to use your libmetis.a static library

3) mpirun -np 4 ./partitioner.mpi -p 16 ../hc_mixed_ph.1.l8.ugrid 


Limitations and things that are hard coded:

1) Outputs files for NSU3D mcell format only. 

2) if you want a partition for each cell then output cell_partition(cell_id) which is filled in after partition_mesh is called

3) little/big endian is not specified because the TMR grid maker checks your system

4) max lines and boundary_tag is hard coded into line_module.F90 





