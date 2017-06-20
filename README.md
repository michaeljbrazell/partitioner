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

3) type make

4) mpirun -np 4 ./partitioner.mpi -p 16 ../hc_mixed_ph.1.l8.ugrid 

5) number of ranks does not have to match number of partitions. In the above example 4 ranks are creating 16 partitions


Limitations and things that are hard coded:

1) Outputs files for NSU3D mcell format only. 

2) little/big endian is not specified because the TMR grid maker checks your system. Change this in load_mesh in two spots or compile with big/little endian option

3) max lines and boundary_tag is hard coded into line_module.F90 

4) lines from hemisphere cylinder grid maker are not used since lines are detected automatically for prisms and hexes. For a mesh with tetrahedral elements and lines this will need to be modified to use the lines output from TMR hemisphere cylinder grid maker.





