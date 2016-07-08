# partitioner

Thanks for using the mixed element line based mesh partitioner!

Some features include:

-Input an unstructured mixed element mcell.unf or gmsh .msh file

-Detects cell based lines (prisms and hexes) and collapse the graph

-Partitions the collapsed graph using Metis (in serial) with a simple weighting for element type and line length

-Outputs the partitioned mesh into separate files with mpi schedules

-Designed for cell-based discretizations

-Some support for high order elements (p=2 tetrahedra) use gmsh manual to add new high order element types


How to use:

1) install metis first http://glaros.dtc.umn.edu/gkhome/metis/metis/download

2) modify makeall.sh to use your libmetis.a static library

3) comment/uncomment to use appropiate compiler in makeall.sh


