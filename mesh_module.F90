module mesh_module

	use my_kinddefs
	implicit none

	integer(i4) :: ntet,npyramid,nprism,nhex,ncell
	integer(i4) :: nnodes,ntri_bdry,nquad_bdry,npatch
	integer(i4) :: ntri,nquad,nface_bdry,nface
	integer(i4) :: pdegree
	
	integer(i4), allocatable :: cell_nodes(:,:)
	integer(i4), allocatable :: cell_type(:)
	integer(i4), allocatable :: cell2face(:,:)
	integer(i4), allocatable :: cell2cell(:,:)
	integer(i4), allocatable :: cell_g2l(:)
	integer(i4), allocatable :: cell_partition(:)

	integer(i4), allocatable :: face_nodes(:,:)
	integer(i4), allocatable :: face_type(:)
	integer(i4), allocatable :: face2cell(:,:)	
	integer(i4), allocatable :: face_patch(:)
		
	real(dp),    allocatable :: xgeom(:,:)
	
	real(dp) :: mode_tetra(3,56)
	real(dp) :: mode_pyr(3,14)
	real(dp) :: mode_prism(3,18)
	real(dp) :: mode_hex(3,125)
	
	contains
	
	
	subroutine mpi_broadcast_mesh()

#ifdef USE_MPI
	
		use mpi_schedule_module
		
		integer(i4) :: host,mpi_ierr


		host = 0

		call MPI_Bcast(ntet,      1, MPI_INTEGER,host,MPI_COMM_WORLD,mpi_ierr)
		call MPI_Bcast(npyramid,  1, MPI_INTEGER,host,MPI_COMM_WORLD,mpi_ierr)
		call MPI_Bcast(nprism,    1, MPI_INTEGER,host,MPI_COMM_WORLD,mpi_ierr)
		call MPI_Bcast(nhex,      1, MPI_INTEGER,host,MPI_COMM_WORLD,mpi_ierr)
		call MPI_Bcast(ncell,     1, MPI_INTEGER,host,MPI_COMM_WORLD,mpi_ierr)

		call MPI_Bcast(nnodes,     1, MPI_INTEGER,host,MPI_COMM_WORLD,mpi_ierr)
		call MPI_Bcast(ntri_bdry,  1, MPI_INTEGER,host,MPI_COMM_WORLD,mpi_ierr)
		call MPI_Bcast(nquad_bdry, 1, MPI_INTEGER,host,MPI_COMM_WORLD,mpi_ierr)
		call MPI_Bcast(npatch,     1, MPI_INTEGER,host,MPI_COMM_WORLD,mpi_ierr)

		call MPI_Bcast(ntri,       1, MPI_INTEGER,host,MPI_COMM_WORLD,mpi_ierr)
		call MPI_Bcast(nquad,      1, MPI_INTEGER,host,MPI_COMM_WORLD,mpi_ierr)
		call MPI_Bcast(nface_bdry, 1, MPI_INTEGER,host,MPI_COMM_WORLD,mpi_ierr)
		call MPI_Bcast(nface,      1, MPI_INTEGER,host,MPI_COMM_WORLD,mpi_ierr)
		call MPI_Bcast(pdegree,    1, MPI_INTEGER,host,MPI_COMM_WORLD,mpi_ierr)
		
		if(rank /= host) then
	
			allocate(cell_nodes(125,ncell))
			allocate(cell_type(ncell))
			allocate(cell_g2l(ncell))
			allocate(cell_partition(ncell))

			allocate(face_nodes(4,nface))
			allocate(face_type(nface))
			allocate(face2cell(2,nface))
			allocate(face_patch(nface))

			allocate(xgeom(3,nnodes))
	
		end if
		
		call MPI_Bcast(cell_nodes,     125*ncell,  MPI_INTEGER,host,MPI_COMM_WORLD,mpi_ierr)
		call MPI_Bcast(cell_type,      ncell,    MPI_INTEGER,host,MPI_COMM_WORLD,mpi_ierr)
		call MPI_Bcast(cell_g2l,       ncell,    MPI_INTEGER,host,MPI_COMM_WORLD,mpi_ierr)
		call MPI_Bcast(cell_partition, ncell,    MPI_INTEGER,host,MPI_COMM_WORLD,mpi_ierr)

		call MPI_Bcast(face_nodes,     4*nface,  MPI_INTEGER,host,MPI_COMM_WORLD,mpi_ierr)
		call MPI_Bcast(face_type,      nface,    MPI_INTEGER,host,MPI_COMM_WORLD,mpi_ierr)
		call MPI_Bcast(face2cell,      2*nface,  MPI_INTEGER,host,MPI_COMM_WORLD,mpi_ierr)
		call MPI_Bcast(face_patch,     nface,    MPI_INTEGER,host,MPI_COMM_WORLD,mpi_ierr)

		call MPI_Bcast(xgeom,          3*nnodes, MPI_DOUBLE,host,MPI_COMM_WORLD,mpi_ierr)
		
#endif	
		
	end subroutine
	
	subroutine deallocate_cell2stuff()
		deallocate(cell2face)
		deallocate(cell2cell)
	end subroutine deallocate_cell2stuff
	
end module mesh_module
