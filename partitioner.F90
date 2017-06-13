program main

    use my_kinddefs
    
    implicit none
    
    character(len=255) :: arg,inputfile
    integer(i4) :: count,i,partition_flag,mesh_type
    logical :: is_numeric 
      
   
    ! set default partition number
    partition_flag = 16
    
    ! find number of command arguments
    count = command_argument_count()

    ! read in commands arguments
    do i=1,count
        call get_command_argument(i,arg)
        !write(*,*) trim(arg)

        if(arg == '-h' .or. arg == '-help' .or. arg == 'help') then
            write(*,*) ' the following command arguments are allowed'
            write(*,*) ' blah.mcell.unf : specify mesh file name'
            write(*,*) ' -p 8 : run partition routines'
            stop
        end if

        
        ! look for partition flag and number
        if(arg == '-p' ) then 
                        
            if(i < count) then
                call get_command_argument(i+1,arg)
            end if
            
            if(is_numeric(arg) .eqv. .true.) then            
                Read( arg, '(i10)' ) partition_flag                
            end if
            
        end if 
        
        ! look for input file
        arg = adjustr(arg)
        if(arg(247:255) == 'mcell.unf') then
            inputfile = adjustl(arg) 
            mesh_type = 1           
        else if (arg(250:255) == '.ugrid') then
            inputfile = adjustl(arg)  
            mesh_type = 2
        else if (arg(252:255) == '.msh') then
            inputfile = adjustl(arg)  
            mesh_type = 3
  end if
        
    end do

    call go(inputfile,partition_flag,mesh_type)
    
    
    
end program


  subroutine go(inputfile,nparts,mesh_type)

    use my_kinddefs
    use line_module
    use partition_module
    use mpi_schedule_module
    
    implicit none
    
    character(len=255), intent(in) :: inputfile
    integer(i4), intent(in) :: nparts
    integer(i4), intent(in) :: mesh_type
    
#ifdef USE_MPI  
    integer(i4) :: mpi_ierr
    call MPI_INIT(mpi_ierr)    
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_ierr) 
    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_proc, mpi_ierr)

!    if(num_proc /= nparts) then
!    	print*,'ahh partitions do not match mpi ranks'
!	 end if   
	
#else
    rank = 0
    num_proc = 1 
#endif

   
    if(rank==0) then
      print*,'loading mesh'
      call load_mesh(inputfile,mesh_type)
 
      print*,'setting up mesh'
      call setup_mesh()
 
      print*,'finding lines'
      call setup_lines()    
 
      print*,'partitioning into ', nparts
      call partition_mesh(nparts) 
 
    !    print*,'output mesh'
    !    call output_mesh()
 
      print*,'global to local'
      call setup_global_to_local(nparts)
 
      print*,'deallocate stuff'
      call deallocate_lines()
      call deallocate_cell2stuff()
    end if
    
    if(rank==0) print*,'mpi broadcast mesh'
    call mpi_broadcast_mesh()
  
    if(rank==0) print*,'output partition mesh'    
    call output_partition_mesh_files(nparts)   
    

#ifdef USE_MPI  
    call mpi_barrier(MPI_COMM_WORLD,mpi_ierr)
    if(rank==0) print*,'done'
    call MPI_FINALIZE(mpi_ierr)
#else
    print*,'done'
#endif
    

    
end subroutine go


FUNCTION is_numeric(string)
  IMPLICIT NONE
  CHARACTER(len=*), INTENT(IN) :: string
  LOGICAL :: is_numeric
  REAL :: x
  INTEGER :: e
  READ(string,*,IOSTAT=e) x
  is_numeric = e == 0
END FUNCTION is_numeric


