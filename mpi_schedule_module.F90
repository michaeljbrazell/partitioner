module mpi_schedule_module

	use my_kinddefs
	
	implicit none
#ifdef USE_MPI
	include 'mpif.h'
#endif

	integer(i4) :: rank
	integer(i4) :: num_proc
	
	
	type mpi_schedule

		integer(i4)          :: nproc_send
		integer(i4)          :: nproc_recv
		integer(i4)          :: nbuff_send
		integer(i4)          :: nbuff_recv

		integer(i4), pointer,contiguous :: iproc_send(:)  => null()
		integer(i4), pointer,contiguous :: iproc_recv(:)  => null()
		integer(i4), pointer,contiguous :: ipntr_send(:)  => null()
		integer(i4), pointer,contiguous :: ipntr_recv(:)  => null()

		integer(i4), pointer,contiguous :: ibuff_send(:)  => null()
		integer(i4), pointer,contiguous :: ibuff_recv(:)  => null()
		real(dp),    pointer,contiguous :: fbuff_send(:)  => null()
		real(dp),    pointer,contiguous :: fbuff_recv(:)  => null()

		complex(cp), pointer,contiguous :: cbuff_send(:)  => null()
		complex(cp), pointer,contiguous :: cbuff_recv(:)  => null()

		integer(i4), pointer,contiguous :: ilocal_send(:) => null()
		integer(i4), pointer,contiguous :: ilocal_recv(:) => null()

		integer(i4), pointer,contiguous ::  msgid(:)      => null()
		integer(i4), pointer,contiguous ::  istatus(:)    => null()

		integer(i4)          :: ierr_send,ierr_recv,ierr
		integer(i4)          :: itype_send,itype_recv
		integer(i4)          :: nproc_all

    end type mpi_schedule
    
    
    !----------------------------------------------------------------------------
	contains

	subroutine mpi_schedule_out(mpi0,iunit)


		type(mpi_schedule), intent(in) :: mpi0
		integer(i4),        intent(in) :: iunit
		
		!-Tmps
		integer(i4) :: i
		integer(i4) :: idum

		idum = 0

		!if (id_proc == 0) print*, 'Writing out MPI Schedule ...'

		write(iunit) mpi0%nproc_send,mpi0%nproc_recv

		do i=1,mpi0%nproc_send
		write(iunit) mpi0%iproc_send(i),mpi0%ipntr_send(i)
		enddo
		write(iunit) mpi0%ipntr_send(mpi0%nproc_send+1)

		do i=1,mpi0%nproc_recv
		write(iunit) mpi0%iproc_recv(i),mpi0%ipntr_recv(i)
		enddo
		write(iunit) mpi0%ipntr_recv(mpi0%nproc_recv+1)

		do i=1,mpi0%nbuff_send
		write(iunit) mpi0%ilocal_send(i)
		enddo
		write(iunit) idum    !Required to match pre_nsu3d

		do i=1,mpi0%nbuff_recv
		write(iunit) mpi0%ilocal_recv(i)
		enddo
		write(iunit) idum    !Required to match pre_nsu3d

	end subroutine mpi_schedule_out
  
  
  	subroutine mpi_schedule_out_reverse(mpi0,iunit)


		type(mpi_schedule), intent(in) :: mpi0
		integer(i4),        intent(in) :: iunit
		
		!-Tmps
		integer(i4) :: i
		integer(i4) :: idum

		idum = 0

		!if (id_proc == 0) print*, 'Writing out MPI Schedule ...'

		write(iunit) mpi0%nproc_recv,mpi0%nproc_send

		do i=1,mpi0%nproc_recv
		write(iunit) mpi0%iproc_recv(i),mpi0%ipntr_recv(i)
		enddo
		write(iunit) mpi0%ipntr_recv(mpi0%nproc_recv+1)

		do i=1,mpi0%nproc_send
		write(iunit) mpi0%iproc_send(i),mpi0%ipntr_send(i)
		enddo
		write(iunit) mpi0%ipntr_send(mpi0%nproc_send+1)

		do i=1,mpi0%nbuff_recv
		write(iunit) mpi0%ilocal_recv(i)
		enddo
		write(iunit) idum    !Required to match pre_nsu3d

		do i=1,mpi0%nbuff_send
		write(iunit) mpi0%ilocal_send(i)
		enddo
		write(iunit) idum    !Required to match pre_nsu3d

	end subroutine mpi_schedule_out_reverse
  
  
  
  	subroutine check_mpi(mpi_cell,part)


        type(mpi_schedule), intent(in) :: mpi_cell
        integer(i4), intent(in) :: part
        
        integer(i4) :: i,j,ind,proc_recv,proc_send
        
        proc_recv = 3
        proc_send = 3
        
        print*,' use global numbering for this test' 
        do i=1,mpi_cell%nproc_send
!        	if(mpi_cell%iproc_send(i) == proc_send) then
        	    print*,'this is part',part,'sending to',mpi_cell%iproc_send(i)
        	    ind = 1
        		do j=mpi_cell%ipntr_send(i),mpi_cell%ipntr_send(i+1)-1
        			print*,ind,mpi_cell%ilocal_send(j)
        			ind = ind + 1
        		end do
!			end if
		end do


        do i=1,mpi_cell%nproc_recv
!        	if(mpi_cell%iproc_recv(i) == proc_recv) then
        	    print*,'this is part',part,'receiving from',mpi_cell%iproc_recv(i)
        	    ind = 1
        		do j=mpi_cell%ipntr_recv(i),mpi_cell%ipntr_recv(i+1)-1
        			print*,ind,mpi_cell%ilocal_recv(j)
        			ind = ind + 1
        		end do
!			end if
		end do
		

  	
  	end subroutine check_mpi


end module mpi_schedule_module