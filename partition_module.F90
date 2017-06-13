module partition_module

  use my_kinddefs

  implicit none

  
  contains


  subroutine partition_mesh(nparts)

    use line_module
    use mesh_module
    use connection_module, only : node_per_element

    implicit none

    integer, intent(in) :: nparts

    integer(i4) :: i,j,last_cell,num_nodes,ind,cellnode    
    integer :: ne,nn,ncommon,objval
    integer, allocatable :: eptr(:),eind(:),epart(:),npart(:)
    !        integer(i4) :: opts(40)
    integer,pointer :: opts
    real,pointer :: tpwgts
    !        integer,allocatable :: vwgt(:),vsize(:)
    integer,pointer :: vwgt(:),vsize(:)=>null()

    integer :: part_count(nparts)

    !        vwgt=>null()
    !        vsize=>null()

    !        print*,'made it to partition mesh module'

    tpwgts=>null()
    opts=>null()


    if(nparts < 2) then
        return
    end if
 
    nn = max(nnodes,8*ncell) 
    ne = num_lines
        
!   print*,'number of lines',ne
!   print*,'number of nodes',nnodes
!   print*,'number of cells',ncell
        
    allocate(eptr(ne+1))
    allocate(eind(nn))
    allocate(vwgt(ne)) 
    allocate(vsize(ne))
!   allocate(tpwgts(nparts))
    allocate(epart(ne))
    allocate(npart(nn))
!   allocate(opts(1:40))
        
    print*,'made it past metis allocation'

    ind = 1
    eptr(1) = 0
    do i = 1,num_lines
        
      ! find the last cell in a line
      last_cell = line(i)%cell(line(i)%num_cells)
      ! num nodes on last cell
      num_nodes = node_per_element(cell_type(last_cell))
      ! pointer array to cell nodes
      eptr(i+1) = eptr(i) + num_nodes
      
      ! cell to node connectivity
      do j=1,num_nodes
          eind(ind) = cell_nodes(j,last_cell)
          ind = ind+1
      end do
      
      ! weight is based on number of cell in the line and number of modes of each cell
      ! this assumes all elements are the same within a line
      vwgt(i) = line(i)%num_cells*num_nodes!*(numfield*cell(last_cell)%sol_nmode)**2
      vsize(i) = 1!line(i)%num_cells!*numfield*cell(last_cell)%sol_nmode
      
    end do

!   print*,'length of eind',ind-1
        
    ncommon = 3
        
!        tpwgts = 1.0_dp/real(nparts,dp)

!        call METIS_SetDefaultOptions(opts)
        
!        opts(1) =  1 ! use kway
!        opts(12) = 1 ! force contiguous
!        opts(18) = 1! use fortran numbering


    call METIS_PartMeshDual(ne,nn,eptr,eind,vwgt,vsize,ncommon,nparts,tpwgts,opts,objval,epart,npart)

    print*,'made it past metis'

        
!        epart = epart+1
        
        

    do i=1,num_lines
      do j=1,line(i)%num_cells
        cell_partition(line(i)%cell(j)) =  epart(i)+1
      end do
    end do        

    part_count(:) = 0
    do i=1,ncell
      part_count(cell_partition(i)) = part_count(cell_partition(i)) + 1
    end do
       
!       print*, 100.0_dp*real(part_count,dp)/real(ncell,dp)

    deallocate(eptr,eind,epart,npart)
  
  
  end subroutine partition_mesh



  subroutine setup_global_to_local(nparts)

    use mesh_module
    use line_module
  
    integer(i4),intent(in) :: nparts
  
    integer(i4) :: i,j,e,p,ind,t,etype,n
    integer(i4) :: partition_count(nparts)
  
    partition_count(:) = 0		
  
  
    do t=4,7
      do i = 1,num_lines            
        do j=1,line(i)%num_cells
          e = line(i)%cell(j)
          etype = cell_type(e)
          if(etype == t) then
            p = cell_partition(e)
            partition_count(p) = partition_count(p) + 1            	
            cell_g2l(e) = partition_count(p)
          end if
        end do            
      end do 
    end do       
        
!        do i=1,ncell
!        	if(cell_type(i) == 6 .and. cell_partition(i) == 1) then
!        	print*,i,cell_type(i),cell_partition(i),cell_g2l(i)
!        	end if
!        end do       


  end subroutine setup_global_to_local




  subroutine output_partition_mesh_files(nparts)

    use mpi_schedule_module, only : rank, num_proc
  
    integer(i4), intent(in) :: nparts
  
    integer(i4) :: part
    integer(i4) :: parts_per_rank
    integer(i4) :: remainder_parts
  
    parts_per_rank = int(real(nparts,dp)/real(num_proc,dp),i4)
  
    do part=parts_per_rank*rank+1,parts_per_rank*(rank+1)
      if(rank==0) print*,int(real(part,dp)/real(parts_per_rank,dp)*100,i4),'%'
      call output_partition_mesh_file(part,nparts)
    end do
  
    remainder_parts = nparts-parts_per_rank*num_proc
    if(rank < remainder_parts) then
      if(rank==0) print*,'101%'
      part = nparts-rank
      call output_partition_mesh_file(part,nparts)
    end if


    
  end subroutine output_partition_mesh_files




  subroutine output_partition_mesh_file(part,nparts)

    use mesh_module
    use mpi_schedule_module
    use connection_module
  
    integer(i4), intent(in) :: part,nparts
  
    integer(i4) :: i,j,t,n,ind,e,el,er,pl,pr,etype,pc,tm
    integer(i4) :: unitnum,dummy
  
    integer(i4) :: ncell_real,ncell_ghost
    integer(i4) :: ncell_type_real(7), ncell_type_ghost(7)		
    integer(i4) :: cell_ghost(ncell)
    integer(i4) :: cell_real(ncell) ! could probably allocate less
    integer(i4) :: cell_ghost_g2l(ncell)
  
    integer(i4) :: part_face_tag(nface)
    integer(i4) :: npart_face
    integer(i4) :: part_face_sorted(nface) ! could allocate this later as length of npart_face
    integer(i4) :: part_face_sorted_ptr(nparts) 
    integer(i4) :: part_face_procs(nparts) 
  
    integer(i4) :: bdry_face_tag(nface)
    integer(i4) :: nbdry_face,nbdry_tri,nbdry_quad
    
    integer(i4) :: ntris,nquads
  
    integer(i4) :: cell_tag(ncell) ! could be changed to logical
    integer(i4) :: node_g2l(nnodes)
    integer(i4) :: node_l2g(nnodes)
    integer(i4) :: nnodes_real,nnodes_ghost
  
    integer(i4) :: proc_recv(nparts)
    integer(i4) :: proc_send(nparts)
  
    integer(i4) :: proc_recv_reorder(nparts)
    integer(i4) :: proc_send_reorder(nparts)
  
    character(len=255) :: part_string,nparts_string,filename
    logical :: file_exists
  
    type(mpi_schedule) :: mpi_cell
  
  
    cell_tag(:) = 0	
    part_face_tag(:) = 0
    npart_face = 0
    bdry_face_tag(:) = 0
    nbdry_face = 0
  
    nbdry_tri = 0
    nbdry_quad = 0
    ntris = 0
    nquads = 0
      
    ! tag any cell that is connected to this partition
    ! also tag partition faces and store them
    do i=1,nface
  
      el = face2cell(1,i)
      pl = cell_partition(el)
      er = face2cell(2,i)
      etype = face_type(i)

      ! if boundary face and on partition then tag			
      if(er < 1) then
        if(pl == part) then
          cell_tag(el) = 1				
        
          if(etype == 2) then
            nbdry_tri = nbdry_tri+1
            ntris = ntris+1
          else if(etype ==3) then
            nbdry_quad = nbdry_quad+1
            nquads = nquads+1
          end if
        
          nbdry_face = nbdry_face + 1
          bdry_face_tag(nbdry_face) = i
        end if 
        ! cycle loop because it is a boundary face
        cycle
      end if
    
      pr = cell_partition(er)
      if(pl == part .or. pr == part) then
        cell_tag(el) = 1
        cell_tag(er) = 1
      
        if(etype == 2) then
          ntris = ntris+1
        else if(etype ==3) then
          nquads = nquads+1
        end if
        
        ! partition boundary
        if(pl /= pr) then
          npart_face = npart_face + 1
          part_face_tag(npart_face) = i
        else
      
        end if
      end if
    
    end do
  
    if(nbdry_tri+nbdry_quad /= nbdry_face) then
      print*,'boundary face error'
    end if
  

    ! count number of partition faces for each proc
    part_face_procs(:) = 0		
    do i=1,npart_face
      e = part_face_tag(i)
      el = face2cell(1,e)
      pl = cell_partition(el)
      er = face2cell(2,e)
      pr = cell_partition(er)
    
      if(pl == part) then
        pc = pr
      else
        pc = pl
      end if

      part_face_procs(pc) = part_face_procs(pc) + 1
    
    end do
    
    ! use this array to add faces in order of proc
    part_face_sorted_ptr(:) = 0	
    do i=1,nparts-1			
      part_face_sorted_ptr(i+1) = part_face_sorted_ptr(i) + part_face_procs(i)
    end do
  
    ! need to sort partition faces by processor first
    part_face_sorted(:) = 0 
    do i=1,npart_face
      e = part_face_tag(i)
      el = face2cell(1,e)
      pl = cell_partition(el)
      er = face2cell(2,e)
      pr = cell_partition(er)
    
      if(pl == part) then
        pc = pr
      else
        pc = pl
      end if
      
      part_face_sorted_ptr(pc) = part_face_sorted_ptr(pc) + 1
      part_face_sorted(part_face_sorted_ptr(pc)) = e
    end do
  

    ! move into temp storage
    do i=1,npart_face
      part_face_tag(i) = part_face_sorted(i)			
    end do
  
    ind = 1
    do i=1,nparts
      if(part_face_procs(i) > 0) then
        call insort(part_face_procs(i),part_face_sorted(ind))
        ind = ind + part_face_procs(i)
      end if
    end do

    do i=1,npart_face
      if(part_face_tag(i)-part_face_sorted(i) /= 0) then
        print*,'a face was sorted'
        print*,part_face_tag(i)-part_face_sorted(i)
      end if
    end do
  

      
    ! find all real and ghost cells
    ncell_real = 0
    ncell_type_real(:) = 0

    do i=1,ncell
      if(cell_tag(i) == 1) then	
        etype = cell_type(i)
        if(cell_partition(i) == part) then
          ! add to real cells
          ncell_type_real(etype) = ncell_type_real(etype) + 1
          ncell_real = ncell_real + 1	
        
          ! use the global to local numbering so that the 
          ! real cells are ordered correctly
          cell_real(cell_g2l(i)) = i
        end if
      end if
    end do
  
  
    ncell_ghost = 0
    ncell_type_ghost(:) = 0
    cell_ghost_g2l(:) = 0
    do t=4,7
      do i=1,ncell
        etype = cell_type(i)
        if(cell_tag(i) == 1 .and. etype == t) then	
          if(cell_partition(i) /= part) then
            ! add to ghost cells
            ncell_type_ghost(etype) = ncell_type_ghost(etype) + 1
            ncell_ghost = ncell_ghost + 1
            ! could re-arrange these so they are aligned with the buffer
            cell_ghost(ncell_ghost) = i 
            cell_ghost_g2l(i) = ncell_ghost
          end if
        end if
      end do
    end do
      
    ! add in real cells to this count
    cell_ghost_g2l(:) = cell_ghost_g2l(:) + ncell_real
  
  
    if(2*(ntris+nquads) /= 4*ncell_type_real(4)+5*ncell_type_real(5)+5*ncell_type_real(6)+6*ncell_type_real(7)+nbdry_face+npart_face) then
      print*,'ahh total faces to not match'
      stop
    end if
  
  
    if(ncell_real /= sum(ncell_type_real)) then
      print*,'real cell sum error'
      stop
    end if
    if(ncell_ghost /= sum(ncell_type_ghost)) then
      print*,'ghost cell sum error'
      stop
    end if
  
    ! tag ghost nodes and processors that ghosts are on
    node_g2l(:) = 0	
    node_l2g(:) = 0
    nnodes_real = 0
    do i=1,ncell_real		
      e = cell_real(i)
      etype = cell_type(e)			
      do j=1,node_per_element(etype)
        if(node_g2l(cell_nodes(j,e)) == 0) then	
          nnodes_real = nnodes_real+1		
          node_g2l(cell_nodes(j,e)) = nnodes_real
          node_l2g(nnodes_real) = cell_nodes(j,e)					
        end if
      end do
    end do
  
    nnodes_ghost = 0
    do i=1,ncell_ghost		
      e = cell_ghost(i)
      etype = cell_type(e)			
      do j=1,node_per_element(etype)	
        if(node_g2l(cell_nodes(j,e)) == 0) then								
          nnodes_ghost = nnodes_ghost+1		
          node_g2l(cell_nodes(j,e)) = nnodes_real+nnodes_ghost
          node_l2g(nnodes_ghost+nnodes_real) = cell_nodes(j,e)
        end if	
      end do				
    end do
  

    ! tag all procs that this partition is touching
    proc_recv(:) = 0
    do i=1,ncell_ghost		
      e = cell_ghost(i)				
      proc_recv(cell_partition(e)) = 1			
    end do
      
  
    ! the number of procs that recv/send is equal
    mpi_cell%nproc_recv = sum(proc_recv)		
    mpi_cell%nproc_send = mpi_cell%nproc_recv
  
    allocate(mpi_cell%iproc_send(mpi_cell%nproc_send))
    allocate(mpi_cell%ipntr_send(mpi_cell%nproc_send+1))
    allocate(mpi_cell%iproc_recv(mpi_cell%nproc_recv))
    allocate(mpi_cell%ipntr_recv(mpi_cell%nproc_recv+1))
    mpi_cell%ipntr_send(:) = 0
    mpi_cell%ipntr_recv(:) = 0
  
        ! shift procs down
    ind = 0
    do i=1,nparts
      if(proc_recv(i) == 1) then	
        ind = ind + 1	
        mpi_cell%iproc_recv(ind) = i
        mpi_cell%iproc_send(ind) = i
      end if
    end do
  
    if(ind /= mpi_cell%nproc_recv) then
      print*,'nproc recv error'
      stop
    end if

  
    ! tag all cells along partition faces
    ! for sending and receiving
    proc_recv(:) = 0
    proc_send(:) = 0

    ind = 0
    do n=1,nparts
  
      if(part_face_procs(n) == 0) then
        cycle
      end if
    
      cell_tag(:) = 0
      i = ind
      do j=1,part_face_procs(n)
        i = i+1	
        el = face2cell(1,part_face_tag(i))
        pl = cell_partition(el)
        er = face2cell(2,part_face_tag(i))
        pr = cell_partition(er)
    
        if(pl == part) then
          if(cell_tag(er) == 0) then
            cell_tag(er) = 1
            proc_recv(pr) = proc_recv(pr) + 1
          end if
        else if(pr == part) then
          if(cell_tag(el) == 0) then
            cell_tag(el) = 1
            proc_recv(pl) = proc_recv(pl) + 1
          end if
        else
          print*,'partition face loop error'
          stop
        end if
      end do
    
      cell_tag(:) = 0
      i = ind
      do j=1,part_face_procs(n)
        i = i+1
        el = face2cell(1,part_face_tag(i))
        pl = cell_partition(el)
        er = face2cell(2,part_face_tag(i))
        pr = cell_partition(er)
    
        if(pl == part) then	
          if(cell_tag(el) == 0) then
            cell_tag(el) = 1
            proc_send(pr) = proc_send(pr) + 1
          end if				
        else if(pr == part) then
          if(cell_tag(er) == 0) then
            cell_tag(er) = 1
            proc_send(pl) = proc_send(pl) + 1
          end if
        else
          print*,'partition face loop error'
          stop
        end if
      end do
    
      ind = ind + part_face_procs(n)
    
    end do
    
  
  
    ! sum to get size of buffer
    mpi_cell%nbuff_send = sum(proc_send)
    mpi_cell%nbuff_recv = sum(proc_recv)
  
    allocate(mpi_cell%ilocal_send(mpi_cell%nbuff_send))
    allocate(mpi_cell%ilocal_recv(mpi_cell%nbuff_recv))
    mpi_cell%ilocal_send(:) = 0
    mpi_cell%ilocal_recv(:) = 0
  

    ! fill in proc pointer arrays
    mpi_cell%ipntr_send(:) = 1
    do j = 1,mpi_cell%nproc_send
      do i = j+1,mpi_cell%nproc_send+1
        mpi_cell%ipntr_send(i) = mpi_cell%ipntr_send(i) + proc_send(mpi_cell%iproc_send(j))
      end do
    end do
  
    if(sum(proc_send)+1 /= mpi_cell%ipntr_send(mpi_cell%nproc_send+1)) then
      print*, ' mpi_cell%ipntr_send error '
      stop
    end if

    mpi_cell%ipntr_recv(:) = 1
    do j = 1,mpi_cell%nproc_recv
      do i = j+1,mpi_cell%nproc_recv+1
        mpi_cell%ipntr_recv(i)	 = 	mpi_cell%ipntr_recv(i) + proc_recv(mpi_cell%iproc_recv(j))
      end do
    end do
  
    if(sum(proc_send)+1 /= mpi_cell%ipntr_send(mpi_cell%nproc_send+1)) then
      print*, ' mpi_cell%ipntr_recv error '
      stop
    end if
  
  
    ! use proc_send_reorder to find proc numbers
    ! use proc_send to add elements
    do i=1,mpi_cell%nproc_send
      proc_send_reorder(mpi_cell%iproc_send(i)) = i
      proc_send(i) = mpi_cell%ipntr_send(i)
    end do		
    do i=1,mpi_cell%nproc_recv
      proc_recv_reorder(mpi_cell%iproc_recv(i)) = i
      proc_recv(i) = mpi_cell%ipntr_recv(i)
    end do
  
  


    ! tag all cells along partition faces
    ! for sending and receiving
    ind = 0
    do n=1,nparts
  
      if(part_face_procs(n) == 0) then
        cycle
      end if
    
      cell_tag(:) = 0
      i = ind
      do j=1,part_face_procs(n)
        i = i+1	
        el = face2cell(1,part_face_tag(i))
        pl = cell_partition(el)
        er = face2cell(2,part_face_tag(i))
        pr = cell_partition(er)
    
        if(pl == part) then
          if(cell_tag(er) == 0) then
            cell_tag(er) = 1
            mpi_cell%ilocal_recv(proc_recv(proc_recv_reorder(pr))) = cell_ghost_g2l(er)
            proc_recv(proc_recv_reorder(pr)) = proc_recv(proc_recv_reorder(pr)) + 1
          end if
        else if(pr == part) then
          if(cell_tag(el) == 0) then
            cell_tag(el) = 1
            mpi_cell%ilocal_recv(proc_recv(proc_recv_reorder(pl))) = cell_ghost_g2l(el)
            proc_recv(proc_recv_reorder(pl)) = proc_recv(proc_recv_reorder(pl)) + 1
          end if
        else
          print*,'partition face loop error'
          stop
        end if
      end do
    
      cell_tag(:) = 0
      i = ind
      do j=1,part_face_procs(n)
        i = i+1
        el = face2cell(1,part_face_tag(i))
        pl = cell_partition(el)
        er = face2cell(2,part_face_tag(i))
        pr = cell_partition(er)
    
        if(pl == part) then	
          if(cell_tag(el) == 0) then
            cell_tag(el) = 1
            mpi_cell%ilocal_send(proc_send(proc_send_reorder(pr))) = cell_g2l(el)					
            proc_send(proc_send_reorder(pr)) = proc_send(proc_send_reorder(pr)) + 1	
          end if				
        else if(pr == part) then
          if(cell_tag(er) == 0) then
            cell_tag(er) = 1
            mpi_cell%ilocal_send(proc_send(proc_send_reorder(pl))) = cell_g2l(er)					
            proc_send(proc_send_reorder(pl)) = proc_send(proc_send_reorder(pl)) + 1	
          end if
        else
          print*,'partition face loop error'
          stop
        end if
      end do
    
      ind = ind + part_face_procs(n)
    
    end do
  
    do i=1,mpi_cell%nbuff_send
      if(mpi_cell%ilocal_send(i) > ncell_real) then
        print*,' something wrong with ilocal send, too big'
      end if
    end do
  
    do i=1,mpi_cell%nbuff_recv
      if(mpi_cell%ilocal_recv(i)  <= ncell_real .or. mpi_cell%ilocal_recv(i)> ncell_real+ncell_ghost) then
        print*,' something wrong with ilocal recv, too big or too small'
      end if
    end do
  
  
    ! call check_mpi(mpi_cell,part)

    ! write out the partition mesh file
  
    write(unit=part_string, fmt='(I6)') part-1
    write(unit=nparts_string, fmt='(I6)') nparts
    unitnum = 35
    filename = 'part.'//trim(adjustl(nparts_string))
  
    inquire(file=filename, exist=file_exists)
    !inquire(directory=filename, exist=file_exists) ! use this for intel
  
    if(file_exists .eqv. .false.) then
      call system('mkdir '//filename)
    end if
  
    filename = 'part.'//trim(adjustl(nparts_string))//'/cell.'//trim(adjustl(part_string))


    open (unit=unitnum,file=filename, form='unformatted',access='STREAM',convert='big_endian')
  
    ! global totals
    write(unitnum) ntet,npyramid,nprism,nhex
    write(unitnum) nnodes,ntri_bdry,nquad_bdry,npatch,1,1
    ! local totals
    write(unitnum) ncell_type_real(4),ncell_type_real(5),ncell_type_real(6),ncell_type_real(7),nbdry_tri,nbdry_quad,nnodes_real
    write(unitnum) ncell_type_ghost(4),ncell_type_ghost(5),ncell_type_ghost(6),ncell_type_ghost(7),nnodes_ghost
    write(unitnum) ntris,nquads
  
    ! write out real cell to node data
    do j=1,ncell_real
      e = cell_real(j)
      etype = cell_type(e)
      do i=1,node_per_element(etype)
        write(unitnum) node_g2l(cell_nodes(i,e))
      end do
    end do
  
    do j=1,ncell_ghost
      e = cell_ghost(j)
      etype = cell_type(e)
      do i=1,node_per_element(etype)
        write(unitnum) node_g2l(cell_nodes(i,e))
      end do
    end do

  

    ! write out bdry triangles
    do j=1,nbdry_face
      e = bdry_face_tag(j)
      etype = face_type(e)
      if(etype == 2) then
        do i=1,node_per_element(etype)
          write(unitnum) node_g2l(face_nodes(i,e))
        end do
      end if
    end do
    
    ! write out bdry quads
    do j=1,nbdry_face
      e = bdry_face_tag(j)
      etype = face_type(e)
      if(etype == 3) then
        do i=1,node_per_element(etype)
          write(unitnum) node_g2l(face_nodes(i,e))
        end do
      end if
    end do
  
    ! write out bdry triangle patch
    do j=1,nbdry_face
      e = bdry_face_tag(j)
      if(etype == 2) then
      write(unitnum) face_patch(e)
      end if
    end do
  
    ! write out bdry quad patch
    do j=1,nbdry_face
      e = bdry_face_tag(j)
      if(etype == 3) then
        write(unitnum) face_patch(e)
      end if
    end do
  
    ! print out all nodes real and then ghost
    do j=1,nnodes_real+nnodes_ghost
      do i=1,3
        write(unitnum) xgeom(i,node_l2g(j))
      end do
    end do
    
    call mpi_schedule_out_reverse(mpi_cell,unitnum)
    call mpi_schedule_out_reverse(mpi_cell,unitnum)

    if(ntris > 0) then
      write(unitnum) ((dummy,j=1,5),i=1,ntris)
    end if
    if(nquads > 0) then
      write(unitnum) ((dummy,j=1,6),i=1,nquads)
    end if    


    ! print out global index for each element type separately
    ! need to subtract off global amounts
    e = 0
    do j=1,ncell_type_real(4)
      e = e+1
      write(unitnum) cell_real(e)
    end do
  
    do j=1,ncell_type_real(5)
      e = e+1
      write(unitnum) cell_real(e)-ntet
    end do
  
    do j=1,ncell_type_real(6)
      e = e+1
      write(unitnum) cell_real(e)-ntet-npyramid
    end do
  
    do j=1,ncell_type_real(7)
      e = e+1
      write(unitnum) cell_real(e)-ntet-npyramid-nprism
    end do
  
    close(unitnum)
  
    
    if(pdegree == 1) then
      return
    end if
  
    if(pdegree > 2) then
      print*,'ahh not tested for high p, aborting'
      stop
    end if
  
    node_g2l(:) = 0	
    node_l2g(:) = 0
    nnodes_real = 0
  
    do i=1,ncell_real
      e = cell_real(i)
      etype = cell_type(e)
    
      select case(etype)
        case(4)
          tm = (pdegree+1)*(pdegree+2)*(pdegree+3)/6
        case(5)
          print*,'ahh not tested aborting'
          stop
          tm = (pdegree+1)*(pdegree+2)*(2*pdegree+3)/6
        case(6)
          print*,'ahh not tested aborting'
          stop
          tm = (pdegree+1)*(pdegree+1)*(pdegree+2)/2
        case(7)
          print*,'ahh not tested aborting'
          stop
          tm = (pdegree+1)*(pdegree+1)*(pdegree+1)
      end select	

      do j=1,tm
        if(node_g2l(cell_nodes(j,e)) == 0) then
          nnodes_real = nnodes_real+1
          node_g2l(cell_nodes(j,e)) = nnodes_real
          node_l2g(nnodes_real) = cell_nodes(j,e)
        end if
      end do
          
    end do
  
  
    filename = 'part.'//trim(adjustl(nparts_string))//'/curved_cell.'//trim(adjustl(part_string))

    open (unit=unitnum,file=filename, form='unformatted',convert='big_endian')

    write(unitnum) pdegree,nnodes_real,ncell_type_real(4),ncell_type_real(5),ncell_type_real(6),ncell_type_real(7)

    write(unitnum) (xgeom(1,node_l2g(i)),i=1,nnodes_real)
    write(unitnum) (xgeom(2,node_l2g(i)),i=1,nnodes_real)
    write(unitnum) (xgeom(3,node_l2g(i)),i=1,nnodes_real)
    
    ! always contains p1 vertex 
    mode_tetra(1,1) = -1.0_dp
    mode_tetra(1,2) =  1.0_dp
    mode_tetra(1,3) = -1.0_dp
    mode_tetra(1,4) = -1.0_dp
    
    mode_tetra(2,1) = -1.0_dp
    mode_tetra(2,2) = -1.0_dp
    mode_tetra(2,3) =  1.0_dp
    mode_tetra(2,4) = -1.0_dp
    
    mode_tetra(3,1) = -1.0_dp
    mode_tetra(3,2) = -1.0_dp
    mode_tetra(3,3) = -1.0_dp
    mode_tetra(3,4) =  1.0_dp
    
    
    ! p=2 middle of edges
    mode_tetra(1,5) =  0.0_dp
    mode_tetra(1,6) =  0.0_dp
    mode_tetra(1,7) = -1.0_dp
    mode_tetra(1,8) = -1.0_dp
    mode_tetra(1,9) = -1.0_dp
    mode_tetra(1,10) = 0.0_dp
  
    mode_tetra(2,5) = -1.0_dp
    mode_tetra(2,6) =  0.0_dp
    mode_tetra(2,7) =  0.0_dp
    mode_tetra(2,8) = -1.0_dp
    mode_tetra(2,9) =  0.0_dp
    mode_tetra(2,10) =-1.0_dp
   
    mode_tetra(3,5) = -1.0_dp
    mode_tetra(3,6) = -1.0_dp
    mode_tetra(3,7) = -1.0_dp
    mode_tetra(3,8) =  0.0_dp
    mode_tetra(3,9) =  0.0_dp
    mode_tetra(3,10) = 0.0_dp
    
    if(ncell_type_real(4) > 0) then
      tm = (pdegree+1)*(pdegree+2)*(pdegree+3)/6
      write(unitnum) (mode_tetra(1,j),j=1,tm)     !xi
      write(unitnum) (mode_tetra(2,j),j=1,tm)     !eta
      write(unitnum) (mode_tetra(3,j),j=1,tm)     !zeta			
      write(unitnum) ((node_g2l(cell_nodes(j,cell_real(i))),j=1,tm),i=1,ncell_type_real(4))
    end if   

    if(ncell_type_real(5) > 0) then
      print*,'ahh not working yet abort'
      stop
      tm  = (pdegree+1)*(pdegree+2)*(2*pdegree+3)/6
      write(unitnum) (mode_pyr(1,j),j=1,tm)     !xi
      write(unitnum) (mode_pyr(2,j),j=1,tm)     !eta
      write(unitnum) (mode_pyr(3,j),j=1,tm)     !zeta			
      write(unitnum) ((node_g2l(cell_nodes(j,cell_real(i+ncell_type_real(4)))),j=1,tm),i=1,ncell_type_real(5))
    end if

    if(ncell_type_real(6) > 0) then
      print*,'ahh not working yet abort'
      stop
      tm = (pdegree+1)*(pdegree+1)*(pdegree+2)/2
      write(unitnum) (mode_prism(1,j),j=1,tm)     !xi
      write(unitnum) (mode_prism(2,j),j=1,tm)     !eta
      write(unitnum) (mode_prism(3,j),j=1,tm)     !zeta			
      write(unitnum) ((node_g2l(cell_nodes(j,cell_real(i+ncell_type_real(4)+ncell_type_real(5)))),j=1,tm),i=1,ncell_type_real(6))
    end if

    if(ncell_type_real(7) > 0) then
      print*,'ahh not working yet abort'
      stop
      tm = (pdegree+1)*(pdegree+1)*(pdegree+1)
      write(unitnum) (mode_hex(1,j),j=1,tm)     !xi
      write(unitnum) (mode_hex(2,j),j=1,tm)     !eta
      write(unitnum) (mode_hex(3,j),j=1,tm)     !zeta			
      write(unitnum) ((node_g2l(cell_nodes(j,cell_real(i+ncell_type_real(4)+ncell_type_real(5)+ncell_type_real(6)))),j=1,tm),i=1,ncell_type_real(7))
    end if

          
    close(unitnum)
  

  
  
  end subroutine output_partition_mesh_file


	
  subroutine output_mesh()
    
    use my_kinddefs
    use mesh_module
    use connection_module

    implicit none
    
    !character(len=255), intent(in) :: output_mesh_filename
    character(len=255) :: output_mesh_filename,string_id_proc       
    integer(i4) :: i,j,offset,unitnum
    integer(i4) :: vtu_tags(7) = [3,5,9,10,14,13,12]
    
    output_mesh_filename = 'mesh.vtu'
     
    unitnum = 33
        
    open (unit=unitnum, file=output_mesh_filename, status='unknown')

    write(unitnum,*) '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
    write(unitnum,*) '<UnstructuredGrid>'
    write(unitnum,*) '<Piece NumberOfPoints="',nnodes,'" NumberOfCells="',ncell, '" >'
   
    write(unitnum,*) '<PointData>'
    write(unitnum,*) '</PointData>'
   
    write(unitnum,*) '<CellData>'
    write(unitnum,*) '<DataArray type="Int32" Name="partition" format="ascii">'
    do i=1,ncell
      write(unitnum,*) cell_partition(i)
    end do

    write(unitnum,*) '</DataArray>'

    write(unitnum,*) '</CellData>'

    
    
    write(unitnum,*) '<Points>'
      write(unitnum,*) '<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
      do i = 1, nnodes
        write(unitnum,*) xgeom(1,i),xgeom(2,i),xgeom(3,i)
      end do
      write(unitnum,*) '</DataArray>'
    write(unitnum,*) '</Points>'
   
    write(unitnum,*) '<Cells>'
    write(unitnum,*) '<DataArray type="Int32" Name="connectivity" format="ascii">'
    do i=1,ncell
      write(unitnum,*) ( cell_nodes(j,i)-1, j=1,node_per_element(cell_type(i)) )
    end do

    write(unitnum,*) '</DataArray>'
    
    write(unitnum,*) '<DataArray type="Int32" Name="offsets" format="ascii">'
    offset = 0
    do i=1,ncell
      offset = offset + node_per_element(cell_type(i))
      write(unitnum,*) offset
    end do
    write(unitnum,*) '</DataArray>'
    
    write(unitnum,*) '<DataArray type="UInt8" Name="types" format="ascii">'
    
    do i=1,ncell
      write(unitnum,*) vtu_tags(cell_type(i)) 
    end do
  
    write(unitnum,*) '</DataArray>'
                        
    write(unitnum,*) '</Cells>'
      
    write(unitnum,*) '</Piece>'
    write(unitnum,*) '</UnstructuredGrid>'
    write(unitnum,*) '</VTKFile>'
    close(unitnum)



    
    
end subroutine output_mesh


end module partition_module



