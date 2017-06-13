module line_module

  use my_kinddefs
  use mesh_module
      
  implicit none
  
  type line_struct

    integer(i4) :: num_cells      
    integer(i4), allocatable :: cell(:)

  end type line_struct

  integer(i4) :: num_lines
      
  type(line_struct), allocatable :: line(:)  
 
  integer(i4), parameter :: tag_boundary = 1
  integer(i4), parameter :: max_line_length = 10000000
  
  contains
  
  subroutine setup_lines
    call fill_lines(0)
    call fill_lines(1)
    call fill_lines(2)
  end subroutine setup_lines
    
 
  ! if line flag is 0 then allocate lines
  ! if line flag is 1 then allocate cells
  ! if line flag is 2 then add cells    
  subroutine fill_lines(line_flag)
  
    use connection_module
    
    integer(i4), intent(in) :: line_flag    	
    integer(i4) :: i,f,f2,e,e2,etype,etype2,num_cells
    integer(i4) :: cell_tag(ncell)
    
    cell_tag(:) = 0
    
    num_lines = 0
    
    do f=1,nface_bdry
!     if(face_patch(f) < 172 .or. face_patch(f) > 177) then
!     if(face_patch(f) < 4 ) then
      if(face_patch(f) == tag_boundary) then
!     if(face_patch(f) == 1) then
!     if(face_patch(f) < 4) then
    ! find first cell on boundary face
        e = face2cell(1,f)
        etype = cell_type(e)

        ! if a prism and triangle face
        if(cell_type(e) == 6 .and. face_type(f) == 2) then
          
          num_lines = num_lines + 1
          num_cells = 0
          
          if(line_flag > 0) then
            line(num_lines)%num_cells = 0
          end if
          
          f2 = f
          e2 = e
          etype2 = etype
          
          do while(etype2 == 6) 
          
            if(num_cells > max_line_length) then
              exit
            end if
            num_cells = num_cells + 1
          
            if(line_flag > 0) then    
              line(num_lines)%num_cells = line(num_lines)%num_cells + 1
            end if
            
            cell_tag(e2) = cell_tag(e2) + 1
          
            if(line_flag == 2) then
              line(num_lines)%cell(line(num_lines)%num_cells) = e2
            end if

            if(cell2face(4,e2) == f2) then
              f2 = cell2face(5,e2)
              e2 = cell2cell(5,e2)
              etype2 = cell_type(e2)
            else if(cell2face(5,e2) == f2) then
                f2 = cell2face(4,e2)
                e2 = cell2cell(4,e2)
              etype2 = cell_type(e2)
            end if
            
          end do
        end if
      
        ! if hex and quad face
        if(cell_type(e) == 7 .and. face_type(f) == 3) then
      
          num_lines = num_lines + 1
          num_cells = 0

          if(line_flag > 0) then
            line(num_lines)%num_cells = 0
          end if
        
          f2 = f
          e2 = e
          etype2 = etype
          
          do while(etype2 == 7)

            if(num_cells > max_line_length) then
              exit
            end if
            num_cells = num_cells + 1
            
            if(line_flag > 0) then					
              line(num_lines)%num_cells = line(num_lines)%num_cells + 1
            end if
            
            cell_tag(e2) = cell_tag(e2) + 1

            if(line_flag == 2) then
              line(num_lines)%cell(line(num_lines)%num_cells) = e2
            end if
            
            do i=1,6
              if(cell2face(i,e2) == f2) then
                f2 = cell2face(opposite_face(7,i),e2)
                e2 = cell2cell(opposite_face(7,i),e2)
                if(e2 > 1) then
                  etype2 = cell_type(e2)
                else
                  etype2 = -1
                end if
                exit
              end if
            end do
            
          end do
        end if
      end if
    end do
    
    
      ! find all cells not tagged and make a length 1 line
    do e=1,ncell
    
      if(cell_tag(e) == 0) then
        
        num_lines = num_lines+1
        
        ! allocate one cell 
        if(line_flag > 0) then
          line(num_lines)%num_cells = 1
        end if   
        ! fill in the first cell
        if(line_flag == 2) then
          line(num_lines)%cell(1) = e
        end if
      end if
    
      if(cell_tag(e) > 1) then
        print*,'cell contains more than one line',e,cell_tag(e)
      end if
    
    end do
    
    ! allocate all lines
    if(line_flag == 0) then		
      allocate(line(num_lines))
    end if
  
    ! allocate cells on lines
    if(line_flag == 1) then
      do i=1,num_lines
        allocate(line(i)%cell(line(i)%num_cells))
      end do
    end if
    
    
  end subroutine fill_lines
    
    
    
    
    ! dont need anymore
!      subroutine count_num_lines()
!   
!     integer(i4) :: f,e
!     
!     num_lines = 0
!     
!     do f=1,nface_bdry
!       if(face_patch(f) == tag_boundary) then
!         ! find first cell on boundary face
!         e = face2cell(1,f)
!       
!         if(e < 1) then
!           print*,'error e is less than 1'
!           stop
!       end if
!       
!         ! if a prism and triangle face
!         if(cell_type(e) == 6 .and. face_type(f) == 2) then
!           num_lines = num_lines + 1
!       end if
!       
!       ! if hex and quad face
!       if(cell_type(e) == 7 .and. face_type(f) == 3) then
!           num_lines = num_lines + 1
!       end if
!       
!     end if
!     end do
!     
!     
!     allocate(line(num_lines))
!     
!   end subroutine count_num_lines
   
   
  subroutine deallocate_lines()
   
    integer(i4) :: i

    do i=1,num_lines
      deallocate(line(i)%cell)
    end do

    deallocate(line)
    
  end subroutine deallocate_lines
   
end module line_module
