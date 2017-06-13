subroutine setup_mesh()

	call setup_face()
	call setup_cell2cell()
	
end subroutine setup_mesh


subroutine setup_face()
    
    use my_kinddefs
    use mesh_module
    use connection_module
    implicit none
    
    integer(i4) :: i,j,k,m,lcl,ind(4),minv,find,eind,findprev,nf,ftype,ftype2,etype,istat
    integer(i4) :: v(4), a(4),vstore(4),astore(4)
   
    
    integer(i4) :: minv_list(nnodes)
    integer(i4) :: face_list(nface)
    
    ! number of faces equals only boundary faces for now
    nf = nface_bdry
 
    ! initialize linked lists
    minv_list(:) = -1
    face_list(:) = -1

    ! fill in face data for boundary faces
    do i=1,nface_bdry
        
        if(face_type(i) == 2) then
            do j=1,3
                v(j) = face_nodes(j,i)
            end do
            minv = min(v(1),v(2),v(3))
        else
            do j=1,4
                v(j) = face_nodes(j,i)
            end do
            minv = min(v(1),v(2),v(3),v(4))
        end if
        
        find = minv_list(minv)
        
        do while (find > 0)
            findprev = find
            find = face_list(find)
        end do
        
        if(minv_list(minv) < 0) then
            minv_list(minv) = i
        else
            face_list(findprev) = i
        end if
        
        npatch = max(npatch,face_patch(i))
        
    end do

     
    ! loop through cells and find all faces
    do i=1,ncell

        ! element type
        etype = cell_type(i)

        ! loop through faces and find if it exists or not
        do j=1,face_per_element(etype)      

            ! face type triangle or quad
            ftype = element_type_face(etype,j)

            ! vertex for face
            do k=1,node_per_element(ftype)
                v(k) = cell_nodes(face_vertex(etype,j,k),i)
            end do     
            
            vstore = v

            ! sort vertex from small to big
            call insort(node_per_element(ftype),v(1:node_per_element(ftype)))

            ! min vertex number
            minv = v(1)

            ! find face to start search
            find = minv_list(minv)
       
            ! if find not -1 then search for face
            do while(find > 0)

                ftype2 = face_type(find)
                
                if(ftype == ftype2) then

                    do k=1,node_per_element(ftype2)
                        a(k) = face_nodes(k,find)
                    end do
                    astore = a

                    ! sort vertex from small to big
                    call insort(node_per_element(ftype2),a)

                    lcl = 0
                    do k=1,node_per_element(ftype2)
                        lcl = lcl + abs(a(k)-v(k))
                    end do

                    ! match found
                    if(lcl == 0) then
                        
                        if(face2cell(1,find) < 0) then
                            face2cell(1,find) = i
                            cell2face(j,i) = find
                        else
                            ! interior face found second time
                            face2cell(2,find) = i 
                            cell2face(j,i) = find 
                        end if

                        goto 123
                    end if

                end if

                findprev = find
                find = face_list(find)                              

            end do

            ! interior face not found create new face
            nf = nf+1
            ! specify face type
            face_type(nf) = ftype
            ! left cell connected to face
            face2cell(1,nf) = i
            cell2face(j,i) = nf

            ! fill in nodes for face nf            
            do k=1,node_per_element(ftype)
                face_nodes(k,nf) = vstore(k)
            end do

            ! update linked lists
            if(minv_list(minv) < 0) then
                minv_list(minv) = nf
            else
                face_list(findprev) = nf
            end if
           
        123 end do
 
        end do
       
       
        if(nf /= nface) then
            print*,'number of faces does not match',nf,nface
            stop
        end if

     
    
end subroutine


subroutine setup_cell2cell()

  use my_kinddefs
  use mesh_module
  use connection_module
  implicit none
  
  integer(i4) :: i,j,etype,e1,e2,find
      
    ! fill in the cell2cell data structure    	
      do i=1,ncell
        etype = cell_type(i)
        do j=1,face_per_element(etype)
          find = cell2face(j,i) 
          e1 = face2cell(1,find)
          e2 = face2cell(2,find)
          if(i==e1) then
            cell2cell(j,i) = e2
          else if(i==e2) then
            cell2cell(j,i) = e1
          end if
        end do
    end do
    
end subroutine setup_cell2cell




subroutine insort(n,array)

  !..Insertion Sort for Array of n elements
  !..Warning: Method is n**2
  !..Only use with small n : ie n < 10

  use my_kinddefs
    
  implicit none

  !..dummy
  integer(i4), intent(IN)    :: n
  integer(i4), intent(INOUT) :: array(n)

  !..local
  integer(i4) :: i,j,a


  do j = 2,n

    a = array(j)

    do i=j-1, 1, -1

      if (array(i) <= a) go to 10
      array(i+1) = array(i)

    enddo

    i = 0
  10 array(i+1) = a
  enddo

end subroutine
