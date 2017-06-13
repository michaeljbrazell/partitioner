subroutine load_mesh(filename,mesh_type)

  use my_kinddefs

  implicit none

  character(len=255), intent(in) :: filename
  integer(i4),intent(in) :: mesh_type


  select case(mesh_type)
    case(1)
      call read_mcell_mesh_header(filename)
      call allocate_mesh()
      call read_mcell_mesh(filename)
    case(2)
      call read_ugrid_mesh_header(filename)
      call allocate_mesh()
      call read_ugrid_mesh(filename)
    case(3)
      call read_gmsh_header(filename)
      call allocate_mesh()
      call read_gmsh(filename)
    case default
      print*,'mesh type unknown'
      stop
    
  end select

end subroutine load_mesh


subroutine read_gmsh_header(filename)

  use my_kinddefs
  use mesh_module
    
  implicit none
  
  character(len=255), intent(in) :: filename


  ! obsolete stuff
  integer(i4) :: ncomponent,nbody,nfrpts,dummy
  integer(i4) :: i,j,inmesh,nel,etype,ninfo
  integer(i4) :: info(100)
  inmesh = 34

    ! open mesh .unf mesh file
    open (unit=inmesh,file=filename, form='formatted')

  read(inmesh,*) 
  read(inmesh,*) 
  read(inmesh,*) 
  read(inmesh,*) 

  ntet=0
  npyramid=0
  nprism=0
  nhex=0
  ntri_bdry=0
  nquad_bdry=0
  npatch=0
  
  read(inmesh,*) nnodes
  
  do i=1,nnodes
    read(inmesh,*)
  end do
  
  read(inmesh,*)
  read(inmesh,*)
  read(inmesh,*) nel
  
  do i=1,nel
  
    read(inmesh,*) dummy, etype, ninfo, (info(j),j=1,ninfo)
    
    select case(etype)
      case(2,9,21,23)
        ntri_bdry = ntri_bdry+1
        npatch = max(npatch,info(1))
      case(3,10)
        nquad_bdry = nquad_bdry+1
        npatch = max(npatch,info(1))
      case(4,11,29,30,31)
        ntet=ntet+1
      case(5,12,92,93)
        nhex=nhex+1
      case(6,13)
        nprism=nprism+1
      case(7,14)
        npyramid=npyramid+1
    end select


  end do


  print*,'ntet',ntet
  print*,'npyramid',npyramid
  print*,'nprism',nprism
  print*,'nhex',nhex

  print*,'nnodes',nnodes
  print*,'ntri_bdry',ntri_bdry
  print*,'nquad_bdry',nquad_bdry
  print*,'npatch',npatch

  close(inmesh)


end subroutine read_gmsh_header


subroutine read_mcell_mesh_header(filename)

  use my_kinddefs
  use mesh_module
  
  implicit none
  
  character(len=255), intent(in) :: filename


  ! obsolete stuff
  integer(i4) :: ncomponent,nbody,nfrpts,dummy
  integer(i4) :: i,j,inmesh
  
  inmesh = 34

  ! open mesh .unf mesh file
  open (unit=inmesh,file=filename, form='unformatted',convert='big_endian')

  ! read in number of tet, pyramid, prism, and hex
  read(inmesh) ntet,npyramid,nprism,nhex

  print*,'ntet',ntet
  print*,'npyramid',npyramid
  print*,'nprism',nprism
  print*,'nhex',nhex

  ! read in number of nodes, boundary tris, boundary quads, and other stuff
  read(inmesh) nnodes,ntri_bdry,nquad_bdry,npatch,ncomponent,nbody,nfrpts

  print*,'nnodes',nnodes
  print*,'ntri_bdry',ntri_bdry
  print*,'nquad_bdry',nquad_bdry
  print*,'npatch',npatch
    
  close(inmesh)

end subroutine read_mcell_mesh_header



subroutine read_ugrid_mesh_header(filename)

  use my_kinddefs
  use mesh_module
  
  implicit none
  
  character(len=255), intent(in) :: filename


  ! obsolete stuff
  integer(i4) :: ncomponent,nbody,nfrpts,dummy
  integer(i4) :: i,j,inmesh,os
  
  inmesh = 34

  ! open mesh ugrid mesh file
  open (unit=inmesh,file=filename, form='unformatted',access="stream", &
                                   status='unknown', iostat=os)! , &
                                   !convert='little_endian')



  read(inmesh) nnodes,ntri_bdry,nquad_bdry,ntet,npyramid,nprism,nhex

  print*,'ntet',ntet
  print*,'npyramid',npyramid
  print*,'nprism',nprism
  print*,'nhex',nhex

  npatch = 1

  print*,'nnodes',nnodes
  print*,'ntri_bdry',ntri_bdry
  print*,'nquad_bdry',nquad_bdry
  ! 	print*,'npatch',npatch
    
  close(inmesh)

end subroutine read_ugrid_mesh_header



subroutine allocate_mesh

  use my_kinddefs
  use mesh_module

  implicit none

  ncell = ntet+npyramid+nprism+nhex

  pdegree = 1

  allocate(cell_nodes(125,ncell))
  cell_nodes = -1
  allocate(cell_type(ncell))
  cell_type = -1
  allocate(cell2face(6,ncell))
  cell2face = -1
  allocate(cell2cell(6,ncell))
  cell2cell = -1	
  allocate(cell_g2l(ncell))
  cell_g2l = -1	
  allocate(cell_partition(ncell))
  cell_partition(:) = 1
  
  nface_bdry = ntri_bdry + nquad_bdry

  ! total number of triangles in mesh
    ntri = (4*ntet + 4*npyramid + 2*nprism + ntri_bdry)/2
    ! total number of quads in mesh
    nquad = (npyramid + 3*nprism + 6*nhex + nquad_bdry)/2

  ! total number of faces
  nface = ntri+nquad

  allocate(face_nodes(4,nface))
  face_nodes = -1
  allocate(face_type(nface))
  face_type = -1
  allocate(face_patch(nface))
  face_patch = -1	
  allocate(face2cell(2,nface))
  face2cell = -1

  print*,'ntri',ntri
  print*,'nquad',nquad

  allocate(xgeom(3,nnodes))

end subroutine allocate_mesh


subroutine read_gmsh(filename)

  use my_kinddefs
  use mesh_module
    
  implicit none
  
  character(len=255), intent(in) :: filename


  ! obsolete stuff
  integer(i4) :: ncomponent,nbody,nfrpts,dummy
  integer(i4) :: i,j,k,etype,inmesh,nel,ninfo
  integer(i4) :: vert(100),nvert(100),info(100)
  integer(i4) :: itet,ipyr,ipri,ihex,itri,iqua
  
  inmesh = 34

  ! open mesh .unf mesh file
  open (unit=inmesh,file=filename, form='formatted')

  read(inmesh,*)
  read(inmesh,*)
  read(inmesh,*)
  read(inmesh,*)
  read(inmesh,*)
  
  do i=1,nnodes
    read(inmesh,*) dummy, (xgeom(j,i),j=1,3)
  end do

  read(inmesh,*)
  read(inmesh,*)
  read(inmesh,*) nel

  nvert(:) = 0
  
  ! tri and quad
  nvert(2) = 3
  nvert(3) = 4
  
  ! high order tris/quads but only load 3,4
  nvert(9) = 3
  nvert(21) = 3
  nvert(23) = 3		
  nvert(25) = 3		
  nvert(10) = 4
  
  ! p1 cells
  nvert(4) = 4
  nvert(5) = 8
  nvert(6) = 6
  nvert(7) = 5
  
  ! high order cells
  nvert(11) = 10
  nvert(12) = 27
  nvert(13) = 18
  nvert(14) = 14
  nvert(29) = 20
  nvert(30) = 35
  nvert(31) = 56
  nvert(92) = 64
  nvert(93) = 125


  itet=0
  ipyr=0
  ipri=0
  ihex=0
  itri=0
  iqua=0
            
  do i=1,nel
  
    read(inmesh,*) dummy, etype, ninfo, (info(j),j=1,ninfo), (vert(k),k=1,nvert(etype))
    
    select case(etype)
      case(8,9,10,11,12,13,14)
        pdegree = 2
      case(26,21,29,92)
        pdegree = 3
      case(27,23,30,93)
        pdegree = 4
      case(28,25,31)
        pdegree = 5
    end select
    
    select case(etype)
    
      case(1,8,26,27,28)
        ! this is an edge don't need it for DG
        
      case(2,9,21,23,25)
        itri = itri + 1
        do j=1,nvert(etype) 
          face_nodes(j,itri) = vert(j)
        end do
        face_type(itri) = 2
        face_patch(itri) = info(1)
        
      case(3,10) 
        iqua = iqua + 1
        do j=1,nvert(etype)
          face_nodes(j,iqua+ntri_bdry) = vert(j)
        end do
        face_type(iqua+ntri_bdry) = 3
        face_patch(iqua+ntri_bdry) = info(1)
        
      case(4,11,29,30,31)
        itet = itet+1
        do j=1,nvert(etype)
          cell_nodes(j,itet) = vert(j)
        end do
        cell_type(itet) = 4

      case(5,12,92,93)
        ihex = ihex+1
        do j=1,nvert(etype)
          cell_nodes(j,ntet+npyramid+nprism+ihex) = vert(j)
        end do
        cell_type(ntet+npyramid+nprism+ihex) = 7
      
      case(6,13)
        ipri = ipri+1
        do j=1,nvert(etype)
          cell_nodes(j,ntet+npyramid+ipri) = vert(j)
        end do
        cell_type(ntet+npyramid+ipri) = 6
          
      case(7,14)
        ipyr = ipyr+1
        do j=1,nvert(etype)
          cell_nodes(j,ntet+ipyr) = vert(j)
        end do
        cell_type(ntet+ipyr) = 5
        
      case(15)
        ! this is a point don't need it for DG

      case default
        print*,'ahh dont have this element type figured out yet'
        stop
        
    end select

    
  end do
    
  close(inmesh)


end subroutine read_gmsh





subroutine read_mcell_mesh(filename)

  use my_kinddefs
  use mesh_module
    
  implicit none
  
  character(len=255), intent(in) :: filename


  ! obsolete stuff
  integer(i4) :: ncomponent,nbody,nfrpts,dummy
  integer(i4) :: i,j,inmesh
  
  inmesh = 34

  ! open mesh .unf mesh file
  open (unit=inmesh,file=filename, form='unformatted',convert='big_endian')

  ! read in number of tet, pyramid, prism, and hex
  read(inmesh) ntet,npyramid,nprism,nhex

  ! read in number of nodes, boundary tris, boundary quads, and other stuff
  read(inmesh) nnodes,ntri_bdry,nquad_bdry,npatch,ncomponent,nbody,nfrpts


  ! load tet nodes
  if (ntet > 0) then
    do i=1,4        
      read(inmesh)(cell_nodes(i,j), j=1, ntet)
    end do
    do j=1,ntet
      cell_type(j) = 4
    end do
  end if

  ! load pyramid nodes
  if (npyramid > 0) then
    do i=1,5
      read(inmesh)(cell_nodes(i,ntet+j), j=1, npyramid)
    end do
    do j=1,npyramid
      cell_type(ntet+j) = 5
    end do
  end if

  ! load prism nodes
  if (nprism > 0) then
    do i=1,6
      read(inmesh)(cell_nodes(i,ntet+npyramid+j), j=1, nprism)
    end do
    do j=1,nprism
      cell_type(ntet+npyramid+j) = 6
    end do
  end if

    
  ! load hex nodes
  if (nhex > 0) then
    do i=1,8
      read(inmesh)(cell_nodes(i,ntet+npyramid+nprism+j), j=1, nhex)
    end do        
    do j=1,nhex
      cell_type(ntet+npyramid+nprism+j) = 7
    end do
  end if
    
  ! load x,y,z location of nodes
  do i=1,3
    read(inmesh)(xgeom(i,j),j=1,nnodes)
  end do

  ! load triangle nodes
  if (ntri_bdry > 0) then
    do i=1,3
      read(inmesh)(face_nodes(i,j), j=1, ntri_bdry)
    end do
    read(inmesh)(dummy, j=1,ntri_bdry)
    do j=1,ntri_bdry
      face_type(j) = 2
    end do
  end if

! load quad nodes
  if (nquad_bdry > 0) then
    do i=1,4
      read(inmesh)(face_nodes(i,j+ntri_bdry),j=1,nquad_bdry)
    end do
    read(inmesh)(dummy, j=1,nquad_bdry)
    do j=1,nquad_bdry
      face_type(ntri_bdry+j) = 3
    end do
  end if

   
  ! load tri bdry patch number
  if (ntri_bdry > 0) then
    read(inmesh)(face_patch(j), j=1, ntri_bdry)
    do i=1,6
      read(inmesh)(dummy, j=1, ntri_bdry)
    end do
  end if

  ! load quad patch number
  if (nquad_bdry > 0) then
    read(inmesh)(face_patch(j+ntri_bdry), j=1, nquad_bdry)
    do i=1,8
      read(inmesh)(dummy, j=1, nquad_bdry)
    end do
  end if

    close(inmesh)


end subroutine read_mcell_mesh



subroutine read_ugrid_mesh(filename)

  use my_kinddefs
  use mesh_module
    
  implicit none
  
  character(len=255), intent(in) :: filename


  integer(i4) :: i,j,inmesh,os,dummy
  
  inmesh = 34

  ! open mesh ugrid mesh file
  open (unit=inmesh,file=filename, form='unformatted',access="stream", &
                                   status='unknown', iostat=os) !, &
                                   !convert='little_endian')
                                     
  read(inmesh) nnodes,ntri_bdry,nquad_bdry,ntet,npyramid,nprism,nhex

  ! load x,y,z location of nodes
  read(inmesh) ((xgeom(i,j),i=1,3),j=1,nnodes)
 

  ! load triangle nodes
    if (ntri_bdry > 0) then
      read(inmesh)((face_nodes(i,j),i=1,3), j=1, ntri_bdry)
      do j=1,ntri_bdry
        face_type(j) = 2
      end do
    end if

  ! load quad nodes
    if (nquad_bdry > 0) then
      read(inmesh)((face_nodes(i,j+ntri_bdry),i=1,4),j=1,nquad_bdry)
      do j=1,nquad_bdry
        face_type(ntri_bdry+j) = 3
      end do
    end if
   
    ! load tri bdry patch number
    if (ntri_bdry > 0) then
       read(inmesh)(face_patch(j), j=1, ntri_bdry)
    end if

    ! load quad patch number
    if (nquad_bdry > 0) then
      read(inmesh)(face_patch(j+ntri_bdry), j=1, nquad_bdry)
    end if
    
    
  ! load tet nodes
  if (ntet > 0) then
    read(inmesh)(cell_nodes(1,j),cell_nodes(2,j),cell_nodes(3,j),cell_nodes(4,j), j=1, ntet)
    do j=1,ntet
      cell_type(j) = 4
    end do
  end if

  ! load pyramid nodes
  if (npyramid > 0) then
    read(inmesh)(cell_nodes(1,ntet+j),cell_nodes(2,ntet+j), &
                cell_nodes(3,ntet+j),cell_nodes(4,ntet+j), &
                cell_nodes(5,ntet+j), j=1, npyramid)
    do j=1,npyramid
      cell_type(ntet+j) = 5
    end do
  end if

  ! load prism nodes
  if (nprism > 0) then
    read(inmesh)(cell_nodes(1,ntet+npyramid+j),cell_nodes(2,ntet+npyramid+j), &
                 cell_nodes(3,ntet+npyramid+j),cell_nodes(4,ntet+npyramid+j), &
                 cell_nodes(5,ntet+npyramid+j),cell_nodes(6,ntet+npyramid+j), j=1, nprism)
    do j=1,nprism
      cell_type(ntet+npyramid+j) = 6
    end do
  end if

    
  ! load hex nodes
  if (nhex > 0) then
    read(inmesh)(cell_nodes(1,ntet+npyramid+nprism+j),cell_nodes(2,ntet+npyramid+nprism+j), &
                 cell_nodes(3,ntet+npyramid+nprism+j),cell_nodes(4,ntet+npyramid+nprism+j), &
                 cell_nodes(5,ntet+npyramid+nprism+j),cell_nodes(6,ntet+npyramid+nprism+j), &
                 cell_nodes(7,ntet+npyramid+nprism+j),cell_nodes(8,ntet+npyramid+nprism+j), j=1, nhex)
    do j=1,nhex
      cell_type(ntet+npyramid+nprism+j) = 7
    end do
  end if
    


    close(inmesh)


end subroutine read_ugrid_mesh





