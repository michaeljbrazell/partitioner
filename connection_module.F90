!     
! File:   connection_module.F90
! Author: mbrazell
!
! Created on November 26, 2012, 2:22 PM
!

module connection_module

    use my_kinddefs
    
    ! DON'T TOUCH ANY OF THIS STUFF BELOW
       
    ! number of nodes for each type of element
    integer(i4) :: node_per_element(7) = [2,3,4,4,5,6,8]
    ! number of faces for each type of element
    integer(i4) :: face_per_element(7) = [0,1,1,4,5,5,6]


    ! type of face on each type of volume element
    integer(i4) :: element_type_face(4:7,6)
    data element_type_face(4,1:4)  / 2, 2, 2, 2 /
    data element_type_face(5,1:5)  / 3, 2, 2, 2, 2 /
    data element_type_face(6,1:5)  / 3, 3, 3, 2, 2 /
    data element_type_face(7,1:6)  / 3, 3, 3, 3, 3, 3 /

    ! the standard vertices for each face on each element
    ! first index is element type 
    ! second index is the face on element
    ! third index are the nodes on each face
    integer(i4) :: face_vertex(4:7,6,4)

    ! tetrahedra faces
    data face_vertex(4,1,1:3) / 3, 2, 4 /
    data face_vertex(4,2,1:3) / 1, 3, 4 /
    data face_vertex(4,3,1:3) / 1, 4, 2 /
    data face_vertex(4,4,1:3) / 1, 2, 3 /

    ! pyramid faces 
    data face_vertex(5,1,1:4) / 1, 2, 3,  4 /
    data face_vertex(5,2,1:3) / 5, 2, 1 /
    data face_vertex(5,3,1:3) / 5, 3, 2 /
    data face_vertex(5,4,1:3) / 5, 4, 3 /
    data face_vertex(5,5,1:3) / 5, 1, 4 /

    ! prism faces 
    data face_vertex(6,1,1:4) / 4, 5, 2,  1 /
    data face_vertex(6,2,1:4) / 3, 6, 4,  1 /
    data face_vertex(6,3,1:4) / 2, 5, 6,  3 /
    data face_vertex(6,4,1:3) / 1, 2, 3 /
    data face_vertex(6,5,1:3) / 6, 5, 4 /

    ! hex faces 
    data face_vertex(7,1,1:4) / 5, 8, 7, 6 /
    data face_vertex(7,2,1:4) / 1, 2, 3, 4 /
    data face_vertex(7,3,1:4) / 3, 7, 8, 4 /
    data face_vertex(7,4,1:4) / 2, 1, 5, 6 /
    data face_vertex(7,5,1:4) / 2, 6, 7, 3 /
    data face_vertex(7,6,1:4) / 1, 4, 8, 5 / 

    
    ! opposite faces 
    ! only works for hex and some prism faces
    integer(i4) :: opposite_face(4:7,1:6)
    data opposite_face(4,1:6) / -1,-1,-1,-1,-1,-1 /
    data opposite_face(5,1:6) / -1,-1,-1,-1,-1,-1 /
    data opposite_face(6,1:6) / -1,-1,-1, 5, 4,-1 /
    data opposite_face(7,1:6) /  2, 1, 4, 3, 6, 5 /
    

    
end module connection_module

