!***************************************************
!                 module my_kinddefs:
!     Specifies machine independant precision of 
!     double and single precision real variables 
!      and 4 types of integer variables.  The reals
!      are mad acessible through variables sp and dp 
!      while the integers are i1, i2 i4, i8, each of 
!      which has that many digits behind the decimal place
!**********************************************************      

module my_kinddefs

  implicit none

  public sp, dp, cp, qp, i1, i2, i4, i8
  public PI, TWOPI, PIO2, half, onethird, twothird, threeights, fourthird

  private
  !---> All these a based on the default IEEE 754 format 

  !---> Single Precision Real: IEEE 754
  integer, parameter :: sp=selected_real_kind(6,37)
  
  !---> Double Precision Real: IEEE 754
  integer, parameter :: dp=selected_real_kind(15,307)
  !integer, parameter :: dp=selected_real_kind(33,4931)
  
  !---> Double Precision Complex: IEEE 754
  integer, parameter :: cp=selected_real_kind(15,307)
  !integer, parameter :: cp=selected_real_kind(33,4931)
  
  !---> Quad Precision Real: IEEE 754
  integer, parameter :: qp=selected_real_kind(33, 4931)
  
  !---> Single decimal place integer (1 byte)
  integer, parameter :: i1=selected_int_kind(2)
  
  !---> 2 decimal place integer (2 bytes): 
  integer, parameter :: i2=selected_int_kind(3)

  !---> 4 decimal place integer ( 4 bytes):
  integer, parameter :: i4=selected_int_kind(5)

  !---> 8 decimal place integer (8 bytes):
  integer, parameter :: i8=selected_int_kind(10)

  !---> Additionally we put some very handy constants in here
  !     to make it a bit more like matlab and other tools
  real(dp), parameter :: half = 1._dp/2._dp
  real(dp), parameter :: onethird= 1._dp/3._dp
  real(dp), parameter :: twothird= 2._dp/3._dp
  real(dp), parameter :: threeights= 3._dp/8._dp
  real(dp), parameter :: fourthird = 4._dp/3._dp
  
  real(dp), parameter :: PI=3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: PIO2=1.57079632679489661923132169163975144209858_dp
  real(dp), parameter :: TWOPI=6.283185307179586476925286766559005768394_dp


end module my_kinddefs
