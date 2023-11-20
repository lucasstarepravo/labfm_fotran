module common_parameter
  use kind_parameters
  implicit none 

  integer(ikind) ,parameter :: npar=140000
  !! NOTE: npmax should be npar*nplink
  integer(ikind) ,parameter :: dims = 2

  !! integer,parameter::nmax=500000

  real(rkind), parameter :: pi=3.141592653589793238462643383279502884197d0
  real(rkind), parameter :: sqrt2=dsqrt(2.0d0)
  real(rkind), parameter :: oosqrt2=1.0d0/dsqrt(2.0d0)
!!  real(rkind), parameter :: pi=4.0_rkind * DATAN(1.0_rkind)

end module common_parameter
