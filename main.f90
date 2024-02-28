program main
!========================================================================
!========================================================================
use atmosphere,  only: atmosphere_main
use fv_arrays,   only: fv_atmos_type

implicit none

type(fv_atmos_type) :: Atm

call atmosphere_main(Atm)


end program main
