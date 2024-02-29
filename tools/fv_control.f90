module fv_control
!========================================================================
! Module for data allocation
!
! Reference
! https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere/blob/main/tools/fv_control.F90
!========================================================================
use fv_arrays,  only: fv_grid_bounds_type, fv_grid_type, fv_atmos_type, &
                      point_structure, R_GRID, erad, pi, pio4, pio2, nbfaces
use fv_grid_tools, only: init_grid
use fv_duogrid   , only: init_lagrange
implicit none

contains

!-------------------------------------------------
! define bounds
!-------------------------------------------------
subroutine init_bounds(bd, npx, npy)
   type(fv_grid_bounds_type), intent(INOUT) :: bd
   integer                  , intent(IN)    :: npx, npy

   bd%is  = 1
   bd%ie  = npx
   bd%isd = bd%is - bd%ng
   bd%ied = bd%ie + bd%ng

   bd%js  = 1
   bd%je  = npy
   bd%jsd = bd%js - bd%ng
   bd%jed = bd%je + bd%ng
 
   bd%npx = npx
   bd%npy = npy

   bd%isd1 = bd%isd-1
   bd%jsd1 = bd%jsd-1
   bd%ied1 = bd%ied+1
   bd%jed1 = bd%jed+1

   bd%isd2 = bd%isd-2
   bd%jsd2 = bd%jsd-2
   bd%ied2 = bd%ied+2
   bd%jed2 = bd%jed+2


end subroutine init_bounds
!-------------------------------------------------
! allocate atmos
!-------------------------------------------------
subroutine init_atmos(atm)
   type(fv_atmos_type), intent(INOUT) :: atm
   integer :: is, ie, isd, ied
   integer :: js, je, jsd, jed
   integer :: i
   character (len=60):: n, tc, hord, dp, iadv, mf
   is  = atm%bd%is
   js  = atm%bd%js
   isd = atm%bd%isd
   jsd = atm%bd%jsd
   ie  = atm%bd%ie
   je  = atm%bd%je
   ied = atm%bd%ied
   jed = atm%bd%jed
   atm%npx = atm%bd%npx
   atm%npy = atm%bd%npy

   write(n   ,'(i8)') atm%npx
   write(tc  ,'(i8)') atm%test_case
   write(hord,'(i8)') atm%hord
   write(dp  ,'(i8)') atm%dp
   write(iadv,'(i8)') atm%inner_adv
   write(mf  ,'(i8)') atm%mass_fixer
   atm%simulation_name = "tc"//trim(adjustl(tc))//"_N"//trim(adjustl(n))//"_hord"//&
   trim(adjustl(hord))//"_iadv"//trim(adjustl(iadv))//"_dp"//trim(adjustl(dp))//"_mf"//&
   trim(adjustl(mf))//"_"

   atm%grid_name = "grid_N"//trim(adjustl(n))


   allocate(atm%qa (isd:ied,jsd:jed, 1:nbfaces))
   allocate(atm%qa0(isd:ied,jsd:jed, 1:nbfaces))

   allocate(atm%uc (isd:ied+1,jsd:jed, 1:nbfaces))
   allocate(atm%uc0(isd:ied+1,jsd:jed, 1:nbfaces))

   allocate(atm%vc (isd:ied  ,jsd:jed+1, 1:nbfaces))
   allocate(atm%vc0(isd:ied,  jsd:jed+1, 1:nbfaces))

   allocate(atm%error_qa(is:ie,js:je, 1:nbfaces))
end subroutine init_atmos


!-------------------------------------------------
! init everything
!-------------------------------------------------
subroutine init_model(atm)
   type(fv_atmos_type), intent(inout) :: atm

   call init_bounds(atm%bd, atm%npx, atm%npy)
   call init_grid  (atm%gridstruct, atm%bd)
   call init_lagrange(atm%L, atm%bd)
   call init_atmos (atm)

end subroutine init_model


end module fv_control 
