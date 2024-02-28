module dyn_core
!========================================================================
! This module contains the routine that computes the PPM flux
!
! reference https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere/blob/main/model/dyn_core.F90
!========================================================================
use fv_arrays,  only: fv_grid_bounds_type, fv_grid_type, R_GRID, lagrange_poly, nbfaces
use test_cases, only: calc_winds
use tp_core,    only: fv_tp_2d
use sw_core,    only: time_averaged_cfl, compute_ra_x_and_ra_y, div_mass_fixer
use fv_duogrid, only: ext_scalar_agrid

implicit none

contains
subroutine dy_core(qa, uc, uc_old, vc_old, vc, bd, gridstruct, time, time_centered, dt, dto2, test_case, hord, lim_fac, &
                   dp, adv_scheme, L)
   type(fv_grid_bounds_type), intent(INOUT) :: bd
   type(fv_grid_type), target, intent(INOUT) :: gridstruct
   type(lagrange_poly), intent(INOUT) :: L
   real(R_GRID), intent(in) :: time, dt, dto2
   real(R_GRID), intent(in) :: time_centered
   real(R_GRID), intent(inout) :: lim_fac
   real(R_GRID), intent(inout) :: qa(bd%isd:bd%ied, bd%jsd:bd%jed, 1:nbfaces)
   real(R_GRID), intent(inout) :: uc    (bd%isd:bd%ied+1  , bd%jsd:bd%jed, 1:nbfaces)
   real(R_GRID), intent(inout) :: uc_old(bd%isd:bd%ied+1  , bd%jsd:bd%jed, 1:nbfaces)
   real(R_GRID), intent(inout) :: vc    (bd%isd:bd%ied, bd%jsd:bd%jed+1  , 1:nbfaces)
   real(R_GRID), intent(inout) :: vc_old(bd%isd:bd%ied, bd%jsd:bd%jed+1  , 1:nbfaces)
 
   real(R_GRID), pointer, dimension(:, :) :: dx_u, dy_u
   real(R_GRID), pointer, dimension(:, :) :: dx_v, dy_v

   real(R_GRID), pointer, dimension(:, :) :: sina_u, sina_v

   integer, intent(IN) :: test_case
   integer, intent(IN) :: hord
   integer, intent(IN) :: dp
   integer, intent(IN) :: adv_scheme
   integer :: p

   real(R_GRID) :: dx
   real(R_GRID) :: dy

   real(R_GRID) :: crx(bd%is:bd%ie+1, bd%jsd:bd%jed, 1:nbfaces)
   real(R_GRID) :: cry(bd%isd:bd%ied, bd%js:bd%je+1, 1:nbfaces)

   real(R_GRID) :: xfx(bd%is:bd%ie+1, bd%jsd:bd%jed, 1:nbfaces)
   real(R_GRID) :: yfx(bd%isd:bd%ied, bd%js:bd%je+1, 1:nbfaces)
 
   real(R_GRID) :: flux_x(bd%is:bd%ie+1, bd%js:bd%je  , 1:nbfaces)
   real(R_GRID) :: flux_y(bd%is:bd%ie  , bd%js:bd%je+1, 1:nbfaces)

   real(R_GRID) :: ra_x(bd%is:bd%ie  , bd%jsd:bd%jed, 1:nbfaces)
   real(R_GRID) :: ra_y(bd%isd:bd%ied, bd%js:bd%je  , 1:nbfaces )

   real(R_GRID) :: div(bd%is:bd%ie, bd%js:bd%je, 1:nbfaces)

   integer :: is, ie, isd, ied, ng
   integer :: js, je, jsd, jed
   integer :: i, j

   is  = bd%is
   ie  = bd%ie
   js  = bd%js
   je  = bd%je
 
   isd = bd%isd
   ied = bd%ied
   jsd = bd%jsd
   jed = bd%jed
 
   ng  = bd%ng

   dx_u  => gridstruct%dx_u
   dy_v  => gridstruct%dy_v

   dx_v  => gridstruct%dx_v
   dy_u  => gridstruct%dy_u

   sina_u  => gridstruct%sina_c
   sina_v  => gridstruct%sina_d

   dx = gridstruct%dx
   dy = gridstruct%dy

   ! winds
   call calc_winds(uc_old, vc_old, bd, gridstruct, time         , test_case)
   call calc_winds(uc    , vc    , bd, gridstruct, time_centered, test_case)

   ! fill ghost cells
   call ext_scalar_agrid(qa, bd, L)

   do p = 1, nbfaces
      ! compute time averaged cfl
      call time_averaged_cfl(gridstruct, bd, crx(:,:,p), cry(:,:,p), uc_old(:,:,p), vc_old(:,:,p), uc(:,:,p), vc(:,:,p), dp, dt)
    
      ! compute adv coeffs
      if(adv_scheme==1) then
         xfx(is:ie+1,jsd:jed,p) = crx(is:ie+1,jsd:jed,p)*dx_u(is:ie+1,jsd:jed)*dy_u(is:ie+1,jsd:jed)*sina_u(is:ie+1,jsd:jed)
         yfx(isd:ied,js:je+1,p) = cry(isd:ied,js:je+1,p)*dx_v(isd:ied,js:je+1)*dy_v(isd:ied,js:je+1)*sina_v(isd:ied,js:je+1)
      else
         xfx(is:ie+1,jsd:jed,p) = crx(is:ie+1,jsd:jed,p)*dx*dy
         yfx(isd:ied,js:je+1,p) = cry(isd:ied,js:je+1,p)*dx*dy
      endif 

      call compute_ra_x_and_ra_y(ra_x(:,:,p), ra_y(:,:,p), xfx(:,:,p), yfx(:,:,p), crx(:,:,p), cry(:,:,p), gridstruct, bd)

      call fv_tp_2d(qa(:,:,p), crx(:,:,p), cry(:,:,p), hord, flux_x(:,:,p), flux_y(:,:,p), &
                   xfx(:,:,p), yfx(:,:,p), gridstruct, bd, ra_x(:,:,p), ra_y(:,:,p), lim_fac, adv_scheme)

   enddo

   ! fix the mass
   call div_mass_fixer(bd, gridstruct, div, flux_x, flux_y)

   ! update the solution
   qa(is:ie,js:je,:)  = qa(is:ie,js:je,:) - div(is:ie,js:je,:)
end subroutine dy_core

end module dyn_core 
