module dyn_core
!========================================================================
! This module contains the routine that computes the PPM flux
!
! reference https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere/blob/main/model/dyn_core.F90
!========================================================================
use fv_arrays,  only: fv_grid_bounds_type, fv_grid_type, R_GRID, lagrange_poly, nbfaces
use test_cases, only: calc_winds
use tp_core,    only: fv_tp_2d
use sw_core,    only: time_averaged_cfl, compute_ra_x_and_ra_y, compute_xfx_and_yfx, div_mass_fixer
use fv_duogrid, only: ext_scalar_agrid

implicit none

contains
subroutine dy_core(qa, uc, uc_old, vc_old, vc, bd, gridstruct, time, time_centered, dt, dto2, test_case, hord, lim_fac, &
                   dp, inner_adv, mass_fixer, L)
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
 
   integer, intent(IN) :: test_case
   integer, intent(IN) :: hord
   integer, intent(IN) :: dp
   integer, intent(IN) :: inner_adv
   integer, intent(IN) :: mass_fixer

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

   integer :: is, ie
   integer :: js, je
   integer :: p

   is = bd%is
   ie = bd%ie
   js = bd%js
   je = bd%je
 
   ! winds
   call calc_winds(uc_old, vc_old, bd, gridstruct, time         , test_case)
   call calc_winds(uc    , vc    , bd, gridstruct, time_centered, test_case)

   ! compute time averaged cfl
   call time_averaged_cfl(gridstruct, bd, crx, cry, uc_old, vc_old, uc, vc, dp, dt)

   ! fill ghost cells
   call ext_scalar_agrid(qa, bd, L)

   ! compute coefficients need for fluxes
   call compute_xfx_and_yfx(xfx, yfx, crx, cry, gridstruct, bd, inner_adv)
   call compute_ra_x_and_ra_y(ra_x, ra_y, xfx, yfx, crx, cry, gridstruct, bd)

   ! compute the fluxes
   call fv_tp_2d(qa, crx, cry, hord, flux_x, flux_y, xfx, yfx, gridstruct, bd, ra_x, ra_y, lim_fac, inner_adv)

   ! fix the div mass
   call div_mass_fixer(bd, gridstruct, div, flux_x, flux_y, mass_fixer)

   ! update the solution
   !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
   !$OMP SHARED(qa, div, is, ie, js, je)
   qa(is:ie,js:je,:)  = qa(is:ie,js:je,:) - div(is:ie,js:je,:)
   !$OMP END PARALLEL WORKSHARE
end subroutine dy_core

end module dyn_core 
