module sw_core
!========================================================================
! This module contains the routine that computes the time averaged winds
!
! reference https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere/blob/main/model/sw_core.F90
!========================================================================
use fv_arrays, only: R_GRID, fv_grid_bounds_type, fv_grid_type, nbfaces
implicit none
contains


subroutine time_averaged_cfl(gridstruct, bd, crx, cry, uc_old, vc_old, uc, vc, dp, dt)
    !--------------------------------------------------
    ! Compute the time average CFL needed
    ! for the departure point scheme
    !
    !--------------------------------------------------
    type(fv_grid_bounds_type), intent(IN) :: bd
    type(fv_grid_type), intent(IN), target :: gridstruct

    real(R_GRID), intent(INOUT), dimension(bd%is:bd%ie+1  , bd%jsd:bd%jed  ) :: crx
    real(R_GRID), intent(INOUT), dimension(bd%isd:bd%ied  , bd%js:bd%je+1  ) :: cry

    real(R_GRID), intent(IN)   , dimension(bd%isd:bd%ied+1, bd%jsd:bd%jed  ) :: uc_old
    real(R_GRID), intent(IN)   , dimension(bd%isd:bd%ied  , bd%jsd:bd%jed+1) :: vc_old

    real(R_GRID), intent(IN)   , dimension(bd%isd:bd%ied+1, bd%jsd:bd%jed  ) :: uc
    real(R_GRID), intent(IN)   , dimension(bd%isd:bd%ied  , bd%jsd:bd%jed+1) :: vc

    real(R_GRID), intent(IN) :: dt

    integer, intent(IN):: dp ! inner departute point method !1 - Euler; 2-RK2

    ! aux
    integer :: i, j
    integer :: is,  ie,  js,  je

    is  = bd%is
    ie  = bd%ie
    js  = bd%js
    je  = bd%je

    call departure_cfl(gridstruct, bd, crx, cry, uc_old, vc_old, uc, vc, dp, dt)

end subroutine time_averaged_cfl

subroutine departure_cfl(gridstruct, bd, crx, cry, &
                             uc_old, vc_old, uc, vc, dp, dt)
    !--------------------------------------------------
    ! Compute the departure CFL for a given 
    ! departure point scheme (dp)
    !--------------------------------------------------
    type(fv_grid_bounds_type), intent(IN) :: bd
    type(fv_grid_type), intent(IN), target :: gridstruct 

    real(R_GRID), intent(INOUT), dimension(bd%is:bd%ie+1, bd%jsd:bd%jed  ) :: crx
    real(R_GRID), intent(INOUT), dimension(bd%isd:bd%ied, bd%js:bd%je+1  ) :: cry

    real(R_GRID), intent(IN)   , dimension(bd%isd:bd%ied+1, bd%jsd:bd%jed  ) :: uc_old
    real(R_GRID), intent(IN)   , dimension(bd%isd:bd%ied  , bd%jsd:bd%jed+1) :: vc_old

    real(R_GRID), intent(IN)   , dimension(bd%isd:bd%ied+1, bd%jsd:bd%jed  ) :: uc
    real(R_GRID), intent(IN)   , dimension(bd%isd:bd%ied  , bd%jsd:bd%jed+1) :: vc

    real(R_GRID), dimension(bd%is-1:bd%ie+2  , bd%jsd:bd%jed  ) :: crx_time_centered
    real(R_GRID), dimension(bd%isd:bd%ied  , bd%js-1:bd%je+2  ) :: cry_time_centered
    real(R_GRID), dimension(bd%is:bd%ie+1  , bd%jsd:bd%jed  ) :: crx_old
    real(R_GRID), dimension(bd%isd:bd%ied  , bd%js:bd%je+1  ) :: cry_old


    real(R_GRID), intent(IN) :: dt

    integer, intent(IN):: dp ! departute point method !1 - Euler; 2-RK2

    ! aux
    real(R_GRID) :: a, a1, a2, c1, c2
    integer :: i, j
    integer :: is,  ie,  js,  je
    integer :: isd, ied, jsd, jed

    is  = bd%is
    ie  = bd%ie
    js  = bd%js
    je  = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed

    select case (dp)
      case (1)
         ! CFL for RK1
         call  compute_cfl(gridstruct, bd, crx, cry, uc, vc, dt, 0)

      case (2)
         ! CFL for RK2
         call  compute_cfl(gridstruct, bd, crx_old          , cry_old          , uc_old, vc_old, dt, 0)
         call  compute_cfl(gridstruct, bd, crx_time_centered, cry_time_centered, uc    , vc    , dt, 1)

         ! RK2
         ! cfl for dp in x direction
         do j = jsd, jed
             do i = is, ie+1
                ! Linear interpolation weight
                a = crx_old(i,j)*0.5d0
                ! Upwind linear interpolation
                if (a>0.d0) then
                   c1 = crx_time_centered(i-1,j)
                   c2 = crx_time_centered(i,j)
                   a1 = a
                   a2 = 1.d0-a
                else
                   c1 = crx_time_centered(i,j)
                   c2 = crx_time_centered(i+1,j)
                   a1 = 1.d0+a
                   a2 = -a
                end if
                crx(i,j) = a1*c1 + a2*c2
             end do
          end do

          ! cfl for dp in y direction
          do i = isd, ied
             do j = js, je+1
                ! Linear interpolation weight
                a = cry_old(i,j)*0.5d0
                ! Upwind linear interpolation
                if (a>0.d0) then
                   c1 = cry_time_centered(i,j-1)
                   c2 = cry_time_centered(i,j)
                   a1 = a
                   a2 = 1.d0-a
                else
                   c1 = cry_time_centered(i,j)
                   c2 = cry_time_centered(i,j+1)
                   a1 = 1.d0+a
                   a2 = -a
                end if
                cry(i,j) = a1*c1 + a2*c2
             end do
          end do

    end select

end subroutine departure_cfl

subroutine compute_cfl(gridstruct, bd, crx, cry, uc, vc, dt, h)
   !--------------------------------------------------
   ! Compute CFL in x and y directions at C grid
   !--------------------------------------------------
   type(fv_grid_bounds_type), intent(IN) :: bd
   type(fv_grid_type), intent(IN), target :: gridstruct 
   real(R_GRID), intent(INOUT), dimension(bd%is-h:bd%ie+1+h  , bd%jsd:bd%jed  ) :: crx
   real(R_GRID), intent(INOUT), dimension(bd%isd:bd%ied  , bd%js-h:bd%je+1+h  ) :: cry
   real(R_GRID), intent(IN)   , dimension(bd%isd:bd%ied+1, bd%jsd:bd%jed  ) :: uc
   real(R_GRID), intent(IN)   , dimension(bd%isd:bd%ied  , bd%jsd:bd%jed+1) :: vc
   real(R_GRID), intent(IN) :: dt
   integer, intent(IN):: h
   real(R_GRID), pointer, dimension(:, :) :: dxc, dyc

   ! aux
   integer :: i, j
   integer :: is,  ie,  js,  je
   integer :: isd, ied, jsd, jed
   real(R_GRID) :: dx, dy

   is  = bd%is
   ie  = bd%ie
   js  = bd%js
   je  = bd%je
   isd = bd%isd
   ied = bd%ied
   jsd = bd%jsd
   jed = bd%jed

   dxc => gridstruct%dx_u
   dyc => gridstruct%dy_v
 
   ! Compute CFL at timestep n
   do j=jsd,jed
      do i=is-h,ie+h+1
        crx(i,j) = uc(i,j)*dt/dxc(i,j)
     enddo
   enddo

   do j=js-h,je+h+1
      do i=isd,ied
        cry(i,j) = vc(i,j)*dt/dyc(i,j)
     enddo
   enddo
end subroutine compute_cfl

subroutine compute_ra_x_and_ra_y(ra_x, ra_y, xfx, yfx, crx, cry, gridstruct, bd)
    type(fv_grid_bounds_type), intent(IN) :: bd
    type(fv_grid_type), intent(IN), target :: gridstruct
    real(R_GRID), intent(INOUT), dimension(bd%is:bd%ie  , bd%jsd:bd%jed  ) :: ra_x
    real(R_GRID), intent(IN)   , dimension(bd%is:bd%ie+1, bd%jsd:bd%jed  ) :: xfx
    real(R_GRID), intent(IN)   , dimension(bd%is:bd%ie+1, bd%jsd:bd%jed  ) :: crx

    real(R_GRID), intent(INOUT), dimension(bd%isd:bd%ied, bd%js:bd%je    ) :: ra_y
    real(R_GRID), intent(IN)   , dimension(bd%isd:bd%ied, bd%js:bd%je+1  ) :: yfx
    real(R_GRID), intent(IN)   , dimension(bd%isd:bd%ied, bd%js:bd%je+1  ) :: cry

    ! Local:
    integer :: i, j
    integer :: is, ie, js, je
    integer :: isd, ied, jsd, jed
 
    real(R_GRID), pointer, dimension(:,:)   :: area

    area => gridstruct%area

    is  = bd%is
    ie  = bd%ie
    js  = bd%js
    je  = bd%je
    isd = bd%isd
    ied = bd%ied
    jsd = bd%jsd
    jed = bd%jed


    do j=jsd,jed
       do i=is,ie
          ra_x(i,j) = area(i,j) + xfx(i,j) - xfx(i+1,j)
       enddo
    enddo
    do j=js,je
       do i=isd,ied
          ra_y(i,j) = area(i,j) + yfx(i,j) - yfx(i,j+1)
       enddo
    enddo
end subroutine compute_ra_x_and_ra_y

subroutine div_mass_fixer(bd, gridstruct, div, flux_x, flux_y)
   type(fv_grid_bounds_type), intent(INOUT) :: bd
   type(fv_grid_type), intent(INOUT), target :: gridstruct
   real(R_GRID), intent(INOUT) :: flux_x(bd%is:bd%ie+1, bd%js:bd%je  , 1:nbfaces)
   real(R_GRID), intent(INOUT):: flux_y(bd%is:bd%ie  , bd%js:bd%je+1, 1:nbfaces)
   real(R_GRID), intent(INOUT):: div(bd%is:bd%ie, bd%js:bd%je, 1:nbfaces)
   integer :: p
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

   call  average_flux_at_cube_intefaces(bd, flux_x, flux_y)
   do p = 1, nbfaces
      div(is:ie,js:je,p) =  (flux_x(is+1:ie+1,js:je,p)-flux_x(is:ie,js:je,p))+ &
                            (flux_y(is:ie,js+1:je+1,p)-flux_y(is:ie,js:je,p))

      div(is:ie,js:je,p) = div(is:ie,js:je,p)*gridstruct%rarea(is:ie,js:je)
   enddo

end subroutine div_mass_fixer

subroutine average_flux_at_cube_intefaces(bd, flux_x, flux_y)
   !---------------------------------------------------------------------------------
   ! AVERAGE_FLUX_AT_CUBE_INTERFACES
   !
   ! Mass fixer that average the flux values at cube interfaces
   !--------------------------------------------------------------------------------
   type(fv_grid_bounds_type), intent(INOUT) :: bd
   real(R_GRID), intent(INOUT) :: flux_x(bd%is:bd%ie+1, bd%js:bd%je  , 1:nbfaces)
   real(R_GRID), intent(INOUT) :: flux_y(bd%is:bd%ie  , bd%js:bd%je+1, 1:nbfaces)
   real(R_GRID) :: a, b
   integer :: is, ie, isd, ied
   integer :: js, je, jsd, jed
 
   is  = bd%is
   ie  = bd%ie
   js  = bd%js
   je  = bd%je
 
   isd = bd%isd
   ied = bd%ied
   jsd = bd%jsd
   jed = bd%jed
 
   a = 0.5d0
   b = 0.5d0
   !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
   !$OMP SHARED(flux_x, flux_y, is, ie, js, je, isd, ied, jsd, jed, a, b)
   ! Average panels 1-2,2-3,3-4,4-1
   flux_x(ie+1,js:je,1:3) = a*flux_x(ie+1,js:je,1:3) + b*flux_x(is,js:je,2:4)
   flux_x(is,js:je,2:4)     =flux_x(ie+1,js:je,1:3)

   flux_x(ie+1,js:je,4) = a*flux_x(ie+1,js:je,4) + b*flux_x(is,js:je,1)
   flux_x(is,js:je,1)  =flux_x(ie+1,js:je,4)

   ! Average panels 1-5
   flux_y(is:ie,js,5)   = a*flux_y(is:ie,js,5) + b*flux_y(is:ie,je+1,1)
   flux_y(is:ie,je+1,1) = flux_y(is:ie,js,5)

   ! Average panels 2-5
   flux_x(ie+1,js:je,5) = a*flux_x(ie+1,js:je,5) - b*flux_y(is:ie,je+1,2)
   flux_y(is:ie,je+1,2) = -flux_x(ie+1,js:je,5)

   ! Average panels 3-5
   flux_y(is:ie,je+1,5) = a*flux_y(is:ie,je+1,5) - b*flux_y(ie:is:-1,je+1,3)
   flux_y(is:ie,je+1,3) = -flux_y(ie:is:-1,je+1,5)

   ! Average panels 4-5
   flux_x(is,js:je,5)    = a*flux_x(is,js:je,5) + b*(flux_y(ie:is:-1,je+1,4))
   flux_y(is:ie,je+1,4) = (flux_x(is,je:js:-1,5))

   ! Average panels 1-6
   flux_y(is:ie,je+1,6) = a*flux_y(is:ie,je+1,6) + b*flux_y(is:ie,js,1)
   flux_y(is:ie,js,1)   = flux_y(is:ie,je+1,6)

   ! Average panels 2-6
   flux_y(is:ie,js,2)     = a*flux_y(is:ie,js,2) + b*flux_x(ie+1,je:js:-1,6)
   flux_x(ie+1,js:je,6) = flux_y(ie:is:-1,js,2)

   ! Average panels 3-6
   flux_y(is:ie,js,3) = a*flux_y(is:ie,js,3) - b*flux_y(ie:is:-1,js,6)
   flux_y(is:ie,js,6) = -flux_y(ie:is:-1,js,3)

   ! Average panels 4-6
   flux_y(is:ie,js,4) = -a*flux_x(is,js:je,6) + b*flux_y(is:ie,js,4)
   flux_x(is,js:je,6) = -flux_y(is:ie,js,4)

   !$OMP END PARALLEL WORKSHARE
end subroutine average_flux_at_cube_intefaces

end module sw_core 
