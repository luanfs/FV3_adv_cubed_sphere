module test_cases
!========================================================================
! Module for ICS
!
! Reference
! https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere/blob/main/tools/test_cases.F90
!========================================================================
use fv_arrays,  only: fv_grid_bounds_type, fv_grid_type, fv_atmos_type, &
                      point_structure, R_GRID, pi, erad, eradi, day2sec, twopi, pio2, pio4, &
                      gravi, omega, sec2day, deg2rad, nbfaces
use fv_grid_tools, only: sph2cart
implicit none

contains

!-------------------------------------------------
! compute the ICS
!-------------------------------------------------
subroutine init_case(atm)
   type(fv_atmos_type), intent(INOUT) :: atm

   call init_scalar(atm%qa0, atm%bd, atm%gridstruct, atm%test_case)
   call init_winds (atm%uc0, atm%vc0, atm%bd, atm%gridstruct, atm%test_case)
   call calc_winds (atm%uc , atm%vc , atm%bd, atm%gridstruct, atm%time_centered, atm%test_case)
end subroutine init_case

!-------------------------------------------------
! compute initial scalar field
!-------------------------------------------------
subroutine init_scalar(qa, bd, gridstruct, test_case)
   type(fv_grid_type), target, intent(INOUT) :: gridstruct
   type(fv_grid_bounds_type), intent(IN) :: bd
   type(point_structure), pointer, dimension(:,:,:) :: agrid
   real(R_GRID), intent(OUT) :: qa(bd%isd:bd%ied,bd%jsd:bd%jed,1:nbfaces)
   real(R_GRID) :: lat, lon, latc, lonc
   real(R_GRID) :: u0, h0, alpha, e1(1:3), e2(1:3), pc(1:3), r
   integer, intent(IN) :: test_case
   integer :: is, ie
   integer :: js, je
   integer :: i, j, p
   integer :: isd, ied
   integer :: jsd, jed
 
   is = bd%is
   ie = bd%ie
   js = bd%js
   je = bd%je

   isd = bd%isd
   ied = bd%ied
   jsd = bd%jsd
   jed = bd%jed
 
   agrid => gridstruct%agrid

   if (test_case==1) then
      do p =1, nbfaces
         do i = is, ie
            do j = js, je
               lat = agrid(i,j,p)%lat
               lon = agrid(i,j,p)%lon
               alpha =  45.d0*deg2rad ! Rotation angle
               u0    =  erad*2.d0*pi/12.d0*sec2day ! Wind speed
               h0 = 2.94d0*10000.d0*gravi
               qa(i,j,p) = h0 - gravi*(erad*omega*u0 + u0*u0*0.5d0) &
               *(-dcos(lon)*dcos(lat)*dsin(alpha) + dsin(lat)*dcos(alpha))**2
            enddo
         enddo
      enddo

   else if (test_case==2) then
      ! center of Gaussian hill (a cube corner)
      e1(1) = 1.d0/dsqrt(3.d0)
      e1(2) = e1(1)
      e1(3) = e1(1)

      do p =1, nbfaces
         do i = is, ie
            do j = js, je
               e2 = agrid(i,j,p)%p
               r = norm2(e2-e1)
               qa(i,j,p) = 0.1d0 + 0.9d0*dexp(-10.d0*r*r)
            enddo
         enddo
      enddo

   else if (test_case==3 .or. test_case==4) then
      ! center of first Gaussian hill
      lonc = pio4
      latc = 0.d0
      call sph2cart(lonc, latc, pc(1), pc(2), pc(3))
      print*, pc
      do p =1, nbfaces
         do j=js,je
            do i=is,ie
               e2 = agrid(i,j,p)%p
               r = norm2(pc-e2)
               qa(i,j,p) = 0.1d0 + 0.9d0*(dexp(-5d0*r*r))
            enddo
         enddo
      enddo

      ! Second Gaussian hill
      lonc = -pio4
      latc = 0.d0
      call sph2cart(lonc, latc, pc(1), pc(2), pc(3))
      do p =1, nbfaces
         do j=js,je
            do i=is,ie
               e2 = agrid(i,j,p)%p
               r = norm2(pc-e2)
               qa(i,j,p) = qa(i,j,p) + 0.9d0*(dexp(-5d0*r*r))
            enddo
         enddo
      enddo

   else
      print*, 'Error in init_scalar: invalid test_case'
      stop

   endif


end subroutine init_scalar

!-------------------------------------------------
! compute initial winds
!-------------------------------------------------
subroutine init_winds(uc, vc, bd, gridstruct, test_case)
   type(fv_grid_type), target, intent(INOUT) :: gridstruct
   type(fv_grid_bounds_type), intent(IN) :: bd
   type(point_structure), pointer, dimension(:,:) :: cgrid
   real(R_GRID), intent(OUT) :: uc(bd%isd:bd%ied+1, bd%jsd:bd%jed  )
   real(R_GRID), intent(OUT) :: vc(bd%isd:bd%ied  , bd%jsd:bd%jed+1)
   integer, intent(IN) :: test_case

   call calc_winds(uc, vc, bd, gridstruct, 0.d0, test_case)

end subroutine init_winds


!-------------------------------------------------
! compute winds at a given time
!-------------------------------------------------
subroutine calc_winds(uc, vc, bd, gridstruct, time, test_case)
   type(fv_grid_type), target, intent(INOUT) :: gridstruct
   type(fv_grid_bounds_type), intent(IN) :: bd

   real(R_GRID), intent(INOUT) :: uc(bd%isd:bd%ied+1,bd%jsd:bd%jed  ,1:nbfaces)
   real(R_GRID), intent(INOUT) :: vc(bd%isd:bd%ied  ,bd%jsd:bd%jed+1,1:nbfaces)
   real(R_GRID), intent(IN) :: time
   integer, intent(IN) :: test_case

   type(point_structure), pointer, dimension(:,:,:) :: cgrid, dgrid
   real(R_GRID), pointer, dimension(:,:,:,:,:) :: c_l2contra, d_l2contra
   real(R_GRID):: ulon, vlat

   integer :: is, ie
   integer :: js, je
   integer :: i, j, p
   integer :: isd, ied
   integer :: jsd, jed
 
   is = bd%is
   ie = bd%ie
   js = bd%js
   je = bd%je
 
   isd = bd%isd
   ied = bd%ied
   jsd = bd%jsd
   jed = bd%jed
 
   cgrid => gridstruct%cgrid
   dgrid => gridstruct%dgrid
   c_l2contra => gridstruct%c_l2contra
   d_l2contra => gridstruct%d_l2contra

   ! u at cgrid
   do p = 1, nbfaces
      do i = isd, ied+1
         do j = jsd, jed
            call compute_wind(ulon, vlat, cgrid(i,j,p)%lon, cgrid(i,j,p)%lat, time, test_case)
            uc(i,j,p) = c_l2contra(1,1,i,j,p)*ulon + c_l2contra(1,2,i,j,p)*vlat
         enddo
      enddo
   enddo

   ! v at dgrid
   do p = 1, nbfaces
      do i = isd, ied
         do j = jsd, jed+1
            call compute_wind(ulon, vlat, dgrid(i,j,p)%lon, dgrid(i,j,p)%lat, time, test_case)
            vc(i,j,p) = d_l2contra(2,1,i,j,p)*ulon + d_l2contra(2,2,i,j,p)*vlat
         enddo
      enddo
   enddo

end subroutine calc_winds

!-------------------------------------------------
! compute wind at a given time and position
!-------------------------------------------------
subroutine compute_wind(u, v, lon, lat, t, test_case)
   real(R_GRID), intent(OUT) :: u, v
   real(R_GRID), intent(IN)  :: lon, lat, t
   integer, intent(IN) :: test_case
   real(R_GRID) :: Tf, u0, alpha, lonp

   Tf = 12.d0*day2sec
   select case (test_case)
      case(1,2)
         alpha =  45.d0*deg2rad ! Rotation angle
         u0    =  2.d0*pi*erad/Tf ! Wind speed
         u     =  u0*(dcos(lat)*dcos(alpha) + dsin(lat)*dcos(lon)*dsin(alpha))
         v     = -u0*dsin(lon)*dsin(alpha)

      case(3)
       u0   =  2.d0*pi*erad/Tf ! Wind speed
       lonp = lon-2.d0*pi*t/Tf
       u    = u0*(dsin((lonp))**2)*(dsin(2.*lat))*(dcos(pi*t/Tf))+u0*dcos(lat)
       v    = u0*(dsin(2*(lonp)))*(dcos(lat))*(dcos(pi*t/Tf))

      case(4)
       u0 =  2.d0*pi*erad/Tf ! Wind speed
       u  = -u0*(dsin((lon+pi)/2.d0)**2)*(dsin(2.d0*lat))*(dcos(lat)**2)*(dcos(pi*t/Tf))
       v  = (u0/2.d0)*(dsin((lon+pi)))*(dcos(lat)**3)*(dcos(pi*t/Tf))

      case default
         print*, 'error in compute_wind: invalid testcase, ', test_case
         stop
   end select
end subroutine compute_wind


end module test_cases
