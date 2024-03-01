module fv_grid_tools
!========================================================================
! Module for data allocation
!
! Reference
! https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere/blob/main/tools/fv_control.F90
!========================================================================
use fv_arrays,  only: fv_grid_bounds_type, fv_grid_type, fv_atmos_type, &
                      point_structure, R_GRID, erad, pi, pio4, pio2, nbfaces
implicit none

contains

!-------------------------------------------------
! allocate grid - bd must be filled
!-------------------------------------------------
subroutine init_grid(gridstruct, bd)
   type(fv_grid_type), target, intent(INOUT) :: gridstruct
   type(fv_grid_bounds_type), intent(IN)    :: bd
   type(point_structure), pointer, dimension(:,:,:) :: agrid
   type(point_structure), pointer, dimension(:,:,:) :: bgrid
   type(point_structure), pointer, dimension(:,:,:) :: cgrid
   type(point_structure), pointer, dimension(:,:,:) :: dgrid

   real(R_GRID), pointer, dimension(:,:) :: area
   real(R_GRID), pointer, dimension(:,:) :: rarea
   real(R_GRID), pointer, dimension(:,:) :: mt
   real(R_GRID), pointer, dimension(:,:) :: dx_u, dy_u
   real(R_GRID), pointer, dimension(:,:) :: dx_v, dy_v

   real(R_GRID):: L, aref, Rref, dx, x, y
   real(R_GRID), allocatable :: gridline_a(:)
   real(R_GRID), allocatable :: gridline_b(:)

   real(R_GRID), allocatable :: angle_a(:)
   real(R_GRID), allocatable :: angle_b(:)
   real(R_GRID) :: p1(2), p2(2), p3(2), p4(2)


   integer :: is, ie, isd, ied
   integer :: js, je, jsd, jed
   integer :: i, j, p

   is  = bd%is
   js  = bd%js
   isd = bd%isd
   jsd = bd%jsd
   ie  = bd%ie
   je  = bd%je
   ied = bd%ied
   jed = bd%jed

   gridstruct%npx = bd%npx
   gridstruct%npy = bd%npy

   allocate(gridline_b(isd:ied+1))
   allocate(gridline_a(isd:ied))

   allocate(angle_b(isd:ied+1))
   allocate(angle_a(isd:ied))

   allocate(gridstruct%agrid(isd:ied  , jsd:jed  , 1:nbfaces))
   allocate(gridstruct%bgrid(isd:ied+1, jsd:jed+1, 1:nbfaces))
   allocate(gridstruct%cgrid(isd:ied+1, jsd:jed  , 1:nbfaces))
   allocate(gridstruct%dgrid(isd:ied  , jsd:jed+1, 1:nbfaces))
   allocate(gridstruct% area(isd:ied  , jsd:jed  ))
   allocate(gridstruct%rarea(isd:ied  , jsd:jed  ))
   allocate(gridstruct%   mt(isd:ied  , jsd:jed  ))

   allocate(gridstruct%sina_c(isd:ied+1, jsd:jed  ))
   allocate(gridstruct%cosa_c(isd:ied+1, jsd:jed  ))
   allocate(gridstruct%sina_d(isd:ied  , jsd:jed+1))
   allocate(gridstruct%cosa_d(isd:ied  , jsd:jed+1))
 
   allocate(gridstruct%dx_u(isd:ied+1, jsd:jed  ))
   allocate(gridstruct%dy_v(isd:ied  , jsd:jed+1))
 
   allocate(gridstruct%dx_v(isd:ied  , jsd:jed+1))
   allocate(gridstruct%dy_u(isd:ied+1, jsd:jed  ))

   allocate(gridstruct%c_contra2l(1:2,1:2,isd:ied+1, jsd:jed, 1:nbfaces  ))
   allocate(gridstruct%c_l2contra(1:2,1:2,isd:ied+1, jsd:jed, 1:nbfaces  ))

   allocate(gridstruct%d_contra2l(1:2,1:2,isd:ied, jsd:jed+1, 1:nbfaces))
   allocate(gridstruct%d_l2contra(1:2,1:2,isd:ied, jsd:jed+1, 1:nbfaces))
 
   agrid => gridstruct%agrid
   bgrid => gridstruct%bgrid
   cgrid => gridstruct%cgrid
   dgrid => gridstruct%dgrid

   area  => gridstruct%area
   rarea => gridstruct%rarea
   mt    => gridstruct%mt
   dx_u   => gridstruct%dx_u
   dx_v   => gridstruct%dx_v
   dy_u   => gridstruct%dy_u
   dy_v   => gridstruct%dy_v

   aref = pio4
   Rref = 1.d0
   dx = 2.d0*Rref*aref/bd%npx
   gridstruct%dx = dx
   gridstruct%dy = dx

   ! compute b grid local coordinates
   do i = isd, ied+1
      gridline_b(i) = -aref + (i-1.d0)*dx
      angle_b(i) = dtan(gridline_b(i))
   enddo

   do i = isd, ied
      gridline_a(i) = -aref + (i-0.5d0)*dx 
      angle_a(i) = dtan(gridline_a(i))
   enddo

   !--------------------------------------------------------------------------------------------
   ! compute bgrid
   do p = 1, nbfaces
      do i = isd, ied+1
         x = angle_b(i)
         do j = jsd, jed+1
            y = angle_b(j)
            call equidistant_gnomonic_map(bgrid(i,j,p)%p, x, y, p)
            call cart2sph ( bgrid(i,j,p)%p(1), bgrid(i,j,p)%p(2), bgrid(i,j,p)%p(3), bgrid(i,j,p)%lon, bgrid(i,j,p)%lat)
         enddo
      enddo
   enddo


   !--------------------------------------------------------------------------------------------

   !--------------------------------------------------------------------------------------------
   ! compute agrid
   do p = 1, nbfaces
      do i = isd, ied
         do j = jsd, jed
            agrid(i,j,p)%p = (bgrid(i,j,p)%p + bgrid(i,j+1,p)%p + bgrid(i+1,j,p)%p + bgrid(i+1,j+1,p)%p)*0.25d0 
            agrid(i,j,p)%p = agrid(i,j,p)%p/norm2(agrid(i,j,p)%p)
            call cart2sph ( agrid(i,j,p)%p(1), agrid(i,j,p)%p(2), agrid(i,j,p)%p(3), agrid(i,j,p)%lon, agrid(i,j,p)%lat)
         enddo
      enddo
   enddo

   !--------------------------------------------------------------------------------------------
   ! compute agrid
   do p = 1, nbfaces
      do i = isd, ied
         x = angle_a(i)
         do j = jsd, jed
            y = angle_a(j)
            !call equidistant_gnomonic_map(agrid(i,j,p)%p, x, y, p)
            !call cart2sph ( agrid(i,j,p)%p(1), agrid(i,j,p)%p(2), agrid(i,j,p)%p(3), agrid(i,j,p)%lon, agrid(i,j,p)%lat)
         enddo
      enddo
   enddo


   !--------------------------------------------------------------------------------------------

   !--------------------------------------------------------------------------------------------
   ! compute cgrid
   do p = 1, nbfaces
      do i = isd, ied+1
         do j = jsd, jed
            cgrid(i,j,p)%p = (bgrid(i,j+1,p)%p + bgrid(i,j,p)%p)*0.5d0 
            cgrid(i,j,p)%p = cgrid(i,j,p)%p/norm2(cgrid(i,j,p)%p)
            call cart2sph ( cgrid(i,j,p)%p(1), cgrid(i,j,p)%p(2), cgrid(i,j,p)%p(3), cgrid(i,j,p)%lon, cgrid(i,j,p)%lat)
         enddo
      enddo
   enddo
   !--------------------------------------------------------------------------------------------

   !--------------------------------------------------------------------------------------------
   ! compute dgrid
   do p = 1, nbfaces
      do i = isd, ied
         do j = jsd, jed+1
            dgrid(i,j,p)%p = (bgrid(i+1,j,p)%p + bgrid(i,j,p)%p)*0.5d0 
            dgrid(i,j,p)%p = dgrid(i,j,p)%p/norm2(dgrid(i,j,p)%p)
            call cart2sph ( dgrid(i,j,p)%p(1), dgrid(i,j,p)%p(2), dgrid(i,j,p)%p(3), dgrid(i,j,p)%lon, dgrid(i,j,p)%lat)
         enddo
      enddo
   enddo
   !--------------------------------------------------------------------------------------------

   !--------------------------------------------------------------------------------------------
   ! compute distances
   do i = is-1, ie+2
      do j = jsd, jed
         dx_u(i,j) = arclen(agrid(i,j,1)%p, agrid(i-1,j,1)%p, erad )
      enddo
   enddo

   do j = js-1, jed-1
      do i = isd, ied
         dy_v(i,j) = arclen(agrid(i,j,1)%p, agrid(i,j-1,1)%p, erad )
      enddo
   enddo

   do i = isd, ied+1
      do j = jsd, jed
         dy_u(i,j) = arclen(bgrid(i,j,1)%p, bgrid(i,j+1,1)%p, erad )
      enddo
   enddo

   do j = jsd, jed+1
      do i = isd, ied
         dx_v(i,j) = arclen(bgrid(i,j,1)%p, bgrid(i+1,j,1)%p, erad )
      enddo
   enddo
 
   !--------------------------------------------------------------------------------------------
   ! compute agrid metric term
   do i = isd, ied
      x = gridline_a(i)
      do j = jsd, jed
         y = gridline_a(j)
         call metricterm(2, x, y, mt(i,j), erad)
         ! area(i,j) = mt(i,j)*dx*dx
         !rarea(i,j) = 1.d0/area(i,j)
      enddo
   enddo
   !--------------------------------------------------------------------------------------------
 

   !--------------------------------------------------------------------------------------------
   ! compute agrid area
   do i = isd, ied
      do j = jsd, jed
         p1(1:2) = [bgrid(i  ,j  ,1)%lon, bgrid(i  ,j  ,1)%lat]
         p2(1:2) = [bgrid(i+1,j  ,1)%lon, bgrid(i+1,j  ,1)%lat]
         p3(1:2) = [bgrid(i+1,j+1,1)%lon, bgrid(i+1,j+1,1)%lat]
         p4(1:2) = [bgrid(i  ,j+1,1)%lon, bgrid(i  ,j+1,1)%lat]
         area(i,j) = get_area(p1, p4, p2, p3, erad)
         rarea(i,j) = 1.d0/area(i,j)
      enddo
   enddo

   !--------------------------------------------------------------------------------------------

   call compute_conversion_matrices(bd, gridstruct)

   deallocate(gridline_a)
   deallocate(gridline_b)
   deallocate(angle_a)
   deallocate(angle_b)
end subroutine init_grid

 subroutine equidistant_gnomonic_map(p, x, y, panel)
 !---------------------------------------------------------------
 ! this routine computes the gnomonic mapping based on the equidistant projection
 ! defined by rancic et al (96) for each panel
 ! - x, y are the variables defined in [-a, a].
 ! - the projection is applied on the points (x,y)
 ! - returns the cartesian coordinates of the
 ! projected point p.
 !
 ! references: 
 ! - rancic, m., purser, r.j. and mesinger, f. (1996), a global shallow-water model using an expanded
 !  spherical cube: gnomonic versus conformal coordinates. q.j.r. meteorol. soc., 122: 959-982. 
 !  https://doi.org/10.1002/qj.49712253209
 ! - nair, r. d., thomas, s. j., & loft, r. d. (2005). a discontinuous galerkin transport scheme on the
 ! cubed sphere, monthly weather review, 133(4), 814-828. retrieved feb 7, 2022, 
 ! from https://journals.ametsoc.org/view/journals/mwre/133/4/mwr2890.1.xml
 !
 !---------------------------------------------------------------
 real(kind=R_GRID), intent(in) :: x ! local coordinates
 real(kind=R_GRID), intent(in) :: y ! local coordinates
 real(kind=R_GRID), intent(out) :: p(1:3) ! projected point

 ! panel
 integer, intent(in) :: panel

 ! compute the cartesian coordinates for each panel
 ! with the aid of the auxiliary variables  
 select case(panel)
   case(1)
     p(1) =  1.d0
     p(2) =  x
     p(3) =  y

   case(2)
     p(1) = -x
     p(2) =  1.d0
     p(3) =  y

   case(3)
     p(1) = -1.d0
     p(2) = -x
     p(3) =  y

   case(4)
     p(1) =  x
     p(2) = -1.d0
     p(3) =  y      

   case(5)
     p(1) = -y
     p(2) =  x
     p(3) =  1.d0

   case(6)
     p(1) =  y
     p(2) =  x
     p(3) = -1.d0

   case default
     print*, 'error on equidistant_gnomonic_map: invalid panel'
     stop
 end select
 p = p/norm2(p)
 return
 end subroutine equidistant_gnomonic_map

 subroutine inverse_equidistant_gnomonic_map(p, x, y, panel)
 !---------------------------------------------------------------
 ! Given a panel, this routine computes the inverse of the equidistant gnomonic map
 !---------------------------------------------------------------
 real(kind=8), intent(out) :: x ! cube coordinates
 real(kind=8), intent(out) :: y ! cube coordinates
 real(kind=8), intent(in) :: p(1:3) ! point on the sphere

 ! panel
 integer, intent(in) :: panel

 ! aux vars
 real(kind=8) :: a

 ! compute the local coordinates for each panel
 select case(panel)
   case(1)
     x = p(2)/p(1)
     y = p(3)/p(1)
   case(2)
     x = -p(1)/p(2)
     y =  p(3)/p(2)
   case(3)
     x =  p(2)/p(1)
     y = -p(3)/p(1)
   case(4)
     x = -p(1)/p(2)
     y = -p(3)/p(2)
   case(5)
     x =  p(2)/p(3)
     y = -p(1)/p(3)
   case(6)
     x = -p(2)/p(3)
     y = -p(1)/p(3)
   case default
     print*, 'error on inverse_equidistant_gnomonic_map: invalid panel'
     stop
 end select
 return
 end subroutine inverse_equidistant_gnomonic_map

 subroutine metricterm(grid_type, x, y, mt, erad)
    integer, intent(in):: grid_type
    ! type of grid: must be 0, 1 or 2.
    real(kind=R_GRID), intent(in) :: x, y, erad
    real(kind=R_GRID), intent(out) :: mt
    real(kind=R_GRID) :: aref, Rref
    real(kind=R_GRID) :: tanx, tany, r
    real(kind=R_GRID) :: cosx, cosy

    if(grid_type==0) then
       ! equiedge grid parameters
       Rref = dsqrt(2.d0)
       aref = dasin(1.d0/dsqrt(3.d0))
    else if (grid_type==2) then
       ! equiangular grid parameters
       Rref = 1.d0
       aref = pi*0.25d0
    endif

    tanx = dtan(x)
    tany = dtan(y)
    cosx = dcos(x)
    cosy = dcos(y)
    
    r = dsqrt(1 + Rref*Rref*tanx*tanx + Rref*Rref*tany*tany)
    mt = Rref*Rref*erad*erad/r**3
    mt = mt/(cosx*cosx*cosy*cosy)

 end subroutine metricterm

  subroutine compute_conversion_matrices(bd, gridstruct)
     !---------------------------------------------------------------------
     !
     ! COMPUTE_CONVERSION_MATRICES
     !
     ! This routine computes the matrices that performs
     ! the conversions between contravariant, covariant and latlon representation
     ! of tangent vectors on the sphere
     !---------------------------------------------------------------------
      type(fv_grid_type), intent(IN), target :: gridstruct
      type(fv_grid_bounds_type), intent(IN) :: bd
      ! Integer auxs
      integer :: i, j, p ! 2D grid counters
      integer :: is , js , ie , je
      integer :: isd, jsd, ied, jed

      ! Real aux vars
      real(R_GRID) :: elon(1:3), elat(1:3)
      real(R_GRID) :: ex(1:3), ey(1:3) ! vectors
      real(R_GRID) :: p0(1:3), px(1:3), py(1:3) ! vectors
      real(R_GRID) :: lat, lon
      real(R_GRID) :: a11, a12, a21, a22, det
      real(R_GRID), pointer, dimension(:,:,:,:,:) :: d_contra2l, d_l2contra
      real(R_GRID), pointer, dimension(:,:,:,:,:) :: c_contra2l, c_l2contra
      type(point_structure), pointer, dimension(:,:,:) :: agrid, bgrid, cgrid, dgrid
      real(R_GRID), pointer, dimension(:,:) :: cosa_c, cosa_d
      real(R_GRID), pointer, dimension(:,:) :: sina_c, sina_d

      c_contra2l => gridstruct%c_contra2l
      d_contra2l => gridstruct%d_contra2l

      c_l2contra => gridstruct%c_l2contra
      d_l2contra => gridstruct%d_l2contra

      agrid => gridstruct%agrid
      bgrid => gridstruct%bgrid
      cgrid => gridstruct%cgrid
      dgrid => gridstruct%dgrid

      cosa_c => gridstruct%cosa_c
      cosa_d => gridstruct%cosa_d

      sina_c => gridstruct%sina_c
      sina_d => gridstruct%sina_d


      is = bd%is
      ie = bd%ie
      js = bd%js
      je = bd%je
      isd = bd%isd
      ied = bd%ied
      jsd = bd%jsd
      jed = bd%jed

      do p = 1, nbfaces
         ! C grid
         do i = isd, ied+1
            do j = jsd, jed
               ! get latlon tangent vector
               call unit_vect_latlon_ext([cgrid(i,j,p)%lon, cgrid(i,j,p)%lat], elon, elat)
               if(i<ied+1)then
                  p0 = cgrid(i,j,p)%p
                  px = agrid(i,j,p)%p
                  py = bgrid(i,j+1,p)%p

                  ! unit vector - x direction
                  ex = px-p0
                  ex = proj_vec_sphere(ex, p0)
                  ex = ex/norm2(ex)

                  ! unit vector - y direction
                  ey = py-p0
                  ey = proj_vec_sphere(ey, p0)
                  ey = ey/norm2(ey)
               else
                  p0 = cgrid(i,j,p)%p
                  px = agrid(i-1,j,p)%p
                  py = bgrid(i,j+1,p)%p

                  ! unit vector - x direction
                  ex = -px+p0
                  ex = proj_vec_sphere(ex, p0)
                  ex = ex/norm2(ex)
    
                  ! unit vector - y direction
                  ey = py-p0
                  ey = proj_vec_sphere(ey, p0)
                  ey = ey/norm2(ey)
               endif

               if(p==1) then
                  cosa_c(i,j) = dot_product(ex, ey)
                  sina_c(i,j) = dsqrt(1.d0 - cosa_c(i,j)**2) 
               endif

               a11 = dot_product(ex, elon)  
               a12 = dot_product(ey, elon)  
               a21 = dot_product(ex, elat)  
               a22 = dot_product(ey, elat) 
               det = a11*a22 - a21*a12

               ! Contra to latlon matrix
               c_contra2l(1,1,i,j,p) = a11
               c_contra2l(1,2,i,j,p) = a12
               c_contra2l(2,1,i,j,p) = a21
               c_contra2l(2,2,i,j,p) = a22

               ! latlon to contra matrix
               c_l2contra(1,1,i,j,p) =  a22/det
               c_l2contra(1,2,i,j,p) = -a12/det
               c_l2contra(2,1,i,j,p) = -a21/det
               c_l2contra(2,2,i,j,p) =  a11/det
            enddo
         enddo

         !print*, minval(cosa_u), maxval(cosa_u)


         ! D grid
         do i = isd, ied
            do j = jsd, jed+1
               ! get latlon tangent vector
               call unit_vect_latlon_ext([dgrid(i,j,p)%lon, dgrid(i,j,p)%lat], elon, elat)
               if(j<jed+1)then
                  p0 = dgrid(i,j,p)%p
                  py = agrid(i,j,p)%p
                  px = bgrid(i+1,j,p)%p

                  ! unit vector - x direction
                  ex = px-p0
                  ex = proj_vec_sphere(ex, p0)
                  ex = ex/norm2(ex)

                  ! unit vector - y direction
                  ey = py-p0
                  ey = proj_vec_sphere(ey, p0)
                  ey = ey/norm2(ey)

               else
                  p0 = dgrid(i,j,p)%p
                  py = agrid(i,j-1,p)%p
                  px = bgrid(i+1,j,p)%p

                  ! unit vector - x direction
                  ex = px-p0
                  ex = proj_vec_sphere(ex, p0)
                  ex = ex/norm2(ex)
    
                  ! unit vector - y direction
                  ey = -py+p0
                  ey = proj_vec_sphere(ey, p0)
                  ey = ey/norm2(ey)
               endif

               if(p==1) then
                  cosa_d(i,j) = dot_product(ex, ey)
                  sina_d(i,j) = dsqrt(1.d0 - cosa_d(i,j)**2) 
               endif

               a11 = dot_product(ex, elon)  
               a12 = dot_product(ey, elon)  
               a21 = dot_product(ex, elat)  
               a22 = dot_product(ey, elat) 
               det = a11*a22 - a21*a12

               ! Contra to latlon matrix
               d_contra2l(1,1,i,j,p) = a11
               d_contra2l(1,2,i,j,p) = a12
               d_contra2l(2,1,i,j,p) = a21
               d_contra2l(2,2,i,j,p) = a22

               ! latlon to contra matrix
               d_l2contra(1,1,i,j,p)=  a22/det
               d_l2contra(1,2,i,j,p)= -a12/det
               d_l2contra(2,1,i,j,p)= -a21/det
               d_l2contra(2,2,i,j,p)=  a11/det
            enddo
         enddo
      enddo
  end subroutine compute_conversion_matrices

  function proj_vec_sphere(v, p)
    !-----------------------------------------------------------
    !  Projects a vector 'v' on the plane tangent to a sphere
    !   Uses the the tangent plane relative to the unit sphere's
    !   point 'p', in cartesian coords
    !-----------------------------------------------------------
    real (kind=R_GRID), intent(in), dimension(1:3) :: v
    real (kind=R_GRID), intent(in), dimension(1:3) :: p
    real (kind=R_GRID), dimension(1:3)  :: proj_vec_sphere

    proj_vec_sphere(1:3)=&
         v(1:3)-dot_product(v,p)*p(1:3)/norm2(p)

    return
  end function proj_vec_sphere

!===============================================================================================
!   the following routines were taken from imodel  https://github.com/pedrospeixoto/imodel
!===============================================================================================

subroutine sph2cart (lon, lat, x, y, z )
  !------------------------------------------------------------------------------------
  ! sph2cart 
  !
  !     transforms geographical coordinates (lat,lon) to cartesian coordinates.
  !     similar to stripack's trans
  !
  !    input: lat, latitudes of the node in radians [-pi/2,pi/2]
  !           lon, longitudes of the nodes in radians [-pi,pi]
  !
  !    output:  x, y, z, the coordinates in the range -1 to 1. 
  !                    x**2 + y**2 + z**2 = 1 
  !---------------------------------------------------------------------
  real (R_GRID), intent(in) :: lon
  real (R_GRID), intent(in) :: lat
  real (R_GRID), intent(out) :: x
  real (R_GRID), intent(out) :: y
  real (R_GRID), intent(out) :: z
  real (R_GRID):: coslat

  coslat = dcos (lat)
  x = coslat * dcos (lon)
  y = coslat * dsin (lon)
  z = dsin (lat)

  return
end subroutine sph2cart

subroutine cart2sph ( x, y, z, lon, lat )
  !---------------------------------------------------------------------
  ! cart2sph 
  !     transforms  cartesian coordinates to geographical (lat,lon) coordinates .
  !     similar to stripack's scoord
  !
  !    input:  x, y, z, the coordinates in the range -1 to 1. 
  !
  !    output, lon, longitude of the node in radians [-pi,pi].
  !                       lon=0 if point lies on the z-axis.  
  !    output, lat, latitude of the node in radians [-pi/2,pi/2].
  !                       lat=0 if   x**2 + y**2 + z**2 = 0.
  !------------------------------------------------------------------------------------
  real    (R_GRID), intent(in) :: x
  real    (R_GRID), intent(in) :: y
  real    (R_GRID), intent(in) :: z
  real    (R_GRID), intent(out) :: lat
  real    (R_GRID), intent(out) :: lon
  real    (R_GRID):: pnrm

  pnrm = dsqrt ( x **2 + y**2 + z**2 )
  if ( x /= 0.0d+00 .or. y /= 0.0d+00 ) then
     lon = datan2 ( y, x )
  else
     lon = 0.0d+00
  end if

  if ( pnrm /= 0.0d+00 ) then
     lat = dasin ( z / pnrm )
  else
     print*, "cart2sph warning: point not in the unit sphere. norm= ", pnrm
     lat = 0.0d+00
  end if

  return
end subroutine cart2sph


  function arclen(p, q, radius)
    !-----------------------------------------------------------
    ! ARCLEN from ssrfpack by r. renka            
    !
    !   This function computes the arc-length (angle in radians)            
    !    between a pair of points on the unit sphere. It is similar to
    !    arcdistxyz, but uses a calculation to avoid inverse cosine function
    !    p,q = arrays of length 3 containing the x, y, and z             
    !             coordinates (in that order) of points on the              
    !             unit sphere.                                                 
    !   Returns the angle in radians between the unit vectors              
    !       p and q.  0 <= arclen <= pi.                       
    !---------------------------------------------------------
    real (R_GRID), intent(in) :: p(3)
    real (R_GRID), intent(in) :: q(3)
    real(kind=R_GRID), intent(in), optional :: radius
    real (R_GRID):: arclen

    !Dot product
    real (R_GRID):: d

    d = dot_product(p+q, p+q)

    if (d==0.) then 
       ! p and q are separated by 180 degrees.
       arclen = pi
    elseif (d>=4.) then 
       ! p and q coincide.
       arclen = 0.d0
    else 
       arclen = 2.d0 * datan (dsqrt ( (4.d0 - d) / d) )
    endif

    if ( present(radius) ) then
      arclen = radius * arclen
    endif

   return 
  end function arclen

!-----------------------------------------------------------------------------
! routines from  GFDL cubed sphere code
!-----------------------------------------------------------------------------
  subroutine unit_vect_latlon_ext(pp, elon, elat)
    real(R_GRID), intent(IN)  :: pp(2)
    real(R_GRID), intent(OUT) :: elon(3), elat(3)
    real(R_GRID):: lon, lat
    real(R_GRID):: sin_lon, cos_lon, sin_lat, cos_lat

    lon = pp(1)
    lat = pp(2)

    sin_lon = dsin(lon)
    cos_lon = dcos(lon)
    sin_lat = dsin(lat)
    cos_lat = dcos(lat)

    elon(1) = -sin_lon
    elon(2) = cos_lon
    elon(3) = 0.d0

    elat(1) = -sin_lat*cos_lon
    elat(2) = -sin_lat*sin_lon
    elat(3) = cos_lat

  end subroutine unit_vect_latlon_ext


subroutine latlon2xyz(p, e, id)
!
! Routine to map (lon, lat) to (x,y,z)
!
 real(kind=R_GRID), intent(in) :: p(2)
 real(kind=R_GRID), intent(out):: e(3)
 integer, optional, intent(in):: id   ! id=0 do nothing; id=1, right_hand

 integer n
 real (R_GRID):: q(2)
 real (R_GRID):: e1, e2, e3

    do n=1,2
       q(n) = p(n)
    enddo

    e1 = cos(q(2)) * cos(q(1))
    e2 = cos(q(2)) * sin(q(1))
    e3 = sin(q(2))
!-----------------------------------
! Truncate to the desired precision:
!-----------------------------------
    e(1) = e1
    e(2) = e2
    e(3) = e3

 end subroutine latlon2xyz

 real(kind=R_GRID) function get_area(p1, p4, p2, p3, radius)
!-----------------------------------------------
 real(kind=R_GRID), intent(in), dimension(2):: p1, p2, p3, p4
 real(kind=R_GRID), intent(in), optional:: radius
!-----------------------------------------------
 real(kind=R_GRID) e1(3), e2(3), e3(3)
 real(kind=R_GRID) ang1, ang2, ang3, ang4

! S-W: 1
       call latlon2xyz(p1, e1)   ! p1
       call latlon2xyz(p2, e2)   ! p2
       call latlon2xyz(p4, e3)   ! p4
       ang1 = spherical_angle(e1, e2, e3)
!----
! S-E: 2
!----
       call latlon2xyz(p2, e1)
       call latlon2xyz(p3, e2)
       call latlon2xyz(p1, e3)
       ang2 = spherical_angle(e1, e2, e3)
!----
! N-E: 3
!----
       call latlon2xyz(p3, e1)
       call latlon2xyz(p4, e2)
       call latlon2xyz(p2, e3)
       ang3 = spherical_angle(e1, e2, e3)
!----
! N-W: 4
!----
       call latlon2xyz(p4, e1)
       call latlon2xyz(p3, e2)
       call latlon2xyz(p1, e3)
       ang4 = spherical_angle(e1, e2, e3)

       if ( present(radius) ) then
            get_area = (ang1 + ang2 + ang3 + ang4 - 2.*pi) * radius**2
       else
            get_area = ang1 + ang2 + ang3 + ang4 - 2.*pi
       endif

 end function get_area


 real(kind=R_GRID) function spherical_angle(p1, p2, p3)

!           p3
!         /
!        /
!       p1 ---> angle
!         \
!          \
!           p2

 real(kind=R_GRID) p1(3), p2(3), p3(3)

 real (R_GRID):: e1(3), e2(3), e3(3)
 real (R_GRID):: px, py, pz
 real (R_GRID):: qx, qy, qz
 real (R_GRID):: angle, ddd
 integer n

  do n=1,3
     e1(n) = p1(n)
     e2(n) = p2(n)
     e3(n) = p3(n)
  enddo

!-------------------------------------------------------------------
! Page 41, Silverman's book on Vector Algebra; spherical trigonmetry
!-------------------------------------------------------------------
! Vector P:
   px = e1(2)*e2(3) - e1(3)*e2(2)
   py = e1(3)*e2(1) - e1(1)*e2(3)
   pz = e1(1)*e2(2) - e1(2)*e2(1)
! Vector Q:
   qx = e1(2)*e3(3) - e1(3)*e3(2)
   qy = e1(3)*e3(1) - e1(1)*e3(3)
   qz = e1(1)*e3(2) - e1(2)*e3(1)

   ddd = (px*px+py*py+pz*pz)*(qx*qx+qy*qy+qz*qz)

   if ( ddd <= 0.0d0 ) then
        angle = 0.d0
   else
        ddd = (px*qx+py*qy+pz*qz) / sqrt(ddd)
        if ( abs(ddd)>1.d0) then
             angle = 2.d0*atan(1.0)    ! 0.5*pi
           !FIX (lmh) to correctly handle co-linear points (angle near pi or 0)
           if (ddd < 0.d0) then
              angle = 4.d0*atan(1.0d0) !should be pi
           else
              angle = 0.d0
           end if
        else
             angle = acos( ddd )
        endif
   endif

   spherical_angle = angle
 end function spherical_angle
end module fv_grid_tools
