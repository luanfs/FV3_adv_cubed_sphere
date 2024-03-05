module fv_arrays
!========================================================================
! This module contains the arrays data structure
!
! Reference:
! https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere/blob/main/model/fv_arrays.F90
!========================================================================
public
integer, public, parameter :: R_GRID = selected_real_kind(12,100)
real(R_GRID), public, parameter :: pio4  = datan(1.d0)
real(R_GRID), public, parameter :: pi    = 4.d0*pio4
real(R_GRID), public, parameter :: pio2  = pi*0.5d0
real(R_GRID), public, parameter :: twopi = pi*2.d0
real(R_GRID), public, parameter :: erad  = 6.37122e6
real(R_GRID), public, parameter :: eradi = 1.d0/erad
real(R_GRID), public, parameter :: day2sec = 86400.d0
real(R_GRID), public, parameter :: sec2day = 1.d0/day2sec
character(len=32) :: datadir = "data/"
character(len=32) :: griddir = "grid/"
character(len=32) :: pardir  = "par/"

!Gravitational accerlaration of the Earth (m/s^2)
real(R_GRID), parameter  :: grav    = 9.8066499999999994d0
real(R_GRID), parameter  :: gravity = grav
real(R_GRID), parameter  :: gravi   = 1.d0/grav
real(R_GRID), parameter  :: gravo2   = grav*0.5d0

! Angular velocity of the Earth (rot/s)
real (kind=8), parameter :: omega   = 7.2921d0*10.d0**(-5)
real (kind=8), parameter :: rotatn   = 7.2921d0*10.d0**(-5)


!Degrees to radians coversion (multiply to obtain conversion)
real(kind=8), parameter :: deg2rad = pi / 180.d0

!Radians to Degrees coversion (multiply to obtain conversion)
real(kind=8), parameter :: rad2deg = 1.d0/deg2rad

integer :: nbfaces = 6

!-----------------------------------------------------------------------------
! domain bounds
!-----------------------------------------------------------------------------
type fv_grid_bounds_type
   integer :: is,  ie   ! x direction interior cell indexes (for agrid; interior cgrid ranges from is to ie+1)
   integer :: js,  je   ! y direction interior cell indexes (for agrid; interior dgrid ranges from js to je+1)
   integer :: isd, ied  ! x direction indexes including ghost cells (for agrid; cgrid ranges from isd to ied+1)
   integer :: jsd, jed  ! y direction indexes including ghost cells (for agrid; dgrid ranges from jsd to jed+1)

   integer :: isd1, ied1
   integer :: jsd1, jed1
   integer :: isd2, ied2
   integer :: jsd2, jed2
 
   integer :: ng  = 3   ! ghost cell layers
   integer :: npx       ! number of interior cells in x direction
   integer :: npy       ! number of interior cells in y direction
end type fv_grid_bounds_type


!---------------------------------------------------------
! Type for Lagrange polynomials at ghost cells on the cubed-sphere
!---------------------------------------------------------
type lagrange_poly
    real (kind=R_GRID), allocatable  :: y_Agrid_given(:)    ! A grid points using panel 2 coordinates
    real (kind=R_GRID), allocatable  :: angles_a(:)         ! A grid point angles
    real (kind=R_GRID), allocatable  :: tan_angles_a(:)     ! A grid point tan(angles)
    real (kind=R_GRID), allocatable  :: y_Agrid_target(:,:) ! Nodes where we perform interpolation
    real (kind=R_GRID), allocatable  :: p_Agrid_ydir(:,:,:) ! Lagrange polynomials at nodes

    ! interpolation weights for corner diagonal  
    real (kind=R_GRID), allocatable  :: p_corner1(:,:)
    real (kind=R_GRID), allocatable  :: p_corner2(:,:)

    ! buffers
    real (kind=R_GRID), allocatable  :: E_buffer(:,:,:)
    real (kind=R_GRID), allocatable  :: N_buffer(:,:,:)
    real (kind=R_GRID), allocatable  :: S_buffer(:,:,:)
    real (kind=R_GRID), allocatable  :: W_buffer(:,:,:)

    real (kind=R_GRID), allocatable  :: E_buffer_local(:,:,:)
    real (kind=R_GRID), allocatable  :: N_buffer_local(:,:,:)
    real (kind=R_GRID), allocatable  :: S_buffer_local(:,:,:)
    real (kind=R_GRID), allocatable  :: W_buffer_local(:,:,:)

    ! aux vars
    real (kind=R_GRID), allocatable  :: var(:,:,:)

    real (kind=R_GRID), allocatable  :: dy
    real (kind=R_GRID), allocatable  :: aref
    real (kind=R_GRID), allocatable  :: Rref

    ! stencil
    integer, allocatable :: stencil_start_a(:,:), stencil_end_a(:,:) ! A grid
    integer :: stencil_offset

    ! polynomial degree
    integer :: degree

    ! order (=degree+1)
    integer :: order

    ! N
    integer :: N
end type lagrange_poly

!-------------------------------------------------
! point structure
!-------------------------------------------------
type point_structure
   ! Cartesian coordinates (R^3)
   real(kind=R_GRID), dimension(1:3) :: p

   ! latlon 
   real(kind=R_GRID) :: lat, lon
end type point_structure

!-------------------------------------------------
! grid structure
!-------------------------------------------------
type fv_grid_type
   type(point_structure), allocatable, dimension(:,:,:) :: agrid
   type(point_structure), allocatable, dimension(:,:,:) :: bgrid
   type(point_structure), allocatable, dimension(:,:,:) :: cgrid
   type(point_structure), allocatable, dimension(:,:,:) :: dgrid

   real(R_GRID), allocatable, dimension(:,:,:,:,:) :: c_contra2l, d_contra2l
   real(R_GRID), allocatable, dimension(:,:,:,:,:) :: c_l2contra, d_l2contra

   real(R_GRID), allocatable ::  mt(:,:) ! metric term
   real(R_GRID), allocatable :: rmt(:,:) ! metric term
   real(R_GRID), allocatable ::  area(:,:)
   real(R_GRID), allocatable :: rarea(:,:)
   real(R_GRID), allocatable :: sina_c(:,:), sina_d(:,:)
   real(R_GRID), allocatable :: cosa_c(:,:), cosa_d(:,:)
   real(R_GRID), allocatable :: dx_u(:,:), dy_u(:,:)
   real(R_GRID), allocatable :: dx_v(:,:), dy_v(:,:)
   real(R_GRID), allocatable :: dxa(:,:) , dya(:,:)

   real(R_GRID), allocatable :: rdx_u(:,:), rdy_u(:,:)
   real(R_GRID), allocatable :: rdx_v(:,:), rdy_v(:,:)
   real(R_GRID), allocatable :: rdxa(:,:) , rdya(:,:)


   real(R_GRID) :: dx, dy
   integer :: npx       ! number of interior cells (x direction)
   integer :: npy       ! number of interior cells (y direction)
   integer :: grid_type  ! 0-equiedge grid; 2-equiangular
end type fv_grid_type

!-------------------------------------------------
! atmosphere structure
!-------------------------------------------------
type fv_atmos_type
   type(fv_grid_type) :: gridstruct
   type(fv_grid_bounds_type) :: bd
   type(lagrange_poly) :: L
   integer  :: test_case
   integer  :: hord
   integer  :: dp
   integer  :: inner_adv
   integer  :: mass_fixer
   integer  :: nplot = 0
   integer  :: nplots
   integer  :: plotstep
   integer  :: panel
   logical  :: first_step=.true.
   real(R_GRID) :: lim_fac = 1.d0

   ! time vars
   real(R_GRID) :: dt            ! time step
   real(R_GRID) :: dto2          ! = (dt/2)
   real(R_GRID) :: time          ! current time
   real(R_GRID) :: time_centered ! current time + dto2
   real(R_GRID) :: Tf            ! final time
   integer  :: total_tsteps           ! total of time steps

   ! Fields
   real(R_GRID), allocatable :: qa (:,:,:) ! A grid scalar field
   real(R_GRID), allocatable :: qa0(:,:,:) ! A grid scalar field (IC)

   real(R_GRID), allocatable :: uc (:,:,:) ! (uc) are mostly used as the C grid winds
   real(R_GRID), allocatable :: uc0(:,:,:) ! (uc) are mostly used as the C grid winds (IC)
   real(R_GRID), allocatable :: uc_old(:,:,:) ! 

   real(R_GRID), allocatable :: vc (:,:,:) ! (vc) are mostly used as the C grid winds
   real(R_GRID), allocatable :: vc0(:,:,:) ! (vc) are mostly used as the C grid winds (IC)
   real(R_GRID), allocatable :: vc_old(:,:,:) ! 


   ! erros vars
   real(R_GRID), allocatable :: error_qa(:,:,:)
   real(R_GRID) :: linf_error_qa, l1_error_qa, l2_error_qa

   ! cfl
   real(R_GRID) :: cfl, cfl_x, cfl_y

   ! mass vars
   real(R_GRID) :: mass_qa0, mass_qa, mass_qa_var

   integer :: npx       ! number of interior cells
   integer :: npy       ! number of interior cells
   character(len=128) :: simulation_name
   character(len=128) :: grid_name
end type fv_atmos_type

 
end module fv_arrays
