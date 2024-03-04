module fv_duogrid
!========================================================================
! Module to handle with ghost cell interpolations
!
!========================================================================
use fv_arrays,  only: fv_grid_bounds_type,  &
                      R_GRID, nbfaces, lagrange_poly, pio4, pi, grav
use fv_grid_tools, only: equidistant_gnomonic_map, inverse_equidistant_gnomonic_map
implicit none

contains

!--------------------------------------------------------------------------
subroutine fill_haloes(q, bd, N_buffer, S_buffer, W_buffer, E_buffer)
   type(fv_grid_bounds_type), intent(IN) :: bd
   real(R_GRID), intent(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed, 1:nbfaces)
   real(R_GRID), intent(INOUT) :: E_buffer(bd%ie+1:bd%ied2, bd%jsd2:bd%jed2, 1:nbfaces)
   real(R_GRID), intent(INOUT) :: W_buffer(bd%isd2:bd%is-1, bd%jsd2:bd%jed2, 1:nbfaces)
   real(R_GRID), intent(INOUT) :: N_buffer(bd%isd2:bd%ied2, bd%je+1:bd%jed2, 1:nbfaces)
   real(R_GRID), intent(INOUT) :: S_buffer(bd%isd2:bd%ied2, bd%jsd2:bd%js-1, 1:nbfaces)
   integer :: is, ie, isd, ied, ng
   integer :: js, je, jsd, jed
   integer :: p, east, north, south, west, j

   is = bd%is
   js = bd%js
   ie = bd%ie
   je = bd%je
   isd = bd%isd2
   jsd = bd%jsd2
   ied = bd%ied2
   jed = bd%jed2
   ng = bd%ng+2

   ! --------------------- Panel 1 ----------------------------
   p = 1
   north = 5
   south = 6
   east  = 2
   west  = 4
 
   ! Data of panel 1 from east
   E_buffer(ie+1:ied,js:je,p) = q(is:is+ng-1, is:ie, east) ! Panel 2

   ! Data of panel 1 from  west
   W_buffer(isd:is-1,js:je,p) = q(ie-ng+1:ie, js:je, west) ! Panel 4

   ! Data of panel 1 from north
   N_buffer(is:ie,je+1:jed,p) = q(is:ie, js:js+ng-1, north) ! Panel 5

   ! Data of panel 1 from south
   S_buffer(is:ie,jsd:js-1,p) = q(is:ie, je-ng+1:je, south) ! Panel 6

   ! --------------------- Panel 2 ----------------------------
   p = 2
   north = 5
   south = 6
   east  = 3
   west  = 1

   ! Data of panel 2 from east
   E_buffer(ie+1:ied,js:je,p) = q(is:is+ng-1, js:je, east) ! Panel 3

   ! Data of panel 2 from west
   W_buffer(isd:is-1,js:je,p) = q(ie-ng+1:ie, js:je, west) ! Panel 1

   ! Data of panel 2 from north
   N_buffer(is:ie,je+1:jed,p) = transpose(Q(je:je-ng+1:-1, js:je, north)) ! Panel 5

   ! Data of panel 2 from south
   S_buffer(is:ie,jsd:js-1,p) = transpose(Q(je-ng+1:je, je:js:-1,south)) ! Panel 6

   ! --------------------- Panel 3 ----------------------------
   p = 3
   north = 5
   south = 6
   east  = 4
   west  = 2

   ! Data of panel 3 from east
   E_buffer(ie+1:ied,js:je,p) = q(is:is+ng-1, js:je, east) ! Panel 4

   ! Data of panel 3 from west
   W_buffer(isd:is-1,js:je,p) = q(ie-ng+1:ie, js:je, west) ! Panel 2

   ! Data of panel 3 from north
   N_buffer(is:ie,je+1:jed,p) = q(ie:is:-1, je:je-ng+1:-1, north) ! Panel 5

   ! Data of panel 3 from south
   S_buffer(is:ie,jsd:js-1,p) = q(ie:is:-1, js+ng-1:js:-1, south) ! Panel 6

   ! --------------------- Panel 4 ----------------------------
   p = 4
   north = 5
   south = 6
   east  = 1
   west  = 3

   ! Data of panel 4 from east
   E_buffer(ie+1:ied,js:je,p) = q(is:is+ng-1, js:je, east) ! Panel 1

   ! Data of panel 4 from west
   W_buffer(isd:is-1,js:je,p) = q(ie-ng+1:ie, js:je, west) ! Panel 3

   ! Data of panel 4 from north
   N_buffer(is:ie,je+1:jed,p) = transpose(q(is:is+ng-1, je:js:-1, north)) ! Panel 5

   ! Data of panel 4 from south
   S_buffer(is:ie,jsd:js-1,p)  = transpose(q(is+ng-1:is:-1, js:je, south)) ! Panel 6

   ! --------------------- Panel 5 ----------------------------
   p = 5
   north = 3
   south = 1
   east  = 2
   west  = 4

   ! Data of panel 5 from east
   E_buffer(ie+1:ied,js:je,p) = transpose(Q(is:ie, je:je-ng+1:-1, east)) ! Panel 2

   ! Data of panel 5 from west
   W_buffer(isd:is-1,js:je,p) = transpose(Q(ie:is:-1,je-ng+1:je, west)) ! Panel 4

   ! Data of panel 5 from north
   N_buffer(is:ie,je+1:jed,p) = Q(ie:is:-1, je:je-ng+1:-1, north) ! Panel 3

   ! Data of panel 5 from south
   S_buffer(is:ie,jsd:js-1,p) = Q(is:ie, je-ng+1:je, south) ! Panel 1


   ! --------------------- Panel 6 ----------------------------
   p = 6
   north = 1
   south = 3
   east  = 2
   west  = 4

   ! Data of panel 6 from east
   E_buffer(ie+1:ied,js:je,p) = transpose(Q(ie:is:-1, js:js+ng-1, east)) ! Panel 2

   ! Data of panel 6 from west
   W_buffer(isd:is-1,js:je,p) = transpose(Q(is:ie,js+ng-1:js:-1, west)) ! Panel 4

   ! Data of panel 6 from north
   N_buffer(is:ie,je+1:jed,p) = Q(is:ie, js:js+ng-1, north) ! Panel 3

   ! Data of panel 6 from south
   S_buffer(is:ie,jsd:js-1,p) = Q(ie:is:-1, js+ng-1:js:-1, south) ! Panel 1

end subroutine fill_haloes

!--------------------------------------------------------------------------
subroutine apply_buffer(q, bd, N_buffer, S_buffer, W_buffer, E_buffer, ng)
   type(fv_grid_bounds_type), intent(IN) :: bd
   real(R_GRID), intent(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed, 1:nbfaces)
   real(R_GRID), intent(INOUT) :: E_buffer(bd%ie+1:bd%ied2, bd%jsd2:bd%jed2, 1:nbfaces)
   real(R_GRID), intent(INOUT) :: W_buffer(bd%isd2:bd%is-1, bd%jsd2:bd%jed2, 1:nbfaces)
   real(R_GRID), intent(INOUT) :: N_buffer(bd%isd2:bd%ied2, bd%je+1:bd%jed2, 1:nbfaces)
   real(R_GRID), intent(INOUT) :: S_buffer(bd%isd2:bd%ied2, bd%jsd2:bd%js-1, 1:nbfaces)
   integer :: is, ie, isd, ied, ng
   integer :: js, je, jsd, jed
   integer :: p, east, north, south, west

   is = bd%is
   js = bd%js
   ie = bd%ie
   je = bd%je
   isd = bd%isd
   jsd = bd%jsd
   ied = bd%ied
   jed = bd%jed
   
   q(ie+1:ied,js:je,:) = E_buffer(ie+1:ied,js:je,:)
   q(isd:is-1,js:je,:) = W_buffer(isd:is-1,js:je,:)
   q(isd:ied,je+1:jed,:) = N_buffer(isd:ied,je+1:jed,:)
   q(isd:ied,jsd:js-1,:) = S_buffer(isd:ied,jsd:js-1,:)


end subroutine apply_buffer
!
!
!--------------------------------------------------------------------------
! extend agrid field
subroutine ext_scalar_agrid(q, bd, L)
   type(fv_grid_bounds_type), intent(INOUT) :: bd
   type(lagrange_poly), intent(inout) :: L
   real(R_GRID), intent(INOUT) :: q(bd%isd:bd%ied, bd%jsd:bd%jed, 1:nbfaces)
   integer :: is, ie, isd, ied, ng
   integer :: js, je, jsd, jed
   integer :: i,j

   is = bd%is
   js = bd%js
   ie = bd%ie
   je = bd%je
   isd = bd%isd
   jsd = bd%jsd
   ied = bd%ied
   jed = bd%jed
   ng = bd%ng

   ! fill haloes
   call fill_haloes(q, bd, L%N_buffer, L%S_buffer, L%W_buffer, L%E_buffer)

   ! Fill off corner points
   call cube_rmp_off_corner(bd, L, L%N_buffer, L%S_buffer, L%E_buffer, L%W_buffer, ng)

   ! Fill off diagonal corner points
   call cube_rmp_corner_offdiag(bd, L, L%N_buffer, L%S_buffer, L%E_buffer, L%W_buffer, ng)

   ! Fill diagonal corner points
   call cube_rmp_corner_diag(bd, L, L%N_buffer, L%S_buffer, L%E_buffer, L%W_buffer, ng)

   ! Apply bufffer
   call apply_buffer(q, bd, L%N_buffer, L%S_buffer, L%W_buffer, L%E_buffer, ng)

end subroutine ext_scalar_agrid
!--------------------------------------------------------------------------

subroutine cube_rmp_off_corner(bd, L, N_buffer, S_buffer, E_buffer, W_buffer, ng)
  type(fv_grid_bounds_type), intent(IN) :: bd
  type(lagrange_poly), target, intent(inout) :: L
  integer, intent(in) :: ng
  real(R_GRID), intent(INOUT) :: E_buffer(bd%ie+1:bd%ied2, bd%jsd2:bd%jed2, 1:nbfaces)
  real(R_GRID), intent(INOUT) :: W_buffer(bd%isd2:bd%is-1, bd%jsd2:bd%jed2, 1:nbfaces)
  real(R_GRID), intent(INOUT) :: N_buffer(bd%isd2:bd%ied2, bd%je+1:bd%jed2, 1:nbfaces)
  real(R_GRID), intent(INOUT) :: S_buffer(bd%isd2:bd%ied2, bd%jsd2:bd%js-1, 1:nbfaces)

  real(kind=R_GRID), pointer, dimension(:, :, :) :: poly
  integer, pointer, dimension(:, :) :: stencil_start

  real(R_GRID), pointer, dimension(:, :, :) :: N_buffer_local
  real(R_GRID), pointer, dimension(:, :, :) :: S_buffer_local
  real(R_GRID), pointer, dimension(:, :, :) :: E_buffer_local
  real(R_GRID), pointer, dimension(:, :, :) :: W_buffer_local

  logical :: rmp_s, rmp_e, rmp_n, rmp_w
  integer :: is, ie, js, je, isd, ied, jsd, jed
  integer :: i, j,  n, g, g2, d, h, k, p

  !--- assign parameters
  is = bd%is
  ie = bd%ie
  js = bd%js
  je = bd%je

  poly => L%p_Agrid_ydir
  stencil_start => L%stencil_start_a

  N_buffer_local => L%N_buffer_local
  S_buffer_local => L%S_buffer_local
  E_buffer_local => L%E_buffer_local
  W_buffer_local => L%W_buffer_local

  N_buffer_local = N_buffer
  S_buffer_local = S_buffer
  E_buffer_local = E_buffer
  W_buffer_local = W_buffer

  !$OMP PARALLEL DO &
  !$OMP DEFAULT(NONE) & 
  !$OMP SHARED(poly, L, stencil_start) &
  !$OMP SHARED(S_buffer, N_buffer, W_buffer, E_buffer) &
  !$OMP SHARED(S_buffer_local, N_buffer_local, W_buffer_local, E_buffer_local) &
  !$OMP SHARED(is, ie, js, je, nbfaces, ng) &
  !$OMP PRIVATE(i, j, g, g2, h) &
  !$OMP SCHEDULE(static) 
  do p = 1, nbfaces
     !--- south
     do g = 1, ng
        g2 = ng-g+1
        j  = js-g2
        do i = is, ie
           S_buffer(i, j, p) = 0.d0
           do n = 1, L%order
              S_buffer(i, j, p) = S_buffer(i, j, p) + S_buffer_local(stencil_start(i,g2)+n-1, j, p)*poly(i,g2,n)
           enddo
        enddo
     enddo

     !--- north
     do g = 1, ng
        j = je+g
        do i = is, ie
           N_buffer(i, j, p) = 0.d0
           do n = 1, L%order
              N_buffer(i, j, p) = N_buffer(i, j, p) + N_buffer_local(stencil_start(i,g)+n-1, j, p)*poly(i,g,n)
           enddo
        enddo
     enddo

     !--- west
     do g = 1, ng
        g2 = ng-g+1
        i  = is-g2
        do j = js, je
           W_buffer(i, j, p) = 0.d0
           do n = 1, L%order
              W_buffer(i, j, p) = W_buffer(i, j, p) + W_buffer_local(i, stencil_start(j,g2)+n-1,p)*poly(j,g2,n)
           enddo
        enddo
     enddo

     !--- east
     do g = 1, ng
        i = ie+g
        h = g-1
        do j = js, je
           E_buffer(i, j, p) = 0.d0
           do n = 1, L%order
              E_buffer(i, j, p) = E_buffer(i, j, p) + E_buffer_local(i, stencil_start(j,g)+n-1, p)*poly(j,g,n)
           enddo
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

end subroutine cube_rmp_off_corner

subroutine cube_rmp_corner_offdiag(bd, L, N_buffer, S_buffer, E_buffer, W_buffer, ng)
!--------------------------------------------------------------------------------------------------
!
! This routine does interpolation at ghost cells at corners but off diagonal
!--------------------------------------------------------------------------------------------------
  type(fv_grid_bounds_type), intent(IN) :: bd
  type(lagrange_poly), target, intent(inout) :: L
  integer, intent(in) :: ng 
  real(R_GRID), intent(INOUT) :: E_buffer(bd%ie+1:bd%ied2, bd%jsd2:bd%jed2, 1:nbfaces)
  real(R_GRID), intent(INOUT) :: W_buffer(bd%isd2:bd%is-1, bd%jsd2:bd%jed2, 1:nbfaces)
  real(R_GRID), intent(INOUT) :: N_buffer(bd%isd2:bd%ied2, bd%je+1:bd%jed2, 1:nbfaces)
  real(R_GRID), intent(INOUT) :: S_buffer(bd%isd2:bd%ied2, bd%jsd2:bd%js-1, 1:nbfaces)

  real(kind=R_GRID), pointer, dimension(:, :, :) :: poly
  integer, pointer, dimension(:, :) :: stencil_start

  real(R_GRID), pointer, dimension(:, :, :) :: N_buffer_local
  real(R_GRID), pointer, dimension(:, :, :) :: S_buffer_local
  real(R_GRID), pointer, dimension(:, :, :) :: E_buffer_local
  real(R_GRID), pointer, dimension(:, :, :) :: W_buffer_local


  integer :: is, ie, js, je, isd, ied, jsd, jed
  integer :: i, j, p, jj, ii, h, n, g, g2

  !--- assign parameters
  is = bd%is
  ie = bd%ie
  js = bd%js
  je = bd%je
  isd = bd%isd
  ied = bd%ied
  jsd = bd%jsd
  jed = bd%jed

  poly => L%p_Agrid_ydir
  stencil_start => L%stencil_start_a

  N_buffer_local => L%N_buffer_local
  S_buffer_local => L%S_buffer_local
  E_buffer_local => L%E_buffer_local
  W_buffer_local => L%W_buffer_local


  !---------------------------------------------------------------------------------------------------------------------------
  !  SE corner (diag = 1,5,9)
  !        1--2--3 |
  !        |     | |
  !        4  5  6 |
  !        |     | |
  !        7--8--9 |
  !  ---------------------------------------------------
  !
  !---------------------------------------------------------------------------------------------------------------------------
  !  Interpolates at points  2,3,6  (above diag)
  do p = 1, nbfaces
       do g = 1, ng
          i = ie+g
          h = g-1
          do j = js-h,js-1
             E_buffer(i, j, p) = 0.d0
             do n = 1, L%order
                E_buffer(i, j, p) = E_buffer(i, j, p) + E_buffer_local(i, stencil_start(j,g)+n-1, p)*poly(j,g,n)
             enddo
           S_buffer(i,j,p) = E_buffer(i,j,p)
          enddo
       enddo
  enddo

  !  Interpolates at points 4,7,8     (below diag)
  do p = 1, nbfaces
     do g = 1, ng
        g2 = ng-g+1
        j  = js-g2
        h  = ng-g
        do i = ie+1,ie+h
           S_buffer(i, j, p) = 0.d0
           do n = 1, L%order
              S_buffer(i, j, p) = S_buffer(i, j, p) + S_buffer_local(stencil_start(i,g2)+n-1, j, p)*poly(i,g2,n)
           enddo
           E_buffer(i, j, p) = S_buffer(i, j, p)
        enddo
     enddo
  enddo
  !---------------------------------------------------------------------------------------------------------------------------
  !  NE corner (diag = 3,5,7)
  !        1--2--3 |
  !        |     | |
  !        4  5  6 |
  !        |     | |
  !        7--8--9 |
  !  ---------------------------------------------------
  !---------------------------------------------------------------------------------------------------------------------------
  !  Interpolates at point 8, 9, 6  (below diag)   
  do p = 1, nbfaces
     do g = 1, ng
        i = ie+g
        h = g-1
        do j = je+1, je+h
           E_buffer(i, j, p) = 0.d0
           do n = 1, L%order
              E_buffer(i, j, p) = E_buffer(i, j, p) + E_buffer_local(i, stencil_start(j,g)+n-1, p)*poly(j,g,n)
           enddo
           N_buffer(i,j,p) = E_buffer(i,j,p)
        enddo
     enddo
  enddo

  !  Interpolates at point 1, 2, 4  (above diag)   
  do p = 1, nbfaces
     do g = 1, ng
        j = je+g
        h = g-1
        do i = ie+1, ie+h
           N_buffer(i, j, p) = 0.d0
           do n = 1, L%order
              N_buffer(i, j, p) = N_buffer(i, j, p) + N_buffer_local(stencil_start(i,g)+n-1, j, p)*poly(i,g,n)
           enddo
           E_buffer(i,j,p) = N_buffer(i,j,p)
        enddo
     enddo
  enddo


  !---------------------------------------------------------------------------------------------------------------------------
  !  NW corner (diag = 1, 5, 9)
  !        1--2--3 |
  !        |     | |
  !        4  5  6 |
  !        |     | |
  !        7--8--9 |
  !        ---------
  !  ---------------------------------------------------
  !---------------------------------------------------------------------------------------------------------------------------
  !  Interpolates at points 7, 8, 4   (below diag)   
  do p = 1, nbfaces
     do g = 1, ng
        g2 = ng-g+1
        i  = is-g2
        do j = je+1,je+h
           W_buffer(i, j, p) = 0.d0
           do n = 1, L%order
              W_buffer(i, j, p) = W_buffer(i, j, p) + W_buffer_local(i, stencil_start(j,g2)+n-1,p)*poly(j,g2,n)
           enddo
           N_buffer(i, j, p) = W_buffer(i, j, p)
        enddo
     enddo
  enddo

  !  Interpolates at points  3, 6, 2   (above diag)
  do p = 1, nbfaces
     do g = 1, ng
        j = je+g
        h = g-1
        do i = is-h, is-1
           N_buffer(i, j, p) = 0.d0
           do n = 1, L%order
              N_buffer(i, j, p) = N_buffer(i, j, p) + N_buffer_local(stencil_start(i,g)+n-1, j, p)*poly(i,g,n)
           enddo
           W_buffer(i,j,p) = N_buffer(i,j,p)
        enddo
     enddo
  enddo



  !---------------------------------------------------------------------------------------------------------------------------
  !  SW corner (diag = 3,5,7)
  !        1--2--3 |
  !        |     | |
  !        4  5  6 |
  !        |     | |
  !        7--8--9 |
  !        ---------
  !---------------------------------------------------------------------------------------------------------------------------
  !  Interpolates at points 1,2,4     (above diag)
  do p = 1, nbfaces
     do g = 1, ng
        g2 = ng-g+1
        i  = is-g2
        do j = js-h,js-1
           W_buffer(i, j, p) = 0.d0
           do n = 1, L%order
              W_buffer(i, j, p) = W_buffer(i, j, p) + W_buffer_local(i, stencil_start(j,g2)+n-1,p)*poly(j,g2,n)
           enddo
           S_buffer(i, j, p) = W_buffer(i, j, p)
        enddo
     enddo
  enddo

  !  Interpolates at points 8,6,9     (below diag)
  do p = 1, nbfaces
     do g = 1, ng
        g2 = ng-g+1
        j  = js-g2
        h  = ng-g
        do i = is-h, is-1
           S_buffer(i, j, p) = 0.d0
           do n = 1, L%order
              S_buffer(i, j, p) = S_buffer(i, j, p) + S_buffer_local(stencil_start(i,g2)+n-1, j, p)*poly(i,g2,n)
           enddo
           W_buffer(i, j, p) = S_buffer(i, j, p)
        enddo
     enddo
  enddo

end subroutine cube_rmp_corner_offdiag

subroutine cube_rmp_corner_diag(bd, L, N_buffer, S_buffer, E_buffer, W_buffer, ng)
  type(fv_grid_bounds_type), intent(IN) :: bd
  type(lagrange_poly), target, intent(inout) :: L
  integer, intent(in) :: ng
  real(R_GRID), intent(INOUT) :: E_buffer(bd%ie+1:bd%ied2, bd%jsd2:bd%jed2, 1:nbfaces)
  real(R_GRID), intent(INOUT) :: W_buffer(bd%isd2:bd%is-1, bd%jsd2:bd%jed2, 1:nbfaces)
  real(R_GRID), intent(INOUT) :: N_buffer(bd%isd2:bd%ied2, bd%je+1:bd%jed2, 1:nbfaces)
  real(R_GRID), intent(INOUT) :: S_buffer(bd%isd2:bd%ied2, bd%jsd2:bd%js-1, 1:nbfaces)

  real(R_GRID), pointer, dimension(:, :, :) :: N_buffer_local
  real(R_GRID), pointer, dimension(:, :, :) :: S_buffer_local
  real(R_GRID), pointer, dimension(:, :, :) :: E_buffer_local
  real(R_GRID), pointer, dimension(:, :, :) :: W_buffer_local
  real(kind=R_GRID) :: interp_buff(1:5), interp_buff2(1:5)
  integer :: is, ie, js, je, isd, ied, jsd, jed
  integer :: i, j, ii, n, g, g2, d, h, k, p

  !--- assign parameters
  is = bd%is
  ie = bd%ie
  js = bd%js
  je = bd%je
  isd = bd%isd
  ied = bd%ied
  jsd = bd%jsd
  jed = bd%jed

  !---------------------------------------------------------------------------------------------------------------------------
  ! SE corner
  do p = 1, nbfaces
     ! First we fill the east points that will be used to interpolate the corner points
     do g = 1, ng+1
        interp_buff(g) = 0.d0
        do d = 1, 4
           interp_buff(g) = interp_buff(g) + L%E_buffer_local(ie+d,js,p)*L%p_corner1(g,d)
        enddo
     enddo

     ! now we are ready to fill the east/south corner points
     do g = 1, ng
        interp_buff2(1:3) = L%S_buffer_local(ie-2:ie,js-g,p)
        interp_buff2(4)   = interp_buff(g)
        E_buffer(ie+g,js-g,p) = 0.d0
        do d = 1, 4
           E_buffer(ie+g,js-g,p) = E_buffer(ie+g,js-g,p) + 0.5d0*interp_buff2(d)*L%p_corner2(g,d)
        enddo
     enddo

     ! First we fill the south points that will be used to interpolate the corner points
     do g = 1, ng+1
        interp_buff(g) = 0.d0
        do d = 1, 4
          interp_buff(g) = interp_buff(g) + L%S_buffer_local(ie,js-d,p)*L%p_corner1(g,d)
        enddo
     enddo
  
     ! now we are ready to fill the east/south corner points
     do g = 1, ng
        interp_buff2(1) = L%E_buffer_local(ie+g,js+2,p)
        interp_buff2(2) = L%E_buffer_local(ie+g,js+1,p)
        interp_buff2(3) = L%E_buffer_local(ie+g,js,p)
        interp_buff2(4) = interp_buff(g)
        do d = 1, 4
           E_buffer(ie+g,js-g,p) = E_buffer(ie+g,js-g,p) + 0.5d0*interp_buff2(d)*L%p_corner2(g,d)
        enddo
        S_buffer(ie+g,js-g,p) = E_buffer(ie+g,js-g,p)
     enddo
  
  enddo


  !---------------------------------------------------------------------------------------------------------------------------
  ! NE corner
  do p = 1, nbfaces
     ! First we fill the east points that will be used to interpolate the corner points
     do g = 1, ng+1
        interp_buff(g) = 0.d0
        do d = 1, 4
          interp_buff(g) = interp_buff(g) + L%E_buffer_local(ie+d,je,p)*L%p_corner1(g,d)
        enddo
     enddo

     ! now we are ready to fill the east/north corner points
     do g = 1, ng
        interp_buff2(1:3) = L%N_buffer_local(ie-2:ie,je+g,p)
        interp_buff2(4)   = interp_buff(g)
        E_buffer(ie+g,je+g,p) = 0.d0
        do d = 1, 4
           E_buffer(ie+g,je+g,p) = E_buffer(ie+g,je+g,p) + 0.5d0*interp_buff2(d)*L%p_corner2(g,d)
        enddo
     enddo

     ! First we fill the north points that will be used to interpolate the corner points
     do g = 1, ng+1
        interp_buff(g) = 0.d0
        do d = 1, 4
          interp_buff(g) = interp_buff(g) + L%N_buffer_local(ie,je+d,p)*L%p_corner1(g,d)
        enddo
     enddo
    
     ! now we are ready to fill the east/north corner points
     do g = 1, ng
        interp_buff2(1) = L%E_buffer_local(ie+g,je-2,p)
        interp_buff2(2) = L%E_buffer_local(ie+g,je-1,p)
        interp_buff2(3) = L%E_buffer_local(ie+g,je,p)
        interp_buff2(4) = interp_buff(g)
        do d = 1, 4
           E_buffer(ie+g,je+g,p) = E_buffer(ie+g,je+g,p) + 0.5d0*interp_buff2(d)*L%p_corner2(g,d)
        enddo
        N_buffer(ie+g,je+g,p) = E_buffer(ie+g,je+g,p)
     enddo
  enddo

  !---------------------------------------------------------------------------------------------------------------------------
  ! NW corner
  do p = 1, nbfaces
      ! First we fill the west points that will be used to interpolate the corner points
      do g = 1, ng+1
         interp_buff(g) = 0.d0
         do d = 1, 4
            interp_buff(g) = interp_buff(g) + L%W_buffer_local(is-d,je,p)*L%p_corner1(g,d)
         enddo
      enddo

      ! now we are ready to fill the west/north corner points
      do g = 1, ng
         interp_buff2(1) = L%N_buffer_local(is+2,je+g,p)
         interp_buff2(2) = L%N_buffer_local(is+1,je+g,p)
         interp_buff2(3) = L%N_buffer_local(is  ,je+g,p)
         interp_buff2(4) = interp_buff(g)
         W_buffer(is-g,je+g,p) = 0.d0
         do d = 1, 4
            W_buffer(is-g,je+g,p) = W_buffer(is-g,je+g,p) + 0.5d0*interp_buff2(d)*L%p_corner2(g,d)
         enddo
      enddo

      ! First we fill the north points that will be used to interpolate the corner points
      do g = 1, ng+1
         interp_buff(g) = 0.d0
         do d = 1, 4
           interp_buff(g) = interp_buff(g) + L%N_buffer_local(is,je+d,p)*L%p_corner1(g,d)
         enddo
      enddo
     
      ! now we are ready to fill the east/north corner points
      do g = 1, ng
         interp_buff2(1) = L%W_buffer_local(is-g,je-2,p)
         interp_buff2(2) = L%W_buffer_local(is-g,je-1,p)
         interp_buff2(3) = L%W_buffer_local(is-g,je,p)
         interp_buff2(4) = interp_buff(g)
         do d = 1, 4
            W_buffer(is-g,je+g,p) = W_buffer(is-g,je+g,p) + 0.5d0*interp_buff2(d)*L%p_corner2(g,d)
         enddo
         N_buffer(is-g,je+g,p) = W_buffer(is-g,je+g,p)
      enddo
  enddo


  !---------------------------------------------------------------------------------------------------------------------------
  ! SW corner
  do p = 1, nbfaces
      ! First we fill the south points that will be used to interpolate the corner points
      do g = 1, ng+1
         interp_buff(g) = 0.d0
         do d = 1, 4
            interp_buff(g) = interp_buff(g) + L%W_buffer_local(is-d,js,p)*L%p_corner1(g,d)
         enddo
      enddo

      ! now we are ready to fill the west/south corner points
      do g = 1, ng
         interp_buff2(1) = L%S_buffer_local(is+2,js-g,p)
         interp_buff2(2) = L%S_buffer_local(is+1,js-g,p)
         interp_buff2(3) = L%S_buffer_local(is  ,js-g,p)
         interp_buff2(4) = interp_buff(g)
         W_buffer(is-g,js-g,p) = 0.d0
         do d = 1, 4
            W_buffer(is-g,js-g,p) = W_buffer(is-g,js-g,p) + 0.5d0*interp_buff2(d)*L%p_corner2(g,d)
         enddo
      enddo

      ! First we fill the north points that will be used to interpolate the corner points
      do g = 1, ng+1
         interp_buff(g) = 0.d0
         do d = 1, 4
           interp_buff(g) = interp_buff(g) + L%S_buffer_local(is,js-d,p)*L%p_corner1(g,d)
         enddo
      enddo
     
      ! now we are ready to fill the east/north corner points
      do g = 1, ng
         interp_buff2(1) = L%W_buffer_local(is-g,js+2,p)
         interp_buff2(2) = L%W_buffer_local(is-g,js+1,p)
         interp_buff2(3) = L%W_buffer_local(is-g,js,p)
         interp_buff2(4) = interp_buff(g)
         do d = 1, 4
            W_buffer(is-g,js-g,p) = W_buffer(is-g,js-g,p) + 0.5d0*interp_buff2(d)*L%p_corner2(g,d)
         enddo
         S_buffer(is-g,js-g,p) = W_buffer(is-g,js-g,p)
      enddo
  
  enddo 
end subroutine cube_rmp_corner_diag

subroutine init_lagrange(L, bd, grid_type)
    type(lagrange_poly), intent(inout):: L
    type(fv_grid_bounds_type), intent(INOUT) :: bd
    integer, intent(in) :: grid_type
    integer :: is, ie, js, je, ng2
    integer :: isd2, ied2, jsd2, jed2
    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je

    isd2 = bd%isd2
    ied2 = bd%ied2
    jsd2 = bd%jsd2
    jed2 = bd%jed2

    ng2 = bd%ng+2


    ! Lagrange polynomials
    L%degree = 3!dg%k2e_nord-1
    L%order  = L%degree + 1
    L%stencil_offset = ceiling(L%order*0.5d0)

    if (grid_type==0) then !equiedge
       L%aref = dasin(1.d0/dsqrt(3.d0))
       L%Rref = dsqrt(2.d0)
    else if (grid_type==2) then !equiangular
       L%aref = pi*0.25d0
       L%Rref = 1.d0
    endif
    L%dy = 2.d0*L%aref/bd%npx

    ! Agrid_p1 points and data
    allocate(L%y_Agrid_given(jsd2:jed2))
    allocate(L%angles_a(jsd2:jed2))
    allocate(L%tan_angles_a(jsd2:jed2))
    allocate(L%p_corner1(1:ng2, 1:4))
    allocate(L%p_corner2(1:ng2, 1:4))

    ! Agrid_p2 where we want to interpolate 
    allocate(L%y_Agrid_target(isd2:ied2, 1:ng2))
    allocate(L%p_Agrid_ydir(jsd2:jed2, 1:ng2, 1:L%order))

    ! a grid buffers
    allocate(L%E_buffer(ie+1:ied2, jsd2:jed2, 1:nbfaces))
    allocate(L%W_buffer(isd2:is-1, jsd2:jed2, 1:nbfaces))
    allocate(L%N_buffer(isd2:ied2, je+1:jed2, 1:nbfaces))
    allocate(L%S_buffer(isd2:ied2, jsd2:js-1, 1:nbfaces))

    allocate(L%E_buffer_local(ie+1:ied2, jsd2:jed2, 1:nbfaces))
    allocate(L%W_buffer_local(isd2:is-1, jsd2:jed2, 1:nbfaces))
    allocate(L%N_buffer_local(isd2:ied2, je+1:jed2, 1:nbfaces))
    allocate(L%S_buffer_local(isd2:ied2, jsd2:js-1, 1:nbfaces))

    ! stencils - agrid
    allocate(L%stencil_start_a(isd2:ied2,1:ng2))
    allocate(L%stencil_end_a  (isd2:ied2,1:ng2))

    call compute_lagrange_polynomials_agrid(L, bd)
end subroutine init_lagrange

subroutine compute_lagrange_polynomials_agrid(L, bd)
    !---------------------------------------------------
    ! compute the lagrange polynomials at ghost cell
    ! points of the cubed-sphere
    ! A grid
    !--------------------------------------------------
    type(lagrange_poly), intent(inout):: L
    type(fv_grid_bounds_type), intent(INOUT) :: bd
    real(kind=R_GRID) :: pa(1:3), x_target_p1, y_target_p1, x_target_p2, y_target_p2
    real(kind=R_GRID) :: x_given_p2, y_given_p2, x_target, y_target, x_target2
    real(kind=R_GRID) :: x_given_p6, y_given_p6
    real(kind=R_GRID) :: p_given_p2(1:3), ll_given_p2(2)
    real(kind=R_GRID) :: support_points(1:L%order), lat, lon, c(1:2), angl_g, y, angl_y
    real(kind=R_GRID) :: support_points_corner1(1:4)
    real(kind=R_GRID) :: support_points_corner2(1:4)
    real(kind=R_GRID) :: support_points_corner3(1:L%order)
    integer :: gridtype = 2!0 or 2
    integer ::  i, j, ii, jj, p, g, d, k, h, jnearest, panel1, panel2, panel6
    integer ::  is, ie, js, je, ng
    integer ::  isd, ied, jsd, jed

    is = bd%is
    ie = bd%ie
    js = bd%js
    je = bd%je

    isd = bd%isd-2
    ied = bd%ied+2
    jsd = bd%jsd-2
    jed = bd%jed+2

    ng = bd%ng+2
    if(gridtype==0 .or. gridtype==2) then
        ! Compute the Agrid in panel 2 (y coordinates only are needed)
        ! these are the points where the data required by interpolation is given.
        do j = jsd, jed
           L%y_Agrid_given(j) = -L%aref + (j-0.5d0)*L%dy
           L%angles_a(j) = L%y_Agrid_given(j) 
           L%tan_angles_a(j) = dtan(L%angles_a(j))*L%Rref
        end do

        !modify ghost cells at left
        do j = jsd, js - 1
            jj  = 2*js-j-1
            L%y_Agrid_given(j) = -0.5*pi - datan(L%tan_angles_a(jj))
            L%angles_a(j) = L%y_Agrid_given(j)
            L%tan_angles_a(j) = dtan(L%angles_a(j))
        end do

        do g = 1, ng
           L%y_Agrid_given(je+g) = -L%y_Agrid_given(1-g) 
           L%angles_a(je+g)      = -L%angles_a(1-g) 
           L%tan_angles_a(je+g)  = -L%tan_angles_a(1-g) 
        end do

        ! Invert node points using the cs mapping from panel 2
        panel1 = 1
        panel2 = 2
        do g = 1, ng
           x_target_p1 = L%tan_angles_a(je+g)
           do j = jsd, jed
              y_target_p1 = L%tan_angles_a(j)

              ! project target point (x_target, y_target) point using panel 1 mapping
              call equidistant_gnomonic_map(pa, x_target_p1, y_target_p1, panel1)

              ! invert the target point using panel 2 mapping to get its (x,y) coordinates in panel 2 system
              call inverse_equidistant_gnomonic_map(pa, x_target_p2, y_target_p2, panel2)
              x_target_p2 = datan(x_target_p2/L%Rref)
              y_target_p2 = datan(y_target_p2/L%Rref)
              L%y_Agrid_target(j,g) = y_target_p2
           end do
        end do

        L%stencil_start_a(:,:) = 0
        L%stencil_end_a(:,:) = 0

        ! Compute stencils
        do g = 1, ng
           do j = js-g+1, je+g-1
              do i = jsd, jed-1
                 if( (L%y_Agrid_given(i) <= L%y_Agrid_target(j,g)) .and. (L%y_Agrid_target(j,g) <= L%y_Agrid_given(i+1) )) then
                    exit
                 end if
              enddo
             
              L%stencil_end_a(j,g)   = i+L%stencil_offset
              L%stencil_start_a(j,g) = L%stencil_end_a(j,g) - L%order + 1

              if(L%stencil_start_a(j,g)<js)then 
                 L%stencil_start_a(j,g) = js
                 L%stencil_end_a(j,g) = L%stencil_start_a(j,g) + L%order - 1
              else if(L%stencil_end_a(j,g)>je)then 
                 L%stencil_end_a(j,g) = je
                 L%stencil_start_a(j,g) = L%stencil_end_a(j,g) - L%order + 1
              end if
           enddo
        enddo

        ! debug to see if the stencils make sense 
        do g = 1, ng
            do j = js-g+1, je+g-1
                if(L%stencil_end_a(j,g)-L%stencil_start_a(j,g) .ne. L%degree) then
                  print*, 'ERROR in compute_lagrange_polynomials_agrid: stencil size is not correct.'
                  stop
                endif

                if(L%stencil_end_a(j,g)<js .or. L%stencil_start_a(j,g)>je) then
                  print*, 'ERROR in compute_lagrange_polynomials_agrid: stencil bounds are not correct.'
                  stop
                endif
 
                if(L%y_Agrid_given(L%stencil_start_a(j,g)) > L%y_Agrid_target(j,g)) then
                  print*, 'ERROR in compute_lagrange_polynomials_agrid: error in neighboring points. >'
                  stop
                endif

                if(L%y_Agrid_given(L%stencil_end_a(j,g))   < L%y_Agrid_target(j,g)) then
                  print*, 'ERROR in compute_lagrange_polynomials_agrid: error in neighboring points. <'
                  stop
                endif
            enddo
        enddo
        ! Compute the Lagrange Agrid_p2 at halo region
        do g = 1, ng
           do j = js-g+1, je+g-1
              support_points(:) = L%y_Agrid_given(L%stencil_start_a(j,g):L%stencil_end_a(j,g))
              do d = 1, L%order
                 call lagrange_basis(L%y_Agrid_target(j,g), support_points(:), L%degree, d, L%p_Agrid_ydir(j,g,d))
              enddo
           enddo
        enddo

        ! extra polynomials for diagonal corner points
        support_points_corner1(1:4) = L%y_Agrid_given(js:js+3)
        panel6 = 6
        do g = 1, ng
           x_given_p6 =  L%tan_angles_a(je+1)
           y_given_p6 =  L%tan_angles_a(je-g+1)

           ! project target point using panel 6 mapping
           call equidistant_gnomonic_map(pa, x_given_p6, y_given_p6, panel6)

           ! invert the target point using panel 2 mapping to get its (x,y) coordinates in panel 2 system
           call inverse_equidistant_gnomonic_map(pa, x_given_p2, y_given_p2, panel2)

           x_given_p2 = datan(x_given_p2/L%Rref)
           y_given_p2 = datan(y_given_p2/L%Rref)
 
           do d = 1, 4 ! cubic polynomials will be used here
             call lagrange_basis(x_given_p2, support_points_corner1(:), 3, d, L%p_corner1(g,d))
           enddo
        enddo

        support_points_corner2(1:3) = L%y_Agrid_given(je-2:je)
        x_target = L%aref
        do g = 1, ng-1
           x_given_p6 =  L%tan_angles_a(je+1)
           y_given_p6 =  L%tan_angles_a(je-g)
           x_given_p2 = datan(L%tan_angles_a(je+1)/L%Rref)
           support_points_corner2(4) = x_given_p2

           do d = 1, 4 ! cubic polynomials will be used here
              call lagrange_basis(x_target, support_points_corner2(:), 3, d, L%p_corner2(g,d))
           enddo
        enddo
    else
      print*, 'ERROR in compute_lagrange_polynomials_agrid: invalid gridtype.'
      stop
    endif

end subroutine compute_lagrange_polynomials_agrid

subroutine lagrange_basis(x, x_Agrid_p1, N, j, Lj)
    !---------------------------------------------------
    ! Compute the jth Lagrange polynomial of degree N
    !--------------------------------------------------
    real(kind=R_GRID), intent(in) :: x, x_Agrid_p1(1:N+1)
    real(kind=R_GRID), intent(inout) :: Lj
    integer, intent(in) :: N, j
    integer :: i

    Lj = 1.d0
    do i = 1, N+1
       if (i .ne. j) then
          Lj = Lj*(x-x_Agrid_p1(i))/(x_Agrid_p1(j)-x_Agrid_p1(i))
       end if
    end do
end subroutine lagrange_basis

end module fv_duogrid
