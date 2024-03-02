module tp_core
!========================================================================
! This module contains the routine that computes the PPM flux
!
! reference https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere/blob/main/model/tp_core.F90
!========================================================================
use fv_arrays, only: fv_grid_type, fv_grid_bounds_type, R_GRID, nbfaces
 implicit none

 private
 public fv_tp_2d

 real, parameter:: ppm_fac = 1.5   ! nonlinear scheme limiter: between 1 and 2
 real, parameter:: r3 = 1./3.
 real, parameter:: near_zero = 1.E-25
 real, parameter:: ppm_limiter = 2.0
 real, parameter:: r12 = 1./12.

! Suresh & Huynh scheme 2.2 (purtabation form)
! The wave-form is more diffusive than scheme 2.1
 real, parameter:: b1 =   0.0375
 real, parameter:: b2 =  -7./30.
 real, parameter:: b3 =  -23./120.
 real, parameter:: b4 =  13./30.
 real, parameter:: b5 = -11./240.
! scheme 2.1: perturbation form
 !real, parameter:: b1 =   1./30.
 !real, parameter:: b2 = -13./60.
 !real, parameter:: b3 = -13./60.
 !real, parameter:: b4 =  0.45
 !real, parameter:: b5 = -0.05

 real, parameter:: t11 = 27./28., t12 = -13./28., t13=3./7.
 real, parameter:: s11 = 11./14., s14 = 4./7.,    s15=3./14.
!----------------------------------------------------
! volume-conserving cubic with 2nd drv=0 at end point:
!----------------------------------------------------
! Non-monotonic
  real, parameter:: c1 = -2./14.
  real, parameter:: c2 = 11./14.
  real, parameter:: c3 =  5./14.
!----------------------
! PPM volume mean form:
!----------------------
  real, parameter:: p1 =  7./12.     ! 0.58333333
  real, parameter:: p2 = -1./12.
!   q(i+0.5) = p1*(q(i-1)+q(i)) + p2*(q(i-2)+q(i+1))
! integer:: is, ie, js, je, isd, ied, jsd, jed

!
!-----------------------------------------------------------------------
contains
 subroutine fv_tp_2d(q, crx, cry, hord, fx, fy, xfx, yfx, &
                     gridstruct, bd, ra_x, ra_y, lim_fac, inner_adv)
   type(fv_grid_bounds_type), intent(IN) :: bd
   integer, intent(in)::hord
   integer, intent(in)::inner_adv

   real(R_GRID), intent(in)::  crx(bd%is:bd%ie+1,bd%jsd:bd%jed, 1:nbfaces)  !
   real(R_GRID), intent(in)::  cry(bd%isd:bd%ied,bd%js:bd%je+1, 1:nbfaces )  !

   real(R_GRID), intent(in)::  xfx(bd%is:bd%ie+1,bd%jsd:bd%jed, 1:nbfaces)  !
   real(R_GRID), intent(in)::  yfx(bd%isd:bd%ied,bd%js:bd%je+1, 1:nbfaces)  !

   real(R_GRID), intent(in):: ra_x(bd%is:bd%ie,bd%jsd:bd%jed, 1:nbfaces)
   real(R_GRID), intent(in):: ra_y(bd%isd:bd%ied,bd%js:bd%je, 1:nbfaces)
   real(R_GRID), intent(inout):: q(bd%isd:bd%ied,bd%jsd:bd%jed, 1:nbfaces)  ! transported scalar
   real(R_GRID), intent(out)::fx(bd%is:bd%ie+1 ,bd%js:bd%je, 1:nbfaces)    ! Flux in x ( E )
   real(R_GRID), intent(out)::fy(bd%is:bd%ie,   bd%js:bd%je+1, 1:nbfaces)    ! Flux in y ( N )

   type(fv_grid_type), intent(IN), target :: gridstruct

   real(R_GRID), intent(in):: lim_fac
! optional Arguments:
! Local:
   real(R_GRID):: q_i(bd%isd:bd%ied,bd%js:bd%je, 1:nbfaces)
   real(R_GRID):: q_j(bd%is:bd%ie,bd%jsd:bd%jed, 1:nbfaces)
   real(R_GRID):: fxx(bd%is:bd%ie+1,bd%jsd:bd%jed, 1:nbfaces)
   real(R_GRID):: fx2(bd%is:bd%ie+1,bd%jsd:bd%jed, 1:nbfaces)
   real(R_GRID):: fyy(bd%isd:bd%ied,bd%js:bd%je+1, 1:nbfaces)
   real(R_GRID):: fy2(bd%isd:bd%ied,bd%js:bd%je+1, 1:nbfaces)
   real(R_GRID):: qmt(bd%isd:bd%ied,bd%jsd:bd%jed, 1:nbfaces)  ! transported scalar * metricterm
   real(R_GRID), pointer, dimension(:,:) ::  mt_a
   integer i, j, p
   integer:: is, ie, js, je, isd, ied, jsd, jed

   is  = bd%is
   ie  = bd%ie
   isd = bd%isd
   ied = bd%ied

   js  = bd%js
   je  = bd%je
   jsd = bd%jsd
   jed = bd%jed

   mt_a  => gridstruct%mt
   !==================================================================================================================

   if(inner_adv==1) then
      !do p = 1, nbfaces
      !$OMP PARALLEL WORKSHARE DEFAULT(NONE) &
      !$OMP SHARED(qmt, q, isd, ied, jsd, jed)
      qmt(isd:ied,jsd:jed,:) = q(isd:ied,jsd:jed,:)
      !$OMP END PARALLEL WORKSHARE
   else if(inner_adv==2) then
      !$OMP PARALLEL DO &
      !$OMP DEFAULT(NONE) & 
      !$OMP SHARED(qmt, mt_a, q) &
      !$OMP SHARED(isd, ied, jsd, jed, nbfaces) &
      !$OMP SCHEDULE(static)
      do p = 1, nbfaces
         qmt(isd:ied,jsd:jed,p) = q(isd:ied,jsd:jed,p)*mt_a(isd:ied,jsd:jed)
      enddo
      !$OMP END PARALLEL DO
   endif

   call yppm(fy2, qmt, cry, hord, isd,ied,isd,ied, js,je,jsd,jed, lim_fac)

   !$OMP PARALLEL DO &
   !$OMP DEFAULT(NONE) & 
   !$OMP SHARED(fyy, yfx, fy2) &
   !$OMP SHARED(is, ie, js, je, nbfaces) &
   !$OMP SHARED(isd, ied, jsd, jed) &
   !$OMP PRIVATE(i, j) &
   !$OMP SCHEDULE(static)
   do p = 1, nbfaces
      do j=js,je+1
         do i=isd,ied
            fyy(i,j,p) = yfx(i,j,p) * fy2(i,j,p)
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO

   ! inner advection
   if(inner_adv==1) then
      !$OMP PARALLEL DO &
      !$OMP DEFAULT(NONE) & 
      !$OMP SHARED(q_i, q, gridstruct, fyy, ra_y) &
      !$OMP SHARED(is, ie, js, je, nbfaces) &
      !$OMP SHARED(isd, ied, jsd, jed) &
      !$OMP PRIVATE(i, j) &
      !$OMP SCHEDULE(static)
      do p = 1, nbfaces
         do j=js,je
            do i=isd,ied
               q_i(i,j,p) = (q(i,j,p)*gridstruct%area(i,j) + fyy(i,j,p)-fyy(i,j+1,p))/ra_y(i,j,p)
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
   else if(inner_adv==2) then
      !$OMP PARALLEL DO &
      !$OMP DEFAULT(NONE) & 
      !$OMP SHARED(q_i, q, gridstruct, fyy) &
      !$OMP SHARED(is, ie, js, je, nbfaces) &
      !$OMP SHARED(isd, ied, jsd, jed) &
      !$OMP PRIVATE(i, j) &
      !$OMP SCHEDULE(static)
      do p = 1, nbfaces
         do j=js,je
            do i=isd,ied
               q_i(i,j,p) = q(i,j,p) -(fyy(i,j+1,p)-fyy(i,j,p))*gridstruct%rarea(i,j)
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
   endif

   if(inner_adv==2) then
      !$OMP PARALLEL DO &
      !$OMP DEFAULT(NONE) & 
      !$OMP SHARED(q_i, mt_a) &
      !$OMP SHARED(is, ie, js, je, nbfaces) &
      !$OMP SHARED(isd, ied, jsd, jed) &
      !$OMP PRIVATE(i, j) &
      !$OMP SCHEDULE(static)
      do p = 1, nbfaces
         do j=js,je
            do i=isd,ied
               q_i(i,j,p) = q_i(i,j,p)*mt_a(i,j)
            enddo
         enddo
      enddo 
      !$OMP END PARALLEL DO
   endif


   call xppm(fx, q_i, crx, hord, is,ie,isd,ied, js,je,jsd,jed, lim_fac)


   !==================================================================================================================
   call xppm(fx2, qmt, crx, hord, is,ie,isd,ied, jsd,jed,jsd,jed, lim_fac)

   !$OMP PARALLEL DO &
   !$OMP DEFAULT(NONE) & 
   !$OMP SHARED(fxx, xfx, fx2) &
   !$OMP SHARED(is, ie, js, je, nbfaces) &
   !$OMP SHARED(isd, ied, jsd, jed) &
   !$OMP PRIVATE(i, j) &
   !$OMP SCHEDULE(static)
   do p = 1, nbfaces
      do j=jsd,jed
         do i=is,ie+1
            fxx(i,j,p) =  xfx(i,j,p) * fx2(i,j,p)
         enddo
      enddo 
   enddo 
   !$OMP END PARALLEL DO

   ! inner adv
   if(inner_adv==1) then
      !$OMP PARALLEL DO &
      !$OMP DEFAULT(NONE) & 
      !$OMP SHARED(q_j, q, gridstruct, fxx, ra_x) &
      !$OMP SHARED(is, ie, js, je, nbfaces) &
      !$OMP SHARED(isd, ied, jsd, jed) &
      !$OMP PRIVATE(i, j) &
      !$OMP SCHEDULE(static)
      do p = 1, nbfaces
         do j=jsd,jed
            do i=is,ie
               q_j(i,j,p) = (q(i,j,p)*gridstruct%area(i,j) + fxx(i,j,p)-fxx(i+1,j,p))/ra_x(i,j,p)
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
   else if(inner_adv==2) then
      !$OMP PARALLEL DO &
      !$OMP DEFAULT(NONE) & 
      !$OMP SHARED(q_j, q, gridstruct, fxx) &
      !$OMP SHARED(is, ie, js, je, nbfaces) &
      !$OMP SHARED(isd, ied, jsd, jed) &
      !$OMP PRIVATE(i, j) &
      !$OMP SCHEDULE(static)
      do p = 1, nbfaces
         do j=jsd,jed
            do i=is,ie
               q_j(i,j,p) = q(i,j,p)-(fxx(i+1,j,p)-fxx(i,j,p))*gridstruct%rarea(i,j)
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
   endif 
 
   if(inner_adv==2) then
      !$OMP PARALLEL DO &
      !$OMP DEFAULT(NONE) & 
      !$OMP SHARED(q_j, mt_a) &
      !$OMP SHARED(is, ie, js, je, nbfaces) &
      !$OMP SHARED(isd, ied, jsd, jed) &
      !$OMP PRIVATE(i, j) &
      !$OMP SCHEDULE(static)
      do p = 1, nbfaces
         do j=jsd,jed
            do i=is,ie
               q_j(i,j,p) = q_j(i,j,p)*mt_a(i,j)
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO
   endif

   call yppm(fy, q_j, cry, hord, is,ie,isd,ied, js,je,jsd,jed, lim_fac)


   !==================================================================================================================
   !$OMP PARALLEL DO &
   !$OMP DEFAULT(NONE) & 
   !$OMP SHARED(fx, fy, fx2, fy2, xfx, yfx) &
   !$OMP SHARED(is, ie, js, je, nbfaces) &
   !$OMP SHARED(isd, ied, jsd, jed) &
   !$OMP PRIVATE(i, j) &
   !$OMP SCHEDULE(static)
   do p = 1, nbfaces
      do j=js,je
         do i=is,ie+1
            fx(i,j,p) = 0.5*(fx(i,j,p)*xfx(i,j,p) + fx2(i,j,p)*xfx(i,j,p))
         enddo
      enddo
      do j=js,je+1
         do i=is,ie
            fy(i,j,p) = 0.5*(fy(i,j,p)*yfx(i,j,p) + fy2(i,j,p)*yfx(i,j,p))
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO
 end subroutine fv_tp_2d

 subroutine xppm(flux, q, c, iord, is,ie,isd,ied, jfirst,jlast,jsd,jed, lim_fac)
 integer, INTENT(IN) :: is, ie, isd, ied, jsd, jed
 integer, INTENT(IN) :: jfirst, jlast  ! compute domain
 integer, INTENT(IN) :: iord
 real(R_GRID), INTENT(IN) :: q(isd:ied,jfirst:jlast,1:nbfaces)
 real(R_GRID), INTENT(IN) :: c(is:ie+1,jsd:jed,1:nbfaces) ! Courant   N (like FLUX)
 real(R_GRID)   , intent(IN) :: lim_fac
! !OUTPUT PARAMETERS:
 real(R_GRID)  , INTENT(OUT) :: flux(is:ie+1,jfirst:jlast,1:6) !  Flux
! Local
 real(R_GRID), dimension(is-1:ie+1,jfirst:jlast):: bl, br, b0, da1
 real(R_GRID):: q1(isd:ied)
 real(R_GRID), dimension(jfirst:jlast):: fx0, fx1, xt1, a4
 logical, dimension(is-1:ie+1,jfirst:jlast):: ext5, ext6, smt5, smt6
 logical, dimension(is:ie+1):: hi5, hi6
 real(R_GRID)  al(is-1:ie+2,jfirst:jlast)
 real(R_GRID)  dm(is-2:ie+2,jfirst:jlast)
 real(R_GRID)  dq(is-3:ie+2,jfirst:jlast)
 integer:: i, j, ie3, is1, ie1, mord, p
 real(R_GRID):: x0, x1, xt, qtmp, pmp_1, lac_1, pmp_2, lac_2

 !if ( .not. bounded_domain .and. grid_type<3 ) then
 !   is1 = max(3,is-1);  ie3 = min(npx-2,ie+2)
 !                       ie1 = min(npx-3,ie+1)
 !else
    is1 = is-1;         ie3 = ie+2
                        ie1 = ie+1
 !end if

 mord = abs(iord)

do p = 1, nbfaces
 mord = abs(iord)
if (iord==0) then
   do i=is1, ie3
      do j=jfirst,jlast
         al(i,j) = p1*(q(i-1,j,p)+q(i,j,p)) + p2*(q(i-2,j,p)+q(i+1,j,p))
      enddo
   enddo

   do i=is1, ie1
      do j=jfirst,jlast
         bl(i,j) =  al(i,j  )-q(i,j,p)
         br(i,j) =  al(i+1,j)-q(i,j,p)
      enddo
   enddo


   do i=is,ie+1
      do j=jfirst,jlast
         if( c(i,j,p)>0. ) then
            flux(i,j,p) = q(i-1,j,p) + (1.-c(i,j,p))*(br(i-1,j)-c(i,j,p)*(bl(i-1,j)+br(i-1,j)))
         else
            flux(i,j,p) = q(i,j  ,p) + (1.+c(i,j,p))*(bl(i,j  )+c(i,j,p)*(bl(i,j  )+br(i,j  )))
         endif
      enddo
   enddo


else if ( iord < 7 ) then

   do j=jfirst,jlast
      do i=is1, ie3
         al(i,j) = p1*(q(i-1,j,p)+q(i,j,p)) + p2*(q(i-2,j,p)+q(i+1,j,p))
      enddo
   enddo

   !if ( .not. bounded_domain .and. grid_type<3 ) then
   !   if( is==1 ) then
   !     do j=jfirst,jlast
   !        al(i,0) = c1*q(i,-2) + c2*q(i,-1) + c3*q(i,0)
   !        al(i,1) = 0.5*(((2.*dya(i,0)+dya(i,-1))*q(i,0)-dya(i,0)*q(i,-1))/(dya(i,-1)+dya(i,0))   &
   !                +      ((2.*dya(i,1)+dya(i,2))*q(i,1)-dya(i,1)*q(i,2))/(dya(i,1)+dya(i,2)))
   !        al(i,2) = c3*q(i,1) + c2*q(i,2) + c1*q(i,3)
   !     enddo
   !   endif
   !   if( (ie+1)==npy ) then
   !     do j=jfirst,jlast
   !      al(i,npy-1) = c1*q(i,npy-3) + c2*q(i,npy-2) + c3*q(i,npy-1)
   !      al(i,npy) = 0.5*(((2.*dya(i,npy-1)+dya(i,npy-2))*q(i,npy-1)-dya(i,npy-1)*q(i,npy-2))/(dya(i,npy-2)+dya(i,npy-1))  &
   !                +      ((2.*dya(i,npy)+dya(i,npy+1))*q(i,npy)-dya(i,npy)*q(i,npy+1))/(dya(i,npy)+dya(i,npy+1)))
   !      al(i,npy+1) = c3*q(i,npy) + c2*q(i,npy+1) + c1*q(i,npy+2)
   !     enddo
   !   endif
   !endif

   if ( iord<0 ) then
      do j=jfirst,jlast
         do i=is-1, ie+2
            al(i,j) = max(0., al(i,j))
         enddo
      enddo
   endif

   if ( mord==1 ) then
       do j=jfirst,jlast
          do i=is-1,ie+1
             bl(i,j) = al(i,j  ) - q(i,j,p)
             br(i,j) = al(i+1,j) - q(i,j,p)
             b0(i,j) = bl(i,j) + br(i,j)
             smt5(i,j) = abs(lim_fac*b0(i,j)) < abs(bl(i,j)-br(i,j))
          enddo
       enddo

       do j=jfirst,jlast
!DEC$ VECTOR ALWAYSi
          do i=is,ie+1
             if ( c(i,j,p) > 0. ) then
                  fx1(j) = (1.-c(i,j,p))*(br(i-1,j) - c(i,j,p)*b0(i-1,j))
                  flux(i,j,p) = q(i-1,j,p)
             else
                  fx1(j) = (1.+c(i,j,p))*(bl(i,j) + c(i,j,p)*b0(i,j))
                  flux(i,j,p) = q(i,j,p)
             endif
             if (smt5(i-1,j).or.smt5(i,j)) flux(i,j,p) = flux(i,j,p) + fx1(j)
          enddo
       enddo

   elseif ( mord==2 ) then   ! Perfectly linear scheme
! Diffusivity: ord2 < ord5 < ord3 < ord4 < ord6  < ord7

      do j=jfirst,jlast
!DEC$ VECTOR ALWAYS
         do i=is,ie+1
            xt = c(i,j,p)
            if ( xt > 0. ) then
                 qtmp = q(i-1,j,p)
                 flux(i,j,p) = qtmp + (1.-xt)*(al(i,j)-qtmp-xt*(al(i-1,j)+al(i,j)-(qtmp+qtmp)))
            else
                 qtmp = q(i,j,p)
                 flux(i,j,p) = qtmp + (1.+xt)*(al(i,j)-qtmp+xt*(al(i,j)+al(i+1,j)-(qtmp+qtmp)))
            endif
         enddo
      enddo

   elseif ( mord==3 ) then

       do j=jfirst,jlast
         do i=is-1,ie+1
              bl(i,j) = al(i,j  ) - q(i,j,p)
              br(i,j) = al(i+1,j) - q(i,j,p)
              b0(i,j) = bl(i,j) + br(i,j)
                   x0 = abs(b0(i,j))
                   xt = abs(bl(i,j)-br(i,j))
              smt5(i,j) =    x0 < xt
              smt6(i,j) = 3.*x0 < xt
           enddo
        enddo
        do j=jfirst,jlast
           do i=is,ie+1
              xt1(j) = c(i,j,p)
           enddo
           do i=is,ie+1
              if ( xt1(j) > 0. ) then
                   if( smt5(i-1,j) .or. smt6(i,j) ) then
                       flux(i,j,p) = q(i-1,j,p) + (1.-xt1(j))*(br(i-1,j) - xt1(j)*b0(i-1,j))
                   else
                       flux(i,j,p) = q(i-1,j,p)
                   endif
              else
                   if( smt6(i-1,j) .or. smt5(i,j) ) then
                       flux(i,j,p) = q(i,j,p) + (1.+xt1(j))*(bl(i,j) + xt1(j)*b0(i,j))
                   else
                       flux(i,j,p) = q(i,j,p)
                   endif
              endif
           enddo
        enddo

   elseif ( mord==4 ) then

        do j=jfirst,jlast
           do i=is-1,ie+1
              bl(i,j) = al(i,j  ) - q(i,j,p)
              br(i,j) = al(i+1,j) - q(i,j,p)
              b0(i,j) = bl(i,j) + br(i,j)
                   x0 = abs(b0(i,j))
                   xt = abs(bl(i,j)-br(i,j))
              smt5(i,j) =    x0 < xt
              smt6(i,j) = 3.*x0 < xt
           enddo
        enddo
        do j=jfirst,jlast
           do i=is,ie+1
              xt1(j) = c(i,j,p)
              hi5(j) = smt5(i-1,j) .and. smt5(i,j)
              hi6(j) = smt6(i-1,j) .or.  smt6(i,j)
              hi5(j) = hi5(j) .or. hi6(j)
           enddo
!DEC$ VECTOR ALWAYS
           do i=is,ie+1
              xt1(j) = c(i,j,p)
                if ( xt1(j) > 0. ) then
                     fx1(j) = (1.-xt1(j))*(br(i-1,j) - xt1(j)*b0(i-1,j))
                     flux(i,j,p) = q(i-1,j,p)
                else
                     fx1(j) = (1.+xt1(j))*(bl(i,j) + xt1(j)*b0(i,j))
                     flux(i,j,p) = q(i,j,p)
                endif
                if ( hi5(j) ) flux(i,j,p) = flux(i,j,p) + fx1(j)
           enddo
        enddo

   else  ! mord=5,6
       if ( iord==5 ) then
          do j=jfirst,jlast
             do i=is-1,ie+1
                bl(i,j) = al(i,j  ) - q(i,j,p)
                br(i,j) = al(i+1,j) - q(i,j,p)
                b0(i,j) = bl(i,j) + br(i,j)
                smt5(i,j) = bl(i,j)*br(i,j) < 0.
             enddo
          enddo
       elseif ( iord==-5 ) then
          do j=jfirst,jlast
             do i=is-1,ie+1
                bl(i,j) = al(i,j  ) - q(i,j,p)
                br(i,j) = al(i+1,j) - q(i,j,p)
                b0(i,j) = bl(i,j) + br(i,j)
                xt1(j) = br(i,j) - bl(i,j)
                 a4(j) = -3.*b0(i,j)
                smt5(i,j) = bl(i,j)*br(i,j) < 0.
             enddo
             do i=is-1,ie+1
                if( abs(xt1(j)) < -a4(j) ) then
                  if( q(i,j,p)+0.25/a4(j)*xt1(j)**2+a4(j)*r12 < 0. ) then
                    if( .not. smt5(i,j) ) then
                        br(i,j) = 0.
                        bl(i,j) = 0.
                        b0(i,j) = 0.
                    elseif( xt1(j) > 0. ) then
                        br(i,j) = -2.*bl(i,j)
                        b0(i,j) =    -bl(i,j)
                    else
                        bl(i,j) = -2.*br(i,j)
                        b0(i,j) =    -br(i,j)
                    endif
                  endif
                endif
             enddo
          enddo
       else
          do j=jfirst,jlast
             do i=is-1,ie+1
                bl(i,j) = al(i,j  ) - q(i,j,p)
                br(i,j) = al(i+1,j) - q(i,j,p)
                b0(i,j) = bl(i,j) + br(i,j)
                smt5(i,j) = 3.*abs(b0(i,j)) < abs(bl(i,j)-br(i,j))
             enddo
          enddo
       endif

       do j=jfirst,jlast
!DEC$ VECTOR ALWAYS
          do i=is,ie+1
             if ( c(i,j,p) > 0. ) then
                  fx1(i) = (1.-c(i,j,p))*(br(i-1,j) - c(i,j,p)*b0(i-1,j))
                  flux(i,j,p) = q(i-1,j,p)
             else
                  fx1(i) = (1.+c(i,j,p))*(bl(i,j) + c(i,j,p)*b0(i,j))
                  flux(i,j,p) = q(i,j,p)
             endif
             if (smt5(i-1,j).or.smt5(i,j)) flux(i,j,p) = flux(i,j,p) + fx1(j)
          enddo
       enddo

   endif
   return

else
! Monotonic constraints:
! ord = 8: PPM with Lin's PPM fast monotone constraint
! ord > 8: PPM with Lin's modification of Huynh 2nd constraint

  do j=jfirst,jlast
     do i=is-2,ie+2
             xt = 0.25*(q(i+1,j,p) - q(i-1,j,p))
        dm(i,j) = sign(min(abs(xt), max(q(i-1,j,p), q(i,j,p), q(i+1,j,p)) - q(i,j,p),   &
                           q(i,j,p) - min(q(i-1,j,p), q(i,j,p), q(i+1,j,p))), xt)
     enddo
  enddo
  do j=jfirst,jlast
     do i=is1,ie1+1
        al(i,j) = 0.5*(q(i-1,j,p)+q(i,j,p)) + r3*(dm(i-1,j) - dm(i,j))
        !al(i,j) = p1*(q(i-1,j)+q(i,j)) + p2*(q(i-2,j)+q(i+1,j))
     enddo
  enddo

  if ( iord==8 ) then
       do j=jfirst,jlast
          do i=is1,ie1
             xt = 2.*dm(i,j)
             bl(i,j) = -sign(min(abs(xt), abs(al(i,j)-q(i,j,p))),   xt)
             br(i,j) =  sign(min(abs(xt), abs(al(i+1,j)-q(i,j,p))), xt)
             !bl(i,j) =  al(i,j)-q(i,j)
             !br(i,j) =  al(i+1,j)-q(i,j)
          enddo
       enddo
  elseif ( iord==10 ) then
       do j=jfirst,jlast
          do i=is1-2,ie1+1
             dq(i,j) = 2.*(q(i+1,j,p) - q(i,j,p))
          enddo
       enddo
       do j=jfirst,jlast
          do i=is1,ie1
             bl(i,j) = al(i,j  ) - q(i,j,p)
             br(i,j) = al(i+1,j) - q(i,j,p)
             if ( abs(dm(i-1,j))+abs(dm(i,j))+abs(dm(i+1,j)) < near_zero ) then
                  bl(i,j) = 0.
                  br(i,j) = 0.
             elseif( abs(3.*(bl(i,j)+br(i,j))) > abs(bl(i,j)-br(i,j)) ) then
                  pmp_2 = dq(i-1,j)
                  lac_2 = pmp_2 - 0.75*dq(i-2,j)
                  br(i,j) = min(max(0.,pmp_2,lac_2), max(br(i,j), min(0.,pmp_2,lac_2)))
                  pmp_1 = -dq(i,j)
                  lac_1 = pmp_1 + 0.75*dq(i+1,j)
                  bl(i,j) = min(max(0.,pmp_1,lac_1), max(bl(i,j), min(0.,pmp_1,lac_1)))
             endif
          enddo
       enddo
  elseif ( iord==11 ) then
       do j=jfirst,jlast
          do i=is1,ie1
             xt = ppm_fac*dm(i,j)
             bl(i,j) = -sign(min(abs(xt), abs(al(i,j)-q(i,j,p))),   xt)
             br(i,j) =  sign(min(abs(xt), abs(al(i+1,j)-q(i,j,p))), xt)
          enddo
       enddo
  elseif ( iord==7 .or. iord==12 ) then
       do j=jfirst,jlast
          do i=is1,ie1
             bl(i,j) = al(i,j  ) - q(i,j,p)
             br(i,j) = al(i+1,j) - q(i,j,p)
              xt1(j) = br(i,j) - bl(i,j)
               a4(j) = -3.*(br(i,j) + bl(i,j))
              hi5(j) = bl(i,j)*br(i,j) > 0.
              hi6(j) = abs(xt1(j)) < -a4(j)
          enddo
          do i=is1,ie1
             if( hi6(i) ) then
                 if( q(i,j,p)+0.25/a4(j)*xt1(j)**2+a4(j)*r12 < 0. ) then
                    if( hi5(j) ) then
                        br(i,j) = 0.
                        bl(i,j) = 0.
                    elseif( xt1(j) > 0. ) then
                        br(i,j) = -2.*bl(i,j)
                    else
                        bl(i,j) = -2.*br(i,j)
                    endif
                 endif
             endif
          enddo
       enddo
  else
       do j=jfirst,jlast
          do i=is1,ie1
             bl(i,j) = al(i,j  ) - q(i,j,p)
             br(i,j) = al(i+1,j) - q(i,j,p)
          enddo
       enddo
  endif
  if ( iord==9 .or. iord==13 ) then
! Positive definite constraint:
     do i=is1,ie1
        call pert_ppm(jlast-jfirst+1, q(jfirst,j,p), bl(jfirst,j), br(jfirst,j), 0)
     enddo
  endif
  !if (.not. bounded_domain .and. grid_type<3) then
  !  if( is==1 ) then
  !    do j=jfirst,jlast
  !       bl(i,0) = s14*dm(i,-1) + s11*(q(i,-1)-q(i,0))
!
  !        xt = 0.5*(((2.*dya(i,0)+dya(i,-1))*q(i,0)-dya(i,0)*q(i,-1))/(dya(i,-1)+dya(i,0))   &
  !          +      ((2.*dya(i,1)+dya(i,2))*q(i,1)-dya(i,1)*q(i,2))/(dya(i,1)+dya(i,2)))
!        if ( iord==8 .or. iord==10 ) then
  !          xt = max(xt, min(q(i,-1),q(i,0),q(i,1),q(i,2)))
  !          xt = min(xt, max(q(i,-1),q(i,0),q(i,1),q(i,2)))
! !       endif
  !       br(i,0) = xt - q(i,0)
  !       bl(i,1) = xt - q(i,1)

  !       xt = s15*q(i,1) + s11*q(i,2) - s14*dm(i,2)
  !       br(i,1) = xt - q(i,1)
  !       bl(i,2) = xt - q(i,2)
!
  !       br(i,2) = al(i,3) - q(i,2)
  !    enddo
  !    call pert_ppm(3*(jlast-jfirst+1), q(jfirst,0), bl(jfirst,0), br(jfirst,0), 1)
  !  endif
   ! if( (ie+1)==npy ) then
   !   do j=jfirst,jlast
   !      bl(i,npy-2) = al(i,npy-2) - q(i,npy-2)

   !      xt = s15*q(i,npy-1) + s11*q(i,npy-2) + s14*dm(i,npy-2)
   !      br(i,npy-2) = xt - q(i,npy-2)
   !      bl(i,npy-1) = xt - q(i,npy-1)

   !      xt = 0.5*(((2.*dya(i,npy-1)+dya(i,npy-2))*q(i,npy-1)-dya(i,npy-1)*q(i,npy-2))/(dya(i,npy-2)+dya(i,npy-1))  &
   !         +      ((2.*dya(i,npy)+dya(i,npy+1))*q(i,npy)-dya(i,npy)*q(i,npy+1))/(dya(i,npy)+dya(i,npy+1)))
!        if ( iord==8 .or. iord==10 ) then
   !         xt = max(xt, min(q(i,npy-2),q(i,npy-1),q(i,npy),q(i,npy+1)))
   !         xt = min(xt, max(q(i,npy-2),q(i,npy-1),q(i,npy),q(i,npy+1)))
!        endif
   !      br(i,npy-1) = xt - q(i,npy-1)
   !     bl(i,npy  ) = xt - q(i,npy)

   !      br(i,npy) = s11*(q(i,npy+1)-q(i,npy)) - s14*dm(i,npy+1)
   !  enddo
   !  call pert_ppm(3*(jlast-jfirst+1), q(jfirst,npy-2), bl(jfirst,npy-2), br(jfirst,npy-2), 1)
   ! endif
 !end if

endif

  if ( iord==7 ) then
      do j=jfirst,jlast
         do i=is-1,ie+1
              b0(i,j) = bl(i,j) + br(i,j)
            smt5(i,j) = bl(i,j) * br(i,j) < 0.
         enddo
      enddo
      do j=jfirst,jlast
         do i=is,ie+1
            if ( c(i,j,p) > 0. ) then
                 fx1(i) = (1.-c(i,j,p))*(br(i-1,j) - c(i,j,p)*b0(i-1,j))
                 flux(i,j,p) = q(i-1,j,p)
            else
                 fx1(i) = (1.+c(i,j,p))*(bl(i,j) + c(i,j,p)*b0(i,j))
                 flux(i,j,p) = q(i,j,p)
            endif
            if ( smt5(i-1,j).or.smt5(i,j) ) flux(i,j,p) = flux(i,j,p) + fx1(i)
         enddo
      enddo
  else
      do j=jfirst,jlast
         do i=is,ie+1
            if( c(i,j,p)>0. ) then
                flux(i,j,p) = q(i-1,j,p) + (1.-c(i,j,p))*(br(i-1,j)-c(i,j,p)*(bl(i-1,j)+br(i-1,j)))
            else
                flux(i,j,p) = q(i,j  ,p) + (1.+c(i,j,p))*(bl(i,j  )+c(i,j,p)*(bl(i,j)+br(i,j)))
            endif
         enddo
      enddo
  endif
enddo

 end subroutine xppm


 subroutine yppm(flux, q, c, jord, ifirst,ilast, isd,ied, js,je,jsd,jed, lim_fac)
 integer, INTENT(IN) :: ifirst,ilast    ! Compute domain
 integer, INTENT(IN) :: isd,ied, js,je,jsd,jed
 integer, INTENT(IN) :: jord
 real(R_GRID), INTENT(IN) :: q(ifirst:ilast,jsd:jed,1:nbfaces)
 real(R_GRID), intent(in) :: c(isd:ied,js:je+1,1:nbfaces )  ! Courant number
 real(R_GRID), INTENT(OUT):: flux(ifirst:ilast,js:je+1,1:nbfaces)   !  Flux
 real(R_GRID), intent(IN) :: lim_fac
! Local:
 real(R_GRID):: dm(ifirst:ilast,js-2:je+2)
 real(R_GRID):: al(ifirst:ilast,js-1:je+2)
 real(R_GRID), dimension(ifirst:ilast,js-1:je+1):: bl, br, b0
 real(R_GRID):: dq(ifirst:ilast,js-3:je+2)
 real(R_GRID),    dimension(ifirst:ilast):: fx0, fx1, xt1, a4
 logical, dimension(ifirst:ilast,js-1:je+1):: smt5, smt6
 logical, dimension(ifirst:ilast):: hi5, hi6
 real(R_GRID):: x0, xt, qtmp, pmp_1, lac_1, pmp_2, lac_2
 integer:: i, j, js1, je3, je1, mord, p

   !if ( .not.bounded_domain .and. grid_type < 3 ) then
! Cubed-sphere:
   !   js1 = max(3,js-1); je3 = min(npy-2,je+2)
   !                      je1 = min(npy-3,je+1)
   !else
! Bounded_domain grid OR Doubly periodic domain:
      js1 = js-1;        je3 = je+2
                         je1 = je+1
   !endif

do p = 1, nbfaces
 mord = abs(jord)
if (jord==0) then
   do i=ifirst,ilast
      do j=js1, je3
         al(i,j) = p1*(q(i,j-1,p)+q(i,j,p)) + p2*(q(i,j-2,p)+q(i,j+1,p))
      enddo
   enddo

   do i=ifirst,ilast
      do j=js1, je1
         bl(i,j) =  al(i,j  )-q(i,j,p)
         br(i,j) =  al(i,j+1)-q(i,j,p)
      enddo
   enddo


   do i=ifirst,ilast
      do j=js,je+1
         if( c(i,j,p)>0. ) then
            flux(i,j,p) = q(i,j-1,p) + (1.-c(i,j,p))*(br(i,j-1)-c(i,j,p)*(bl(i,j-1)+br(i,j-1)))
         else
            flux(i,j,p) = q(i,j  ,p) + (1.+c(i,j,p))*(bl(i,j  )+c(i,j,p)*(bl(i,j  )+br(i,j  )))
         endif
      enddo
   enddo


else if ( jord < 7 ) then

   do j=js1, je3
      do i=ifirst,ilast
         al(i,j) = p1*(q(i,j-1,p)+q(i,j,p)) + p2*(q(i,j-2,p)+q(i,j+1,p))
      enddo
   enddo

   !if ( .not. bounded_domain .and. grid_type<3 ) then
   !   if( js==1 ) then
   !     do i=ifirst,ilast
   !        al(i,0) = c1*q(i,-2) + c2*q(i,-1) + c3*q(i,0)
   !        al(i,1) = 0.5*(((2.*dya(i,0)+dya(i,-1))*q(i,0)-dya(i,0)*q(i,-1))/(dya(i,-1)+dya(i,0))   &
   !                +      ((2.*dya(i,1)+dya(i,2))*q(i,1)-dya(i,1)*q(i,2))/(dya(i,1)+dya(i,2)))
   !        al(i,2) = c3*q(i,1) + c2*q(i,2) + c1*q(i,3)
   !     enddo
   !   endif
   !   if( (je+1)==npy ) then
   !     do i=ifirst,ilast
   !      al(i,npy-1) = c1*q(i,npy-3) + c2*q(i,npy-2) + c3*q(i,npy-1)
   !      al(i,npy) = 0.5*(((2.*dya(i,npy-1)+dya(i,npy-2))*q(i,npy-1)-dya(i,npy-1)*q(i,npy-2))/(dya(i,npy-2)+dya(i,npy-1))  &
   !                +      ((2.*dya(i,npy)+dya(i,npy+1))*q(i,npy)-dya(i,npy)*q(i,npy+1))/(dya(i,npy)+dya(i,npy+1)))
   !      al(i,npy+1) = c3*q(i,npy) + c2*q(i,npy+1) + c1*q(i,npy+2)
   !     enddo
   !   endif
   !endif

   if ( jord<0 ) then
      do j=js-1, je+2
         do i=ifirst,ilast
            al(i,j) = max(0., al(i,j))
         enddo
      enddo
   endif

   if ( mord==1 ) then
       do j=js-1,je+1
          do i=ifirst,ilast
             bl(i,j) = al(i,j  ) - q(i,j,p)
             br(i,j) = al(i,j+1) - q(i,j,p)
             b0(i,j) = bl(i,j) + br(i,j)
             smt5(i,j) = abs(lim_fac*b0(i,j)) < abs(bl(i,j)-br(i,j))
          enddo
       enddo
       do j=js,je+1
!DEC$ VECTOR ALWAYS
          do i=ifirst,ilast
             if ( c(i,j,p) > 0. ) then
                  fx1(i) = (1.-c(i,j,p))*(br(i,j-1) - c(i,j,p)*b0(i,j-1))
                  flux(i,j,p) = q(i,j-1,p)
             else
                  fx1(i) = (1.+c(i,j,p))*(bl(i,j) + c(i,j,p)*b0(i,j))
                  flux(i,j,p) = q(i,j,p)
             endif
             if (smt5(i,j-1).or.smt5(i,j)) flux(i,j,p) = flux(i,j,p) + fx1(i)
          enddo
       enddo

   elseif ( mord==2 ) then   ! Perfectly linear scheme
! Diffusivity: ord2 < ord5 < ord3 < ord4 < ord6  < ord7

      do j=js,je+1
!DEC$ VECTOR ALWAYS
         do i=ifirst,ilast
            xt = c(i,j,p)
            if ( xt > 0. ) then
                 qtmp = q(i,j-1,p)
                 flux(i,j,p) = qtmp + (1.-xt)*(al(i,j)-qtmp-xt*(al(i,j-1)+al(i,j)-(qtmp+qtmp)))
            else
                 qtmp = q(i,j,p)
                 flux(i,j,p) = qtmp + (1.+xt)*(al(i,j)-qtmp+xt*(al(i,j)+al(i,j+1)-(qtmp+qtmp)))
            endif
         enddo
      enddo

   elseif ( mord==3 ) then

        do j=js-1,je+1
           do i=ifirst,ilast
              bl(i,j) = al(i,j  ) - q(i,j,p)
              br(i,j) = al(i,j+1) - q(i,j,p)
              b0(i,j) = bl(i,j) + br(i,j)
                   x0 = abs(b0(i,j))
                   xt = abs(bl(i,j)-br(i,j))
              smt5(i,j) =    x0 < xt
              smt6(i,j) = 3.*x0 < xt
           enddo
        enddo
        do j=js,je+1
           do i=ifirst,ilast
              xt1(i) = c(i,j,p)
           enddo
           do i=ifirst,ilast
              if ( xt1(i) > 0. ) then
                   if( smt5(i,j-1) .or. smt6(i,j) ) then
                       flux(i,j,p) = q(i,j-1,p) + (1.-xt1(i))*(br(i,j-1) - xt1(i)*b0(i,j-1))
                   else
                       flux(i,j,p) = q(i,j-1,p)
                   endif
              else
                   if( smt6(i,j-1) .or. smt5(i,j) ) then
                       flux(i,j,p) = q(i,j,p) + (1.+xt1(i))*(bl(i,j) + xt1(i)*b0(i,j))
                   else
                       flux(i,j,p) = q(i,j,p)
                   endif
              endif
           enddo
        enddo

   elseif ( mord==4 ) then

        do j=js-1,je+1
           do i=ifirst,ilast
              bl(i,j) = al(i,j  ) - q(i,j,p)
              br(i,j) = al(i,j+1) - q(i,j,p)
              b0(i,j) = bl(i,j) + br(i,j)
                   x0 = abs(b0(i,j))
                   xt = abs(bl(i,j)-br(i,j))
              smt5(i,j) =    x0 < xt
              smt6(i,j) = 3.*x0 < xt
           enddo
        enddo
        do j=js,je+1
           do i=ifirst,ilast
              xt1(i) = c(i,j,p)
              hi5(i) = smt5(i,j-1) .and. smt5(i,j)
              hi6(i) = smt6(i,j-1) .or.  smt6(i,j)
              hi5(i) = hi5(i) .or. hi6(i)
           enddo
!DEC$ VECTOR ALWAYS
           do i=ifirst,ilast
                if ( xt1(i) > 0. ) then
                     fx1(i) = (1.-xt1(i))*(br(i,j-1) - xt1(i)*b0(i,j-1))
                     flux(i,j,p) = q(i,j-1,p)
                else
                     fx1(i) = (1.+xt1(i))*(bl(i,j) + xt1(i)*b0(i,j))
                     flux(i,j,p) = q(i,j,p)
                endif
                if ( hi5(i) ) flux(i,j,p) = flux(i,j,p) + fx1(i)
           enddo
        enddo

   else  ! mord=5,6
       if ( jord==5 ) then
          do j=js-1,je+1
             do i=ifirst,ilast
                bl(i,j) = al(i,j  ) - q(i,j,p)
                br(i,j) = al(i,j+1) - q(i,j,p)
                b0(i,j) = bl(i,j) + br(i,j)
                smt5(i,j) = bl(i,j)*br(i,j) < 0.
             enddo
          enddo
       elseif ( jord==-5 ) then
          do j=js-1,je+1
             do i=ifirst,ilast
                bl(i,j) = al(i,j  ) - q(i,j,p)
                br(i,j) = al(i,j+1) - q(i,j,p)
                b0(i,j) = bl(i,j) + br(i,j)
                xt1(i) = br(i,j) - bl(i,j)
                 a4(i) = -3.*b0(i,j)
                smt5(i,j) = bl(i,j)*br(i,j) < 0.
             enddo
             do i=ifirst,ilast
                if( abs(xt1(i)) < -a4(i) ) then
                  if( q(i,j,p)+0.25/a4(i)*xt1(i)**2+a4(i)*r12 < 0. ) then
                    if( .not. smt5(i,j) ) then
                        br(i,j) = 0.
                        bl(i,j) = 0.
                        b0(i,j) = 0.
                    elseif( xt1(i) > 0. ) then
                        br(i,j) = -2.*bl(i,j)
                        b0(i,j) =    -bl(i,j)
                    else
                        bl(i,j) = -2.*br(i,j)
                        b0(i,j) =    -br(i,j)
                    endif
                  endif
                endif
             enddo
          enddo
       else
          do j=js-1,je+1
             do i=ifirst,ilast
                bl(i,j) = al(i,j  ) - q(i,j,p)
                br(i,j) = al(i,j+1) - q(i,j,p)
                b0(i,j) = bl(i,j) + br(i,j)
                smt5(i,j) = 3.*abs(b0(i,j)) < abs(bl(i,j)-br(i,j))
             enddo
          enddo
       endif

       do j=js,je+1
!DEC$ VECTOR ALWAYS
          do i=ifirst,ilast
             if ( c(i,j,p) > 0. ) then
                  fx1(i) = (1.-c(i,j,p))*(br(i,j-1) - c(i,j,p)*b0(i,j-1))
                  flux(i,j,p) = q(i,j-1,p)
             else
                  fx1(i) = (1.+c(i,j,p))*(bl(i,j) + c(i,j,p)*b0(i,j))
                  flux(i,j,p) = q(i,j,p)
             endif
             if (smt5(i,j-1).or.smt5(i,j)) flux(i,j,p) = flux(i,j,p) + fx1(i)
          enddo
       enddo

   endif
   return

else
! Monotonic constraints:
! ord = 8: PPM with Lin's PPM fast monotone constraint
! ord > 8: PPM with Lin's modification of Huynh 2nd constraint

  do j=js-2,je+2
     do i=ifirst,ilast
             xt = 0.25*(q(i,j+1,p) - q(i,j-1,p))
        dm(i,j) = sign(min(abs(xt), max(q(i,j-1,p), q(i,j,p), q(i,j+1,p)) - q(i,j,p),   &
                           q(i,j,p) - min(q(i,j-1,p), q(i,j,p), q(i,j+1,p))), xt)
     enddo
  enddo
  do j=js1,je1+1
     do i=ifirst,ilast
        al(i,j) = 0.5*(q(i,j-1,p)+q(i,j,p)) + r3*(dm(i,j-1) - dm(i,j))
        !al(i,j) = p1*(q(i,j-1)+q(i,j)) + p2*(q(i,j-2)+q(i,j+1))
     enddo
  enddo

  if ( jord==8 ) then
       do j=js1,je1
          do i=ifirst,ilast
             xt = 2.*dm(i,j)
             bl(i,j) = -sign(min(abs(xt), abs(al(i,j)-q(i,j,p))),   xt)
             br(i,j) =  sign(min(abs(xt), abs(al(i,j+1)-q(i,j,p))), xt)
             !bl(i,j) =  al(i,j)-q(i,j)
             !br(i,j) =  al(i,j+1)-q(i,j)
          enddo
       enddo
  elseif ( jord==10 ) then
       do j=js1-2,je1+1
          do i=ifirst,ilast
             dq(i,j) = 2.*(q(i,j+1,p) - q(i,j,p))
          enddo
       enddo
       do j=js1,je1
          do i=ifirst,ilast
             bl(i,j) = al(i,j  ) - q(i,j,p)
             br(i,j) = al(i,j+1) - q(i,j,p)
             if ( abs(dm(i,j-1))+abs(dm(i,j))+abs(dm(i,j+1)) < near_zero ) then
                  bl(i,j) = 0.
                  br(i,j) = 0.
             elseif( abs(3.*(bl(i,j)+br(i,j))) > abs(bl(i,j)-br(i,j)) ) then
                  pmp_2 = dq(i,j-1)
                  lac_2 = pmp_2 - 0.75*dq(i,j-2)
                  br(i,j) = min(max(0.,pmp_2,lac_2), max(br(i,j), min(0.,pmp_2,lac_2)))
                  pmp_1 = -dq(i,j)
                  lac_1 = pmp_1 + 0.75*dq(i,j+1)
                  bl(i,j) = min(max(0.,pmp_1,lac_1), max(bl(i,j), min(0.,pmp_1,lac_1)))
             endif
          enddo
       enddo
  elseif ( jord==11 ) then
       do j=js1,je1
          do i=ifirst,ilast
             xt = ppm_fac*dm(i,j)
             bl(i,j) = -sign(min(abs(xt), abs(al(i,j)-q(i,j,p))),   xt)
             br(i,j) =  sign(min(abs(xt), abs(al(i,j+1)-q(i,j,p))), xt)
          enddo
       enddo
  elseif ( jord==7 .or. jord==12 ) then
       do j=js1,je1
          do i=ifirst,ilast
             bl(i,j) = al(i,j  ) - q(i,j,p)
             br(i,j) = al(i,j+1) - q(i,j,p)
              xt1(i) = br(i,j) - bl(i,j)
               a4(i) = -3.*(br(i,j) + bl(i,j))
              hi5(i) = bl(i,j)*br(i,j) > 0.
              hi6(i) = abs(xt1(i)) < -a4(i)
          enddo
          do i=ifirst,ilast
             if( hi6(i) ) then
                 if( q(i,j,p)+0.25/a4(i)*xt1(i)**2+a4(i)*r12 < 0. ) then
                    if( hi5(i) ) then
                        br(i,j) = 0.
                        bl(i,j) = 0.
                    elseif( xt1(i) > 0. ) then
                        br(i,j) = -2.*bl(i,j)
                    else
                        bl(i,j) = -2.*br(i,j)
                    endif
                 endif
             endif
          enddo
       enddo
  else
       do j=js1,je1
          do i=ifirst,ilast
             bl(i,j) = al(i,j  ) - q(i,j,p)
             br(i,j) = al(i,j+1) - q(i,j,p)
          enddo
       enddo
  endif
  if ( jord==9 .or. jord==13 ) then
! Positive definite constraint:
     do j=js1,je1
        call pert_ppm(ilast-ifirst+1, q(ifirst,j,p), bl(ifirst,j), br(ifirst,j), 0)
     enddo
  endif
  !if (.not. bounded_domain .and. grid_type<3) then
  !  if( js==1 ) then
  !    do i=ifirst,ilast
  !       bl(i,0) = s14*dm(i,-1) + s11*(q(i,-1)-q(i,0))
!
  !        xt = 0.5*(((2.*dya(i,0)+dya(i,-1))*q(i,0)-dya(i,0)*q(i,-1))/(dya(i,-1)+dya(i,0))   &
  !          +      ((2.*dya(i,1)+dya(i,2))*q(i,1)-dya(i,1)*q(i,2))/(dya(i,1)+dya(i,2)))
!        if ( jord==8 .or. jord==10 ) then
  !          xt = max(xt, min(q(i,-1),q(i,0),q(i,1),q(i,2)))
  !          xt = min(xt, max(q(i,-1),q(i,0),q(i,1),q(i,2)))
! !       endif
  !       br(i,0) = xt - q(i,0)
  !       bl(i,1) = xt - q(i,1)

  !       xt = s15*q(i,1) + s11*q(i,2) - s14*dm(i,2)
  !       br(i,1) = xt - q(i,1)
  !       bl(i,2) = xt - q(i,2)
!
  !       br(i,2) = al(i,3) - q(i,2)
  !    enddo
  !    call pert_ppm(3*(ilast-ifirst+1), q(ifirst,0), bl(ifirst,0), br(ifirst,0), 1)
  !  endif
   ! if( (je+1)==npy ) then
   !   do i=ifirst,ilast
   !      bl(i,npy-2) = al(i,npy-2) - q(i,npy-2)

   !      xt = s15*q(i,npy-1) + s11*q(i,npy-2) + s14*dm(i,npy-2)
   !      br(i,npy-2) = xt - q(i,npy-2)
   !      bl(i,npy-1) = xt - q(i,npy-1)

   !      xt = 0.5*(((2.*dya(i,npy-1)+dya(i,npy-2))*q(i,npy-1)-dya(i,npy-1)*q(i,npy-2))/(dya(i,npy-2)+dya(i,npy-1))  &
   !         +      ((2.*dya(i,npy)+dya(i,npy+1))*q(i,npy)-dya(i,npy)*q(i,npy+1))/(dya(i,npy)+dya(i,npy+1)))
!        if ( jord==8 .or. jord==10 ) then
   !         xt = max(xt, min(q(i,npy-2),q(i,npy-1),q(i,npy),q(i,npy+1)))
   !         xt = min(xt, max(q(i,npy-2),q(i,npy-1),q(i,npy),q(i,npy+1)))
!        endif
   !      br(i,npy-1) = xt - q(i,npy-1)
   !     bl(i,npy  ) = xt - q(i,npy)

   !      br(i,npy) = s11*(q(i,npy+1)-q(i,npy)) - s14*dm(i,npy+1)
   !  enddo
   !  call pert_ppm(3*(ilast-ifirst+1), q(ifirst,npy-2), bl(ifirst,npy-2), br(ifirst,npy-2), 1)
   ! endif
 !end if

endif

  if ( jord==7 ) then
      do j=js-1,je+1
         do i=ifirst,ilast
              b0(i,j) = bl(i,j) + br(i,j)
            smt5(i,j) = bl(i,j) * br(i,j) < 0.
         enddo
      enddo
      do j=js,je+1
         do i=ifirst,ilast
            if ( c(i,j,p) > 0. ) then
                 fx1(i) = (1.-c(i,j,p))*(br(i,j-1) - c(i,j,p)*b0(i,j-1))
                 flux(i,j,p) = q(i,j-1,p)
            else
                 fx1(i) = (1.+c(i,j,p))*(bl(i,j) + c(i,j,p)*b0(i,j))
                 flux(i,j,p) = q(i,j,p)
            endif
            if ( smt5(i,j-1).or.smt5(i,j) ) flux(i,j,p) = flux(i,j,p) + fx1(i)
         enddo
      enddo
  else
      do j=js,je+1
         do i=ifirst,ilast
            if( c(i,j,p)>0. ) then
                flux(i,j,p) = q(i,j-1,p) + (1.-c(i,j,p))*(br(i,j-1)-c(i,j,p)*(bl(i,j-1)+br(i,j-1)))
            else
                flux(i,j,p) = q(i,j  ,p) + (1.+c(i,j,p))*(bl(i,j  )+c(i,j,p)*(bl(i,j)+br(i,j)))
            endif
         enddo
      enddo
  endif
enddo
 end subroutine yppm


 subroutine pert_ppm(im, a0, al, ar, iv)
 integer, intent(in):: im
 integer, intent(in):: iv
 real(R_GRID), intent(in)   :: a0(im)
 real(R_GRID), intent(inout):: al(im), ar(im)
! Local:
 real(R_GRID) a4, da1, da2, a6da, fmin
 integer i

!-----------------------------------
! Optimized PPM in perturbation form:
!-----------------------------------

 if ( iv==0 ) then
! Positive definite constraint
    do i=1,im
     if ( a0(i) <= 0. ) then
          al(i) = 0.
          ar(i) = 0.
     else
        a4 = -3.*(ar(i) + al(i))
       da1 =      ar(i) - al(i)
      if( abs(da1) < -a4 ) then
         fmin = a0(i) + 0.25/a4*da1**2 + a4*r12
         if( fmin < 0. ) then
             if( ar(i)>0. .and. al(i)>0. ) then
                 ar(i) = 0.
                 al(i) = 0.
             elseif( da1 > 0. ) then
                 ar(i) = -2.*al(i)
             else
                 al(i) = -2.*ar(i)
             endif
         endif
      endif
     endif
    enddo
 else
! Standard PPM constraint
    do i=1,im
       if ( al(i)*ar(i) < 0. ) then
            da1 = al(i) - ar(i)
            da2 = da1**2
            a6da = 3.*(al(i)+ar(i))*da1
! abs(a6da) > da2 --> 3.*abs(al+ar) > abs(al-ar)
            if( a6da < -da2 ) then
                ar(i) = -2.*al(i)
            elseif( a6da > da2 ) then
                al(i) = -2.*ar(i)
            endif
       else
! effect of dm=0 included here
            al(i) = 0.
            ar(i) = 0.
       endif
  enddo
 endif

 end subroutine pert_ppm
end module tp_core 
