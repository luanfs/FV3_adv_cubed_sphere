module atmosphere
!========================================================================
!========================================================================
use fv_control, only: init_model
use fv_duogrid, only: ext_scalar_agrid
use fv_arrays , only: fv_atmos_type, datadir, griddir, pardir, R_GRID, day2sec, sec2day, nbfaces
use test_cases, only: init_case
use dyn_core  , only: dy_core

implicit none

contains
!--------------------------------------------------------------
! main atmosphere routine
!--------------------------------------------------------------
subroutine atmosphere_main(atm)
   type(fv_atmos_type), intent(inout) :: atm
   integer :: tstep
   logical :: first_step=.false.

   ! init everything
   call atmosphere_init(atm)

   ! loop over time
   do tstep = 1, atm%total_tsteps
      ! compute one time step
      call atmosphere_timestep(atm)

      ! compute diagnostics
      call atmosphere_diag(atm, first_step)

      ! output data
      call atmosphere_output(atm, tstep)
   enddo

   ! compute error
   call atmosphere_end(atm)
end subroutine atmosphere_main

!--------------------------------------------------------------
! init the model
!--------------------------------------------------------------
subroutine atmosphere_init(atm)
   type(fv_atmos_type), intent(inout) :: atm
   logical :: first_step=.true.
   integer :: p
   integer :: is, js, ie, je
   integer :: isd, jsd, ied, jed
   integer :: i, j, g
   real(R_GRID), allocatable :: q(:,:,:)

   call atmosphere_input(atm)

   ! initialize all variables
   call init_model(atm)

   is = atm%bd%is
   ie = atm%bd%ie
   js = atm%bd%js
   je = atm%bd%je

   isd = atm%bd%isd
   ied = atm%bd%ied
   jsd = atm%bd%jsd
   jed = atm%bd%jed

   ! time vars
   atm%time          = 0.d0
   atm%time_centered = atm%time + atm%dto2

   ! compute ICs
   call init_case(atm)

   ! cfl
   atm%cfl_x = 0.d0
   atm%cfl_y = 0.d0
   do p = 1, nbfaces
      atm%cfl_x = max(atm%cfl_x, maxval(abs(atm%uc0(is:ie+1,js:je  ,p))/atm%gridstruct%dx_u(is:ie+1,js:je  ))*atm%dt)
      atm%cfl_y = max(atm%cfl_y, maxval(abs(atm%vc0(is:ie  ,js:je+1,p))/atm%gridstruct%dy_v(is:ie  ,js:je+1))*atm%dt)
   enddo
   atm%cfl = max(atm%cfl_x, atm%cfl_y)

   ! update qa and uc
   atm%qa(is:ie,js:je,:)= atm%qa0(is:ie,js:je,:)
   atm%uc_old = atm%uc0
   atm%vc_old = atm%vc0
 
   !call ext_scalar_agrid(atm%qa, atm%bd, atm%L)
   !print*, maxval(abs(atm%qa - atm%qa0))/maxval(abs(atm%qa0))


   print*,"test case  :", atm%test_case
   print*,"npx        :", atm%npx
   print*,"dt         :", atm%dt
   print*,"hord       :", atm%hord
   print*,"dp         :", atm%dp
   print*,"inner_adv  :", atm%inner_adv
   print*,"nplots     :", atm%nplots
   print*,"adjusted dt:", atm%dt
   print*,"cfl        :", atm%cfl
   print*,'------------------------------------------------------------------'

   ! compute initial diagnostics
   call atmosphere_diag(atm, first_step)

   ! plot IC
   call atmosphere_output(atm, 0)
end subroutine atmosphere_init

!--------------------------------------------------------------
! compute on timestep of atmosphere dynamics
!--------------------------------------------------------------
subroutine atmosphere_timestep(atm)
   type(fv_atmos_type), intent(inout) :: atm
   logical :: first_step=.false.

   ! solves dynamics
   call dy_core(atm%qa, atm%uc, atm%uc_old, atm%vc, atm%vc_old, atm%bd, atm%gridstruct, atm%time, atm%time_centered,&
                   atm%dt, atm%dto2, atm%test_case, atm%hord, atm%lim_fac, atm%dp, atm%inner_adv, atm%L)

   ! update times
   atm%time          = atm%time + atm%dt
   atm%time_centered = atm%time + atm%dto2

end subroutine atmosphere_timestep


!--------------------------------------------------------------
! atmosphere output
!--------------------------------------------------------------
subroutine atmosphere_output(atm, step)
   type(fv_atmos_type), intent(inout) :: atm
   integer, intent(in) :: step
   integer :: is, ie
   integer :: js, je
   integer :: i, j, iunit, p
   character (len=60):: nplot
   character (len=60):: panel
   character (len=60):: filename

   if(step==0 .or. step==atm%total_tsteps .or. mod(step,atm%plotstep)==0 )then
      write(nplot, '(i8)') atm%nplot
      filename = trim(datadir)//trim(atm%simulation_name)//"t"//trim(adjustl(nplot))//".txt"
      print*, 'saving ', filename

      ! write the q data in a text file
      call getunit(iunit)
      open(iunit, file=filename, status='replace')
      write(iunit,*) atm%time*sec2day
      write(iunit,*) atm%mass_qa_var
      write(iunit,*) atm%cfl
      close(iunit)

      is = atm%bd%is
      ie = atm%bd%ie
      js = atm%bd%js
      je = atm%bd%je
 
      do p = 1, nbfaces
         write(panel, '(i8)') p

         ! save scalar field in binary format
         call getunit(iunit)
         filename = trim(datadir)//trim(atm%simulation_name)//"t"//trim(adjustl(nplot))//"_face"//trim(adjustl(panel))//".dat"
         print*, 'saving ', filename

         !Write whole block to file (much faster)
         open(iunit, file=filename, status='replace', access='stream', form='unformatted')
         write(iunit) atm%qa(is:ie,js:je, p)
         close(iunit)
      enddo
      atm%nplot = atm%nplot + 1
      print*
      print*
   endif

   if(step<0) then
      do p = 1, nbfaces
         is = atm%bd%is
         ie = atm%bd%ie
         js = atm%bd%js
         je = atm%bd%je
         write(panel, '(i8)') p

         ! save scalar field in binary format
         call getunit(iunit)
         filename = trim(griddir)//trim(atm%grid_name)//"_face"//trim(adjustl(panel))//"_lon.dat"
         print*, 'saving ', filename

         !Write whole block to file (much faster)
         open(iunit, file=filename, status='replace', access='stream', form='unformatted')
         write(iunit) atm%gridstruct%bgrid(is:ie+1,js:je+1,p)%lon
         close(iunit)

         call getunit(iunit)
         filename = trim(griddir)//trim(atm%grid_name)//"_face"//trim(adjustl(panel))//"_lat.dat"
         print*, 'saving ', filename

         !Write whole block to file (much faster)
         open(iunit, file=filename, status='replace', access='stream', form='unformatted')
         write(iunit) atm%gridstruct%bgrid(is:ie+1,js:je+1,p)%lat
         close(iunit)
 
      enddo
   endif
 
end subroutine atmosphere_output

!--------------------------------------------------------------
! Compute diagnostics
!--------------------------------------------------------------
subroutine atmosphere_diag(atm, first_step)
   type(fv_atmos_type), intent(inout) :: atm
   logical, intent(in) :: first_step
   integer :: is, ie
   integer :: js, je
   integer :: p

   is = atm%bd%is
   ie = atm%bd%ie
   js = atm%bd%js
   je = atm%bd%je
 
  if(first_step) then
     atm%mass_qa0 = 0.d0
     do p = 1, nbfaces
        atm%mass_qa0 = atm%mass_qa0 + sum(atm%qa0(is:ie,js:je,p) * atm%gridstruct%area(is:ie,js:je))
     enddo
  else
     atm%mass_qa = 0.d0
     do p = 1, nbfaces
        atm%mass_qa = atm%mass_qa + sum(atm%qa(is:ie,js:je,p) * atm%gridstruct%area(is:ie,js:je))
     enddo
  endif

   if(.not. first_step) then
      atm%mass_qa_var = (atm%mass_qa0-atm%mass_qa)/atm%mass_qa0
   endif
end subroutine atmosphere_diag


!--------------------------------------------------------------
! end the model
!--------------------------------------------------------------
subroutine atmosphere_end(atm)
   type(fv_atmos_type), intent(inout) :: atm
   integer :: i, iunit 
   character (len=60):: filename
   integer :: is, ie
   integer :: js, je

   is = atm%bd%is
   ie = atm%bd%ie
   js = atm%bd%js
   je = atm%bd%je

   ! compute errors of qa
   atm%error_qa(is:ie,js:je,:) = abs(atm%qa(is:ie,js:je,:)-atm%qa0(is:ie,js:je,:))
   atm%linf_error_qa   = maxval(atm%error_qa(is:ie,js:je,:))/maxval(abs(atm%qa0(is:ie,js:je,:)))
   atm%l1_error_qa     = sum   (atm%error_qa(is:ie,js:je,:))/sum   (abs(atm%qa0(is:ie,js:je,:)))
   atm%l2_error_qa     = norm2 (atm%error_qa(is:ie,js:je,:))/norm2 (abs(atm%qa0(is:ie,js:je,:)))

   filename = trim(datadir)//trim(atm%simulation_name)//"errors.txt"
   print*, 'saving ', filename
   call getunit(iunit)
   open(iunit, file=filename, status='replace')
   write(iunit,*) atm%linf_error_qa
   write(iunit,*) atm%l1_error_qa
   write(iunit,*) atm%l2_error_qa
   close(iunit)
   print '(a33, 3e16.8)','(linf, l1, l2) error norms:', &
   atm%linf_error_qa, atm%l1_error_qa, atm%l2_error_qa
   print '(a33, 1e16.8)','mass variation            :', &
   atm%mass_qa_var
   print*,'------------------------------------------------------------------'
end subroutine atmosphere_end

subroutine atmosphere_input(atm)
    !---------------------------------------------------
    ! read parameters from file par/input.par
    !--------------------------------------------------
    type(fv_atmos_type), intent(inout):: atm
    character (len=60):: filename
    character (len=300):: buffer
    integer :: fileunit
    integer:: i
    integer:: n

    !Standard advection parameters file
    filename=trim(pardir)//"input.par"

    print*,"Input parameters: ", trim(filename)
    print*

    fileunit = 7
    !A parameters file must exist 
    open(fileunit,file=filename,status='old')

    read(fileunit,*)  buffer
    read(fileunit,*)  buffer
    read(fileunit,*)  atm%test_case
    read(fileunit,*)  buffer
    read(fileunit,*)  atm%npx
    read(fileunit,*)  buffer
    read(fileunit,*)  atm%dt
    read(fileunit,*)  buffer
    read(fileunit,*)  atm%hord
    read(fileunit,*)  buffer
    read(fileunit,*)  atm%dp
    read(fileunit,*)  buffer
    read(fileunit,*)  atm%inner_adv
    read(fileunit,*)  buffer
    read(fileunit,*)  atm%nplots
    close(fileunit)

    ! Time vars
    atm%Tf   = 12.d0 * day2sec
    atm%dto2 = atm%dt*0.5d0
    atm%total_tsteps  = int(atm%Tf/atm%dt)
    atm%plotstep = atm%total_tsteps/atm%nplots
    ! Readjust time step
    atm%dt  = atm%Tf/atm%total_tsteps
    atm%npy = atm%npx
    return
end subroutine atmosphere_input

  !===============================================================================================
  !    MISCELLANEOUS
  !===============================================================================================

  subroutine getunit ( iunit )
    !----------------------------------------------------------
    ! GETUNIT returns a free FORTRAN unit number.
    !
    !    A "free" FORTRAN unit number is an integer between 1 and 99 which
    !    is not currently associated with an I/O device.  A free FORTRAN unit
    !    number is needed in order to open a file with the OPEN command.
    !
    !    If IUNIT = 0, then no free FORTRAN unit could be found, although
    !    all 99 units were checked (except for units 5, 6 and 9, which
    !    are commonly reserved for console I/O).
    !
    !    Otherwise, IUNIT is an integer between 1 and 99, representing a
    !    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
    !    are special, and will never return those values.
    !
    !    John Burkardt
    !    18 September 2005
    !----------------------------------------------------------------------------
    integer :: i
    integer :: ios
    integer :: iunit
    logical:: lopen

    iunit = 0
    do i = 11, 99
       if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then
          inquire ( unit = i, opened = lopen, iostat = ios )
          if ( ios == 0 ) then
             if ( .not. lopen ) then
                iunit = i
                return
             end if
          end if
       end if
    end do

    return
  end subroutine getunit
end module atmosphere
