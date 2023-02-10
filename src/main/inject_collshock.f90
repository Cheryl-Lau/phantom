!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module inject
!
! Routine for injecting a 1D colliding shock at centre of box
!
! :References: None
!
! :Owner: Cheryl Lau
!
! :Runtime parameters:
!   - collaxis : *direction to divide box and inject shock*
!   - vel_sk   : *velocity magnitude of shock*
!   - time_sk  : *time of shock injection*
!
! :Dependencies: eos, infile_utils, io, part, partinject, physcon
!
 implicit none
 character(len=*), parameter, public :: inject_type = 'shock'

 public :: init_inject,inject_particles,write_options_inject,read_options_inject

 ! Runtime parameters for shock injection to read from input file
 ! in code units of setup_unifdis_collshock
 integer, public :: collaxis = 1     ! direction of shock
 real, public    :: vel_sk   = 76.2  ! code units
 real, public    :: time_sk  = 2E-5

 private
 real    :: vel_sh_vec(3)
 logical :: sk_flag

contains
!-----------------------------------------------------------------------
!+
!  Initialize global variables or arrays needed for injection routine
!+
!-----------------------------------------------------------------------
subroutine init_inject(ierr)
 integer, intent(out) :: ierr
 !
 ! return without error
 !
 ierr = 0
 !
 ! Set velocity vector of shock
 !
 vel_sh_vec(:) = (/ 0.,0.,0. /)
 vel_sh_vec(collaxis) = vel_sk
 !
 ! Flag when shock is injected
 !
 sk_flag = .false.

end subroutine init_inject

!-----------------------------------------------------------------------
!+
!  Main routine handling shock injection
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                            npart,npartoftype,dtinject)
 use io,         only:id,master
 use part,       only:rhoh,massoftype,igas
 use partinject, only:add_or_update_particle
 use timestep,   only:dtmax
 use boundary,   only:xmin,ymin,zmin,xmax,ymax,zmax
 real,    intent(in)    :: time,dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 integer :: i
 real    :: xi
 logical :: inject_sk = .true.

 !
 ! Inject at specified time
 !
 inject_sk = abs(time_sk - time) < dtmax
 !
 ! Inject colliding shocks by changing particle velocities
 !
 if (inject_sk .and. .not.sk_flag) then
    print*,'Injecting shock'

    do i = 1,npart
       xi = xyzh(collaxis,i)
       if (xi <= 0.) then
          vxyzu(:,i) = vxyzu(:,i) + vel_sh_vec(:)
          xyzh(4,i) = 1.
       elseif (xi > 0.) then
          vxyzu(:,i) = vxyzu(:,i) - vel_sh_vec(:)
       endif
    enddo
    sk_flag = .true.
 endif
 !
 ! timestep constraint
 !
 dtinject = huge(dtinject)

end subroutine inject_particles

!-----------------------------------------------------------------------
!+
!  Writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use infile_utils, only: write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options for injecting shock'
 call write_inopt(time_sk,'time_sk','time of injecting shock',iunit)
 call write_inopt(vel_sk,'vel_sk','velocity of shock',iunit)
 call write_inopt(collaxis,'collaxis','direction of shock',iunit)

end subroutine write_options_inject

!-----------------------------------------------------------------------
!+
!  Reads input options from the input file.
!+
!-----------------------------------------------------------------------
subroutine read_options_inject(name,valstring,imatch,igotall,ierr)
 use io,  only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical, intent(out) :: imatch,igotall
 integer,intent(out) :: ierr
 integer, save :: ngot = 0
 integer       :: noptions
 character(len=30), parameter :: label = 'read_options_inject'

 imatch  = .true.
 igotall = .false.

 select case(trim(name))
 case('time_sk')
    read(valstring,*,iostat=ierr) time_sk
    ngot = ngot + 1
    if (time_sk < 0.) call fatal(label,'invalid setting for time_sk (<0)')
case('vel_sk')
    read(valstring,*,iostat=ierr) vel_sk
    ngot = ngot + 1
case('collaxis')
    read(valstring,*,iostat=ierr) collaxis
    ngot = ngot + 1
 case default
    imatch = .false.
 end select

 noptions = 3
 igotall  = (ngot >= noptions)

end subroutine read_options_inject

end module inject
