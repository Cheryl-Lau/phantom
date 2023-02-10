!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module inject
!
! Routine for injecting a simple shock (for testing)
!
! :References: None
!
! :Owner: Cheryl Lau
!
! :Runtime parameters:
!   - x/y/zmin/max_sk  : *boundary of region to apply 1-D shock*
!
! :Dependencies: eos, infile_utils, io, part, partinject, physcon
!
 implicit none
 character(len=*), parameter, public :: inject_type = '1d-shock'

 public :: init_inject,inject_particles,write_options_inject,read_options_inject

 ! in code units of setup_unifdis_shock
 real, public :: xmin_sk = -10.
 real, public :: xmax_sk = 10.
 real, public :: ymin_sk = -10.
 real, public :: ymax_sk = 10.
 real, public :: zmin_sk = -10.
 real, public :: zmax_sk = 10.
 real, public :: time_sk = 1E-7

 private
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

 sk_flag = .false. ! flag when occurred

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
 real,    intent(in)    :: time,dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 integer :: i
 real    :: x,y,z,vxyz_sk(3),u_sk,old_v(4)
 logical :: inject_sk

 !
 ! Inject at specified time
 !
 inject_sk = abs(time_sk - time) < dtmax
 !
 ! timestep constraint
 !
 dtinject = huge(dtinject)

 !
 ! Inject shock by changing KE and TE of particles
 !
 if (inject_sk .and. .not.sk_flag) then
    print*,'Injecting shock'

    do i = 1,npart
       x = xyzh(1,i)
       y = xyzh(2,i)
       z = xyzh(3,i)
       if ((x >= xmin_sk) .and. (x <= xmax_sk) .and. &
           (y >= ymin_sk) .and. (y <= ymax_sk) .and. &
           (z >= zmin_sk) .and. (z <= zmax_sk)) then
!          print*,'updating particle',i
          old_v = vxyzu(:,i)
          vxyz_sk = (/ vxyzu(1,i), 250., vxyzu(3,i) /)
          call add_or_update_particle(igas,xyzh(1:3,i),vxyz_sk(1:3),xyzh(4,i),vxyzu(4,i),i,npart,&
                                     npartoftype,xyzh,vxyzu)
          print*,'updated particle: old new vxyzu',old_v,vxyzu(:,i)
       endif
    enddo
    sk_flag = .true.
 endif

end subroutine inject_particles

!-----------------------------------------------------------------------
!+
!  Writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use physcon,      only: au, solarm, years
 use infile_utils, only: write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options for injecting shock'
 call write_inopt(time_sk,'time_sk','time of injecting shock',iunit)
 call write_inopt(xmin_sk,'xmin_sk','boundary of shock region',iunit)
 call write_inopt(xmax_sk,'xmax_sk','boundary of shock region',iunit)
 call write_inopt(ymin_sk,'ymin_sk','boundary of shock region',iunit)
 call write_inopt(ymax_sk,'ymax_sk','boundary of shock region',iunit)
 call write_inopt(zmin_sk,'zmin_sk','boundary of shock region',iunit)
 call write_inopt(zmax_sk,'zmax_sk','boundary of shock region',iunit)

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
case('xmin_sk')
    read(valstring,*,iostat=ierr) xmin_sk
    ngot = ngot + 1
case('xmax_sk')
    read(valstring,*,iostat=ierr) xmax_sk
    ngot = ngot + 1
case('ymin_sk')
    read(valstring,*,iostat=ierr) ymin_sk
    ngot = ngot + 1
case('ymax_sk')
    read(valstring,*,iostat=ierr) ymax_sk
    ngot = ngot + 1
case('zmin_sk')
    read(valstring,*,iostat=ierr) zmin_sk
    ngot = ngot + 1
case('zmax_sk')
    read(valstring,*,iostat=ierr) zmax_sk
    ngot = ngot + 1
 case default
    imatch = .false.
 end select

 noptions = 7
 igotall  = (ngot >= noptions)

end subroutine read_options_inject

end module inject
