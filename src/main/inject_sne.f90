!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module inject
!
! Routine for injecting one supernova at user-specified time and location
! Modified based on phantom test case of Balsara & Kim (2004)
!
! :References: Balsara & Kim (2004), ApJ 602, 1079
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - x_sn     : *x-coord of supernova*
!   - y_sn     : *y-coord of supernova*
!   - z_sn     : *z-coord of supernova*
!   - time_sn  : *time of supernova*
!   - r_sn     : *supernova blast radius*
!   - pr_sn    : *new pressure within blast radius*
!
! :Dependencies: eos, infile_utils, io, part, partinject, physcon
!
 implicit none
 character(len=*), parameter, public :: inject_type = 'supernovae'

 public :: init_inject,inject_particles,write_options_inject,read_options_inject

 ! Runtime parameters for supernovae injection to read from input file
 ! in code units of setup_sphereinbox
 real, public :: x_sn    = 4.5
 real, public :: y_sn    = 0.
 real, public :: z_sn    = 0.
 real, public :: time_sn = 0.005
 real, public :: r_sn    = 1.
 real, public :: pr_sn   = 1.0d9

 private
 logical  :: sn_flag

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

 sn_flag = .true. ! flag when one sn has occurred

end subroutine init_inject

!-----------------------------------------------------------------------
!+
!  Main routine handling supernovae injection
!  Note that we actually only inject thermal energy, not kinetic energy
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,&
                            npart,npartoftype,dtinject)
 use io,      only:id,master
 use eos,     only:gamma
 use part,    only:rhoh,massoftype,iphase,igas,iunknown
 use partinject, only: updated_particle
 real,    intent(in)    :: time,dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 integer :: i,ipart
 real    :: xyz_sn(3),dx(3),uval,r2,rhoi
 logical :: inject_sn
 !
 ! Inject one sn at specified time
 !
 if (sn_flag) then
    inject_sn = abs(time_sn - time) < 1.e-4
 endif
 !timestep constraint
 dtinject = huge(dtinject)
 !
 !--inject sn by changing internal energy of particles
 !
 if (inject_sn .and. sn_flag) then
    xyz_sn = (/ x_sn, y_sn, z_sn /)
    if (id==master) print*, 'Injecting supernova at ',xyz_sn(1:3)
    print*,' gamma = ',gamma
    ipart = 0
    do i = 1,npart
       dx = xyzh(1:3,i) - xyz_sn(1:3)
       r2 = dot_product(dx,dx)
       if (r2 < r_sn**2) then
          rhoi = rhoh(xyzh(4,i),massoftype(igas))
          uval = pr_sn / ((gamma - 1.)*rhoi)
          print*,'Old & New thermal energy',vxyzu(4,i),uval
          vxyzu(4,i) = uval
          iphase(i)  = iunknown ! flag this particle to update its timestep
          ipart      = ipart + 1
!          dtinject   = min(dtinject,0.01*xyzh(4,i)*sqrt(rhoi/(gamma*pr_sn)))
          updated_particle = .true.
       endif
    enddo
    print*,' Energy injected into ',ipart,' particles'
    sn_flag = .false.
    print*,'--------'
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

 write(iunit,"(/,a)") '# options for injecting supernova'
 call write_inopt(time_sn,'time_sn','time of supernova injection',iunit)
 call write_inopt(x_sn,'x_sn','x-location of supernova injection',iunit)
 call write_inopt(y_sn,'y_sn','y-location of supernova injection',iunit)
 call write_inopt(z_sn,'z_sn','z-location of supernova injection',iunit)
 call write_inopt(r_sn,'r_sn','blast radius of supernova',iunit)
 call write_inopt(pr_sn,'pr_sn','new pressure within supernova blast radius',iunit)

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
 case('time_sn')
    read(valstring,*,iostat=ierr) time_sn
    ngot = ngot + 1
    if (time_sn < 0.) call fatal(label,'invalid setting for time_sn (<0)')
 case('x_sn')
    read(valstring,*,iostat=ierr) x_sn
    ngot = ngot + 1
 case('y_sn')
    read(valstring,*,iostat=ierr) y_sn
    ngot = ngot + 1
 case('z_sn')
    read(valstring,*,iostat=ierr) z_sn
    ngot = ngot + 1
 case('r_sn')
    read(valstring,*,iostat=ierr) r_sn
    ngot = ngot + 1
    if (r_sn < 0.) call fatal(label,'invalid setting for r_sn (<0)')
 case('pr_sn')
    read(valstring,*,iostat=ierr) pr_sn
    ngot = ngot + 1
    if (pr_sn < 0.) call fatal(label,'invalid setting for pr_sn (<0)')
 case default
    imatch = .false.
 end select

 noptions = 6
 igotall  = (ngot >= noptions)

end subroutine read_options_inject

end module inject
