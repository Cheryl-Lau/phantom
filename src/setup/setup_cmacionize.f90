!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: setup
!
!  DESCRIPTION:
!   Setup routine for uniform distribution
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: boundary, dim, dust, io, mpiutils, part, physcon,
!    prompting, setup_params, unifdis
!+
!--------------------------------------------------------------------------
module setup
 implicit none
 public :: setpart

 integer, private :: npartx,ilattice
 real,    private :: polykset
 private

contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use setup_params, only:rhozero,npart_total
 use io,           only:master
 use unifdis,      only:set_unifdis
 use boundary,     only:xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound
 use mpiutils,     only:bcast_mpi
 use part,         only:dustfrac
 use prompting,    only:prompt
 use physcon,      only:pi,pc,solarm,mass_proton_cgs,kboltz

 use units,        only:umass,udist,utime,set_units
 use eos,          only : gmw

 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 real :: totmass,deltax, factor
 integer :: i,maxp,maxvxyzu, ierr
!
!--general parameters
!
 time = 0.
 hfact = 1.2
 gamma = 1.00011
 gmw = 1.

 !number of particles per dimension
 npartx = 64
 !density
 rhozero = 5.21e-21 / solarm * pc*pc*pc
 totmass = npartx*npartx*npartx/1000.

 !Set the units so that the total mass = 1, and the box goes from -0.5 to 0.5
 call set_units(dist=(totmass/rhozero)**(1./3.)*pc,mass=totmass*solarm,G=1.)

!
!--setup particles
!
 maxp = size(xyzh(1,:))
 maxvxyzu = size(vxyzu(:,1))

 ! Do I need these two?
 polyk = 0.
 polykset = 0.


 npart = npartx*npartx*npartx
 npart_total = npart

 ! Open glass distribution file
 open(unit=505,file="../glassCube_64.dat",iostat=ierr)
 if (ierr /= 0) then
    write(*,*) 'ERROR opening glassCube_64.dat'
    return
 endif

 ! Read in glass distribution
 do i=1,npart
   read(505,*) xyzh(1,i),xyzh(2,i),xyzh(3,i),xyzh(4,i)
 end do

 close(505)

 npartoftype(:) = 0
 npartoftype(1) = npart
 print*,' npart = ',npart,npart_total

 !massoftype = totmass/npart_total
 massoftype = 1.0/npart_total
 print*,' particle mass = ',massoftype(1)

 ! Set velocities to 0, and particle energies to correspond to 100 K
 factor = 1.0/(mass_proton_cgs/kboltz * (udist/utime)**2*gmw*(gamma-1))
 do i=1,npart
    vxyzu(1:3,i) = 0.
    if (size(vxyzu(:,i)) >= 4) then
      vxyzu(4,i) = 100. * factor
    end if
 enddo


end subroutine setpart

end module setup
