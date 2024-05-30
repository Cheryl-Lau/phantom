!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module moddump
!
! Routine to manually carve out a conic channel in a spherical HII region
! assumming the HII region is centered at the origin
!
! :References: None
!
! :Owner: Cheryl Lau
!
! :Runtime parameters: 
!   - rad_chnl          : *maximum radius of the channell to clear*
!   - omegafrac_chnl    : *size of channell as a fraction of sphere's total solid angle*
!   - pfrac_chnl        : *fraction of particles to remove within channel*
!
! :Dependencies: None
!
 implicit none

 real   :: rad_chnl = 10.05
 real   :: omegafrac_chnl = 0.1
 real   :: pfrac_chnl = 0.9999

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use random,  only:ran2
 use part,    only:igas
 use physcon, only:pi
 use io,      only:fatal
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,    allocatable   :: xyzh_store(:,:),vxyzu_store(:,:)
 integer :: maxpart,npart_store,ip
 integer :: iseed = -12345
 real    :: omega,omega_chnl,x,y,z,r,radius,rad_circ,area
 logical :: add_particle

 omega_chnl = omegafrac_chnl * 4.*pi

 maxpart = npart
 allocate(xyzh_store(4,maxpart))
 allocate(vxyzu_store(4,maxpart))
 xyzh_store(:,:)  = 0. 
 vxyzu_store(:,:) = 0.
 npart_store = 0

 !- Go through each particle and only store those which are outside the channel
 do ip = 1,npart 
    x = xyzh(1,ip)
    y = xyzh(2,ip)
    z = xyzh(3,ip)
    r = (x**2 + y**2 + z**2)**(1./2.) 
    !- Checks 
    add_particle = .true. 
    if (x > 0. .and. r < rad_chnl) then 
       radius   = x 
       rad_circ = (abs(y)**2 + abs(z)**2)**(1./2.)
       area     = 4.*pi*rad_circ**2
       omega    = area/radius**2
       if (omega < omega_chnl) then
          if (ran2(iseed) < pfrac_chnl) add_particle = .false.
       endif
    endif
    if (add_particle) then 
       npart_store = npart_store + 1
       xyzh_store(:,npart_store)  = xyzh(:,ip)
       vxyzu_store(:,npart_store) = vxyzu(:,ip)
    endif 
 enddo 
 
 print*,'Number of particles before and after carving: ',npart,npart_store
 if (npart == npart_store) call fatal('moddump_clear_channel','no particles removed')

 npart = npart_store
 xyzh(:,:)  = xyzh_store(:,:)
 vxyzu(:,:) = vxyzu_store(:,:)

 npartoftype(:) = 0
 npartoftype(igas) = npart

 return
end subroutine modify_dump

end module moddump
