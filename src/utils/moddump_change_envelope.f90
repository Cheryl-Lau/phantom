!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module moddump
!
! Routine to set up a new envelope 
!
! :References: None
!
! :Owner: Cheryl Lau
!
! :Runtime parameters: 
!
! :Dependencies: None
!
 implicit none

 real :: r_cloud = 24.5 
 integer :: npart_envelope = 3E5
 real :: rho_envelope_cgs = 4e-25
 real :: u_envelope_cgs = 9.6087E+10

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,       only:igas,kill_particle,hfact 
 use physcon,    only:pi
 use io,         only:fatal
 use units,      only:unit_density,unit_ergg
 use partinject, only:add_or_update_particle
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: ip,np_outer,iadd 
 real    :: ran_x,ran_y,ran_z 
 real    :: totmass,vol_envelope,rad_envelope,u_envelope,hguess,r,xyz(3),vxyz(3),x,y,z 

 print*,'npart before delete',npart 

 do ip = 1,npart 
    r = sqrt(dot_product(xyzh(1:3,ip),xyzh(1:3,ip)))
    if (r > r_cloud) then 
       call kill_particle(ip) 
    endif 
 enddo 

 hguess = hfact*(massoftype(igas)/abs(rho_envelope_cgs/unit_density))**(1.d0/3.d0)
 u_envelope = u_envelope_cgs/unit_ergg

 totmass = npart_envelope *massoftype(igas)
 vol_envelope = totmass/(rho_envelope_cgs/unit_density)
 rad_envelope = (vol_envelope/(4./3.*pi))**(1./3.)
 
 iadd = npart 
 np_outer = 0 
 do while (np_outer <= npart_envelope)
    call random_number(ran_x)
    call random_number(ran_y)
    call random_number(ran_z)
    x = (2.*ran_x-1)*rad_envelope
    y = (2.*ran_y-1)*rad_envelope
    z = (2.*ran_z-1)*rad_envelope
    if (sqrt(x**2 + y**2 + z**2) < rad_envelope) then 
       np_outer = np_outer + 1 
       if (sqrt(x**2 + y**2 + z**2) > r_cloud) then 
           iadd = iadd + 1 
           xyz = (/ x, y, z /)
           vxyz = (/ 0., 0., 0. /)
           call add_or_update_particle(igas,xyz,vxyz,hguess,u_envelope,iadd,npart,&
                                       npartoftype,xyzh,vxyzu)
       endif 
    endif 
 enddo 

 npartoftype(:) = 0
 npartoftype(igas) = npart

 print*,'npart after reset',npart 

 return
end subroutine modify_dump

end module moddump
