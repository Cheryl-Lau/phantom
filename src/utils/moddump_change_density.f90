!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module moddump
!
! Routine to change the density within a given radius 
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

 real :: r_moddens = 4.8
 real :: rho_moddens_cgs = 1e-23 

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,       only:igas,kill_particle,hrho
 use partinject, only:add_or_update_particle
 use physcon,    only:pi
 use io,         only:fatal
 use units,      only:unit_density 
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: ip,np_toadd,np_placed,np_killed,iadd
 real    :: ran_x,ran_y,ran_z
 real    :: r,xyz(3),vxyz(3),x,y,z,totmass,umean,hguess,pmass,rho_moddens 

 ! convert to code units 
 rho_moddens = rho_moddens_cgs/unit_density 

 ! Estimate h with rho 
 pmass = massoftype(igas) 
 hguess = hrho(rho_moddens,pmass)
 print*,'h',hguess 

 !- Clear region 
 umean = 0. 
 np_killed = 0
 do ip = 1,npart 
    r = sqrt(dot_product(xyzh(1:3,ip),xyzh(1:3,ip)))
    if (r < r_moddens) then 
       np_killed = np_killed + 1 
       umean = umean + vxyzu(4,ip) 
       call kill_particle(ip) 
    endif 
 enddo 
 umean = umean/np_killed 

 !- Number of particles to place 
 totmass = rho_moddens * (4./3.)*pi*r_moddens**3 
 print*,'params',rho_moddens,pi,r_moddens**3
 print*,'totmass pmass',totmass, pmass 
 np_toadd = nint(totmass/pmass)

 print*,'np_killed,np_toadd',np_killed,np_toadd 

 np_placed = 0
 iadd = npart 
 do while (np_placed <= np_toadd)
    call random_number(ran_x)
    call random_number(ran_y)
    call random_number(ran_z)
    x = (2.*ran_x-1)*r_moddens
    y = (2.*ran_y-1)*r_moddens
    z = (2.*ran_z-1)*r_moddens 
    if (sqrt(x**2 + y**2 + z**2) < r_moddens) then 
        np_placed = np_placed + 1 
        iadd = iadd + 1 
        xyz = (/ x, y, z /)
        vxyz = (/ 0., 0., 0. /)
        call add_or_update_particle(igas,xyz,vxyz,hguess,umean,iadd,npart,&
                                    npartoftype,xyzh,vxyzu)
    endif 
 enddo

 npartoftype(:) = 0
 npartoftype(igas) = npart

 return
end subroutine modify_dump

end module moddump
