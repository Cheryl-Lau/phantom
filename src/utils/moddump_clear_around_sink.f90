!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module moddump
!
! Slightly push away the particles around a sink to ease feedback 
!
! :References: Blomme, R, 1990, A&A, 229, 513-526 (https://ui.adsabs.harvard.edu/abs/1990A%26A...229..513B/abstract)
!
! :Owner: Cheryl Lau
!
! :Runtime parameters: None
!
! :Dependencies: part
!

 implicit none

 integer :: isink = 6            ! index of target sink 
 real    :: rad_max  = 0.1        ! threshold radius 
 real    :: vel_max_cgs = 1e7 !3e8     ! stellar wind vel (Blomme 1990)

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,         only:igas,xyzmh_ptmass,vxyz_ptmass,nptmass,ihsoft,ihacc 
 use units,        only:unit_velocity 
 use centreofmass, only:reset_centreofmass
 implicit none
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: ip
 real    :: xyz_sink(3),dr(3),dr2,rad2_max,vel_max,absv,absdr,r_fac

 xyz_sink = xyzmh_ptmass(1:3,isink)
 rad2_max = rad_max**2
 vel_max = vel_max_cgs/unit_velocity 

 do ip = 1,npart 
    dr = xyzh(1:3,ip)-xyz_sink(1:3)
    dr2 = mag2(dr)
    if (dr2 < rad2_max) then
       r_fac = dr2/rad2_max 
       absv = (1.-r_fac)*vel_max
       absdr = sqrt(dr2)
       vxyzu(1:3,ip) = absv * dr/absdr 
    endif 
 enddo 

 npartoftype(:) = 0
 npartoftype(igas) = npart

 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

 return
end subroutine modify_dump


real function mag2(vec)
 real,   intent(in) :: vec(3)

 mag2 = dot_product(vec,vec)

end function mag2



end module moddump
