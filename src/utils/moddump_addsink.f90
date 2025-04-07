!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module moddump
!
! Add sink particle
!
! :References: None
!
! :Owner: Cheryl Lau
!
! :Runtime parameters: None
!
! :Dependencies: part
!

 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,         only:igas,xyzmh_ptmass,vxyz_ptmass,nptmass,ihsoft,ihacc 
 use centreofmass, only:reset_centreofmass
 implicit none
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: isink 

 nptmass = 1
 isink = 1
 xyzmh_ptmass(1:3,isink) = (/0.,0.,0./)
 xyzmh_ptmass(4,isink)   = 50. 
 xyzmh_ptmass(ihsoft,isink) = 0.005
 xyzmh_ptmass(ihacc,isink) = 0.005

 npartoftype(:) = 0
 npartoftype(igas) = npart

 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

 return
end subroutine modify_dump

end module moddump
