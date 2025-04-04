!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine which computes total mass (column density) along longitudes 
! and latitues 
!
! :References: None
!
! :Owner: Cheryl Lau
!
! :Runtime parameters: None
!
! :Dependencies:
!
 implicit none
 character(len=20), parameter, public :: analysistype = 'column_density'

 public :: do_analysis

 private

 integer, parameter :: ntheta = 360/1
 integer, parameter :: nphi = 180/1
 real    :: maxr = 15.
 integer :: isink_src = 14 
 real    :: xyz_src_in(3) = (/ 0.,0.,0. /)
 logical :: use_sink = .true. 


contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use linklist, only:node,ifirstincell,listneigh,set_linklist
 use kdtree,   only:getneigh
 use dim,      only:maxneigh
 use units,    only:udist,utime,umass,unit_velocity,unit_density,unit_pressure,unit_ergg
 use io,       only:fatal,warning 
 use part,     only:hfact,rhoh,massoftype,igas
 use part,     only:igas,xyzmh_ptmass,vxyz_ptmass,nptmass,ihsoft,ihacc 
 use eos,      only:gamma
 use physcon,  only:gg,pi
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer, parameter :: neighcachesize = 1E5
 integer :: ip,i,nneigh,ixyzcachesize,itheta,iphi
 real    :: xyzcache(3,neighcachesize)
 real    :: xyz_src(3),totmass(ntheta,nphi)
 real    :: theta_min,theta_max,theta,dtheta,phi_min,phi_max,phi,dphi
 real    :: mintheta,maxtheta,minphi,maxphi 
 real    :: r,x,y,z,pmass,vol_element,area_element
 real    :: phi_lowbound(nphi),phi_upbound(nphi),theta_lowbound(ntheta),theta_upbound(ntheta)
 real,   allocatable :: dumxyzh(:,:)
 character(len=70) :: filename

 !- Particle mass
 pmass = massoftype(igas)

 !- Set centre 
 if (use_sink) then 
    if (isink_src > nptmass) stop 'sink not found'
    xyz_src = xyzmh_ptmass(1:3,isink_src)
 else 
    xyz_src = xyz_src_in(1:3)
 endif 
 print*,'Getting column density w.r.t.: ',xyz_src

 !- Set range 
 theta_min = -pi
 theta_max = pi
 dtheta    = (theta_max-theta_min)/ntheta 
 phi_min   = -pi/2.d0
 phi_max   = pi/2.d0 
 dphi      = (phi_max-phi_min)/nphi
 
 do iphi = 1,nphi
    phi_lowbound(iphi) = phi_min + (iphi-1)*dphi
    phi_upbound(iphi) = phi_min + (iphi)*dphi
 enddo 
 do itheta = 1,ntheta 
    theta_lowbound(itheta) = theta_min + (itheta-1)*dtheta 
    theta_upbound(itheta) = theta_min + (itheta)*dtheta 
 enddo 

 !- Init storage 
 totmass = 0.

 minphi = huge(phi)
 maxphi = tiny(phi)
 mintheta = huge(theta)
 maxtheta = tiny(theta) 

 over_parts: do ip = 1,npart 
    r = sqrt(mag2(xyzh(1:3,ip)-xyz_src(1:3)))
    if (r > maxr) cycle 
    x = min(xyzh(1,ip)/r,1.d0)
    y = min(xyzh(2,ip)/r,1.d0)
    z = min(xyzh(3,ip)/r,1.d0)
    phi   = dasin(z)
    theta = atan2(y,x)

    minphi = min(minphi,phi)
    maxphi = max(maxphi,phi)
    mintheta = min(mintheta,theta)
    maxtheta = max(maxtheta,theta)

    iphi = floor((phi-phi_min)/dphi)+1 
    itheta = floor((theta-theta_min)/dtheta)+1
    if (iphi < 1 .or. iphi > nphi) then
       print*, iphi, phi, x, y, z
       stop 'wrong iphi'
    elseif (itheta < 1 .or. itheta > ntheta) then 
       print*, itheta, theta, x, y, z
       stop 'wrong itheta'
    endif 

    !- Accumulate mass 
!    vol_element = 1./3.*maxr**3 * (cos(pi/2.-phi_upbound(iphi))-cos(pi/2.-phi_lowbound(iphi))) &
!                  &* (theta_upbound(itheta)-theta_lowbound(itheta))
    area_element = maxr**2 * (cos(pi/2.-phi_upbound(iphi))-cos(pi/2.-phi_lowbound(iphi))) &
                  &* (theta_upbound(itheta)-theta_lowbound(itheta))
    totmass(itheta,iphi) = totmass(itheta,iphi) + pmass*umass/abs(area_element*udist**2)  ! weighted by area bounded by phi & theta range

 enddo over_parts 

 print*,'theta range',mintheta,maxtheta
 print*,'phi range',minphi,maxphi


 !- Write 2D array to file 
 open(unit=2026,file='theta_phi_'//TRIM(dumpfile)//'.dat',status='replace')
 do iphi = 1,nphi
    write(2026,*) totmass(1:ntheta,iphi)
 enddo 


end subroutine do_analysis


real function mag2(vec)
 real,   intent(in) :: vec(3)

 mag2 = dot_product(vec,vec)

end function mag2

!--------------------------------------------------------------------------
end module analysis
