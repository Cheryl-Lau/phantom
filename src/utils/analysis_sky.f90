!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine which computes total mass, kinetic energy and thermal energy
! per unit surface area for each longitude and latitude seen from a source. 
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
 character(len=20), parameter, public :: analysistype = 'sky'

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
 use units,    only:udist,utime,umass,unit_velocity,unit_density,unit_ergg,unit_energ
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
 real    :: xyz_src(3),vxyz_src(3),dr(3)
 real    :: totmass(ntheta,nphi),totekin(ntheta,nphi),totetherm(ntheta,nphi),totmomen(ntheta,nphi)
 real    :: theta_min,theta_max,theta,dtheta,phi_min,phi_max,phi,dphi
 real    :: mintheta,maxtheta,minphi,maxphi 
 real    :: absdr,dx,dy,dz,pmass,vol_element,area_element,absdv
 real    :: phi_lowbound(nphi),phi_upbound(nphi),theta_lowbound(ntheta),theta_upbound(ntheta)
 real,   allocatable :: dumxyzh(:,:)
 character(len=70) :: filename

 !- Particle mass
 pmass = massoftype(igas)

 !- Set centre 
 if (use_sink) then 
    if (isink_src > nptmass) stop 'sink not found'
    xyz_src  = xyzmh_ptmass(1:3,isink_src)
    vxyz_src = vxyz_ptmass(1:3,isink_src)
    print*,'xyz_src',xyz_src,vxyz_src
 else 
    xyz_src  = xyz_src_in(1:3)
    vxyz_src = (/ 0.,0.,0. /)
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
 totmass   = 0.
 totekin   = 0. 
 totmomen  = 0.
 totetherm = 0. 

 minphi = huge(phi)
 maxphi = tiny(phi)
 mintheta = huge(theta)
 maxtheta = tiny(theta) 

 over_parts: do ip = 1,npart 
    dr = xyzh(1:3,ip)-xyz_src(1:3)
    absdr = sqrt(mag2(dr))
    if (absdr > maxr) cycle 
    dx = min(dr(1)/absdr,1.d0)  !- limit unitvec to [-1,+1]
    dx = max(dr(1)/absdr,-1.d0)
    dy = min(dr(2)/absdr,1.d0)
    dy = max(dr(2)/absdr,-1.d0)
    dz = min(dr(3)/absdr,1.d0)
    dz = max(dr(3)/absdr,-1.d0)
    phi   = dasin(dz)
    theta = atan2(dy,dx)

    minphi = min(minphi,phi)   !- record range of theta & phi
    maxphi = max(maxphi,phi)
    mintheta = min(mintheta,theta)
    maxtheta = max(maxtheta,theta)

    iphi = floor((phi-phi_min)/dphi)+1 
    itheta = floor((theta-theta_min)/dtheta)+1
    if (iphi < 1 .or. iphi > nphi) then
       print*, iphi, phi, dx, dy, dz, absdr
       stop 'wrong iphi'
    elseif (itheta < 1 .or. itheta > ntheta) then 
       print*, itheta, theta, dx, dy, dz, absdr 
       stop 'wrong itheta'
    endif 

    !- Area bounded by this phi & theta cell 
!    vol_element = 1./3.*maxr**3 * (cos(pi/2.-phi_upbound(iphi))-cos(pi/2.-phi_lowbound(iphi))) &
!                  &* (theta_upbound(itheta)-theta_lowbound(itheta))
    area_element = maxr**2 * (cos(pi/2.-phi_upbound(iphi))-cos(pi/2.-phi_lowbound(iphi))) &
                  &* (theta_upbound(itheta)-theta_lowbound(itheta))

    !- Accumulate mass 
    totmass(itheta,iphi) = totmass(itheta,iphi) + pmass*umass/abs(area_element*udist**2)  

    !- Accumulate kinetic energy and momentum
    absdv = sqrt(mag2(vxyzu(1:3,ip)-vxyz_src(1:3)))
    totekin(itheta,iphi) = totekin(itheta,iphi) + 5d-1*pmass*umass*(absdv*unit_velocity)**2/abs(area_element*udist**2)
    totmomen(itheta,iphi) = totmomen(itheta,iphi) + pmass*umass*(absdv*unit_velocity)/abs(area_element*udist**2)

    !- Accumulate thermal energy 
    totetherm(itheta,iphi) = totetherm(itheta,iphi) + vxyzu(4,ip)*pmass*unit_energ/abs(area_element*udist**2)

 enddo over_parts 

 print*,'theta range',mintheta,maxtheta
 print*,'phi range',minphi,maxphi


 !- Write 2D array to file 
 open(unit=2026,file='sky_density_'//TRIM(dumpfile)//'.dat',status='replace')
 do iphi = 1,nphi
    write(2026,*) totmass(1:ntheta,iphi)
 enddo 
 open(unit=2027,file='sky_ekin_'//TRIM(dumpfile)//'.dat',status='replace')
 do iphi = 1,nphi
    write(2027,*) totekin(1:ntheta,iphi)
 enddo 
 open(unit=2028,file='sky_momen_'//TRIM(dumpfile)//'.dat',status='replace')
 do iphi = 1,nphi
    write(2028,*) totmomen(1:ntheta,iphi)
 enddo 
 open(unit=2029,file='sky_etherm_'//TRIM(dumpfile)//'.dat',status='replace')
 do iphi = 1,nphi
    write(2029,*) totetherm(1:ntheta,iphi)
 enddo 


end subroutine do_analysis


real function mag2(vec)
 real,   intent(in) :: vec(3)

 mag2 = dot_product(vec,vec)

end function mag2

!--------------------------------------------------------------------------
end module analysis
