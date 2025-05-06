!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine which computes total kinetic energy, thermal energy and 
! radial momentum for each longitude from a source. 
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
 character(len=20), parameter, public :: analysistype = 'theta_energy'

 public :: do_analysis

 private

 integer, parameter :: ntheta = 360/1
 real    :: maxr = 7.
 integer :: isink_src = 14 
 real    :: xyz_src_in(3) = (/ 0.,0.,0. /)
 logical :: use_sink = .false. 


contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use units,    only:udist,utime,umass,unit_velocity,unit_density,unit_ergg
 use io,       only:fatal,warning 
 use part,     only:hfact,rhoh,massoftype,igas
 use part,     only:igas,xyzmh_ptmass,vxyz_ptmass,nptmass,ihsoft,ihacc 
 use eos,      only:gamma
 use physcon,  only:gg,pi
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer :: ip,i,itheta
 real    :: xyz_src(3),vxyz_src(3),dr(3)
 real    :: totmass(ntheta),totekin(ntheta),totetherm(ntheta),totmomen(ntheta)
 real    :: theta_min,theta_max,theta,dtheta
 real    :: mintheta,maxtheta,minphi,maxphi 
 real    :: absdr,dx,dy,dz,pmass,absdv,r_unitvec(3),radvel
 real    :: theta_lowbound(ntheta),theta_upbound(ntheta)
 character(len=70) :: filename

 !- Particle mass
 pmass = massoftype(igas)

 !- Set centre 
 if (use_sink) then 
    if (isink_src > nptmass) stop 'sink not found'
    xyz_src  = xyzmh_ptmass(1:3,isink_src)
    vxyz_src = (/ 0.,0.,0. /)  !vxyz_ptmass(1:3,isink_src)
    print*,'xyz_src',xyz_src,vxyz_src
 else 
    xyz_src  = xyz_src_in(1:3)
    vxyz_src = (/ 0.,0.,0. /)
 endif 

 !- Set range 
 theta_min = -pi
 theta_max = pi
 dtheta    = (theta_max-theta_min)/ntheta 

 do itheta = 1,ntheta 
    theta_lowbound(itheta) = theta_min + (itheta-1)*dtheta 
    theta_upbound(itheta) = theta_min + (itheta)*dtheta 
 enddo 

 !- Init storage 
 totmass   = 0.
 totekin   = 0. 
 totetherm = 0. 
 totmomen  = 0. 

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

    theta = atan2(dy,dx)

    mintheta = min(mintheta,theta)  !- record range of theta
    maxtheta = max(maxtheta,theta)

    itheta = floor((theta-theta_min)/dtheta)+1
    if (itheta < 1 .or. itheta > ntheta) then 
       print*, itheta, theta, dx, dy, dz, absdr 
       stop 'wrong itheta'
    endif 

    !- Accumulate mass 
    totmass(itheta) = totmass(itheta) + pmass*umass

    !- Accumulate kinetic energy 
    absdv = sqrt(mag2(vxyzu(1:3,ip)-vxyz_src(1:3)))
    totekin(itheta) = totekin(itheta) + 5d-1*pmass*umass*(absdv*unit_velocity)**2

    !- Accumulate thermal energy 
    totetherm(itheta) = totetherm(itheta) + pmass*umass*vxyzu(4,ip)*unit_ergg

    !- Accumulate radial momentum 
    r_unitvec = dr/absdr
    radvel = dot_product(vxyzu(1:3,ip),r_unitvec)
    totmomen(itheta) = totmomen(itheta) + pmass*umass*radvel*unit_velocity

 enddo over_parts 

 print*,'theta range',mintheta,maxtheta


 !- Write array to file 
 open(unit=2026,file='theta_energies/theta_mass_'//TRIM(dumpfile)//'.dat',status='replace')
 do itheta = 1,ntheta
    write(2026,*) totmass(itheta)
 enddo 
 open(unit=2027,file='theta_energies/theta_ekin_'//TRIM(dumpfile)//'.dat',status='replace')
 do itheta = 1,ntheta
    write(2027,*) totekin(itheta)
 enddo 
 open(unit=2028,file='theta_energies/theta_etherm_'//TRIM(dumpfile)//'.dat',status='replace')
 do itheta = 1,ntheta
    write(2028,*) totetherm(itheta)
 enddo 
 open(unit=2029,file='theta_energies/theta_momen_'//TRIM(dumpfile)//'.dat',status='replace')
 do itheta = 1,ntheta
    write(2029,*) totmomen(itheta)
 enddo 

end subroutine do_analysis


real function mag2(vec)
 real,   intent(in) :: vec(3)

 mag2 = dot_product(vec,vec)

end function mag2

!--------------------------------------------------------------------------
end module analysis
