!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine which measures the velocity dispersion of all/high-density particles
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
 character(len=40), parameter, public :: analysistype = 'particle_properties_from_ionizing_source'

 public :: do_analysis

 private

! real    :: xyz_src(3) = (/ 2.18E+00,2.88E+00,-2.58E+00 /)
 integer :: isink_src = 139
 integer, parameter :: nrad  = 8
 real    :: rad_thresh(nrad) = (/ 6e-2,1e-1,2e-1,4e-1,6e-1,1e0,2e0,4e0 /)
 real    :: rholimit_cgs  = 1e-21
 real    :: templimit     = 8000 
 logical :: only_highrho  = .false. 
 logical :: only_hightemp = .true.

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use linklist, only:node,ifirstincell,listneigh,set_linklist
 use kdtree,   only:getneigh
 use dim,      only:maxneigh
 use units,    only:udist,utime,umass,unit_velocity,unit_density,unit_pressure,unit_ergg
 use io,       only:fatal,warning 
 use part,     only:hfact,rhoh,massoftype,igas
 use part,     only:xyzmh_ptmass,vxyz_ptmass,ihacc,nptmass
 use eos,      only:gamma,gmw 
 use physcon,  only:gg,pi,kboltz,mass_proton_cgs
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer, parameter :: neighcachesize = 1E5
 integer :: ip,ineigh,nneigh,ixyzcachesize,n,irad,ip_neigh
 real    :: xyzcache(3,neighcachesize)
 real    :: pmass,xyz_src(3),vxyz_src(3),maxrad,rho,temp,xb_xa(3),vb_va(3),r,vr_wrtsrc
 real    :: rad2,mean_v,dist2,sigma_v,sigma_v_allrad(nrad),dist_cgs,vr_cgs,rho_cgs
 real,   allocatable :: dumxyzh(:,:)
 character(len=100) :: filename

 !- Particle mass
 pmass = massoftype(igas)

 !- Build tree 
 allocate(dumxyzh(4,npart))
 dumxyzh = xyzh
 call set_linklist(npart,npart,dumxyzh,vxyzu)

 open(unit=2026,file='radvel_fromsrc_'//TRIM(dumpfile)//'.dat',status='replace')
 open(unit=2027,file='density_fromsrc_'//TRIM(dumpfile)//'.dat',status='replace')
 open(unit=2028,file='sigma_v_fromsrc_'//TRIM(dumpfile)//'.dat',status='replace')

 ! Properties of ionizing source 
 xyz_src  = xyzmh_ptmass(1:3,isink_src)
 vxyz_src = vxyz_ptmass(1:3,isink_src)
 print*,'sink info',xyz_src,vxyz_src

 ! Max radius for measuring sigma_v 
 maxrad = rad_thresh(nrad)


 !$omp parallel do default(none) shared(npart,pmass,xyzh,vxyzu,rad_thresh,rholimit_cgs) &
 !$omp shared(gamma,gmw,templimit,xyz_src,vxyz_src) &
 !$omp shared(node,hfact,maxrad,xyzcache,ifirstincell,only_highrho,only_hightemp) &
 !$omp shared(umass,udist,unit_velocity,unit_density,unit_ergg) &
 !$omp private(ip,rho,temp,nneigh,irad,rad2,mean_v,sigma_v,n,ineigh,ip_neigh,dist2) &
 !$omp private(sigma_v_allrad,xb_xa,vb_va,r,vr_wrtsrc,dist_cgs,vr_cgs,rho_cgs) &
 !$omp schedule(runtime)
 over_part: do ip = 1,npart
    print*,'ip/npart',ip,'/',npart
    rho  = pmass*(hfact/abs(xyzh(4,ip)))**3
    temp = vxyzu(4,ip)/kboltz*(gmw*mass_proton_cgs*(gamma-1.))*unit_ergg

    if (only_highrho  .and. rho  < rholimit_cgs/unit_density) cycle over_part
    if (only_hightemp .and. temp < templimit)                 cycle over_part

    ! Radial velocity of particle relative to ionizing source 
    xb_xa = xyzh(1:3,ip) - xyz_src(1:3)
    vb_va = vxyzu(1:3,ip) - vxyz_src(1:3)
    r = sqrt(mag2(xb_xa))
    vr_wrtsrc = dot_product(xb_xa,vb_va)/r

    ! Get neighbours 
    call getneigh(node,xyzh(1:3,ip),0.,maxrad,3,listneigh,nneigh,xyzh,xyzcache,neighcachesize,ifirstincell,.false.)
    print*,'getneigh'

    ! Velocity dispersion at each radii
    over_rad: do irad = 1,nrad 
       rad2 = (rad_thresh(irad))**2 
      
       !- Find mean velocity 
       mean_v = 0. 
       n = 0
       do ineigh = 1,nneigh
          ip_neigh = listneigh(ineigh)
          dist2 = mag2(xyzh(1:3,ip)-xyzh(1:3,ip_neigh))
          if (dist2 < rad2) then 
             n = n + 1 
             mean_v = mean_v + sqrt(mag2(vxyzu(1:3,ip_neigh)))
          endif 
       enddo 
       if (n == 0) cycle over_rad
       mean_v = mean_v/real(n) 

       !- Find dispersion (std)
       sigma_v = 0.
       do ineigh = 1,nneigh
          ip_neigh = listneigh(ineigh)
          dist2 = mag2(xyzh(1:3,ip)-xyzh(1:3,ip_neigh)) 
          if (dist2 < rad2) then 
             sigma_v = sigma_v + (sqrt(mag2(vxyzu(1:3,ip_neigh))) - mean_v)**2 
          endif 
       enddo 
       sigma_v = sqrt(sigma_v/real(n)) 

       sigma_v_allrad(irad) = sigma_v*unit_velocity
    enddo over_rad
    print*,'done over rad'

    dist_cgs = r *udist 
    vr_cgs   = vr_wrtsrc *unit_velocity 
    rho_cgs  = rho *unit_density 
    print*,'prop',dist_cgs,vr_cgs,sigma_v_allrad(1)

    !$omp critical 
    write(2026,*) dist_cgs, vr_cgs
    write(2027,*) dist_cgs, rho_cgs
    write(2028,*) dist_cgs, sigma_v_allrad(1:nrad)
    !$omp end critical 

 enddo over_part
 !$omp end parallel do

 close(2026)
 close(2027)
 close(2028) 
 deallocate(dumxyzh)

end subroutine do_analysis


real function mag2(vec)
 real,   intent(in) :: vec(3)

 mag2 = dot_product(vec,vec)

end function mag2

!--------------------------------------------------------------------------
end module analysis
