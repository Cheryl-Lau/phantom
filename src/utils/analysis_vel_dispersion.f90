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
 character(len=20), parameter, public :: analysistype = 'vel_dispersion'

 public :: do_analysis

 private

 integer, parameter :: nrad  = 10
 real    :: rad_thresh(nrad) = (/ 6e-2,8e-2,1e-1,2e-1,4e-1,6e-1,8e-1,1e0,2e0,4e0 /)
 real    :: rholimit_cgs     = 1e-22
 logical :: only_highdenpart = .true. 

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use linklist, only:node,ifirstincell,listneigh,set_linklist
 use kdtree,   only:getneigh
 use dim,      only:maxneigh
 use units,    only:udist,utime,umass,unit_velocity,unit_density,unit_pressure,unit_ergg
 use io,       only:fatal,warning 
 use part,     only:hfact,rhoh,massoftype,igas
 use eos,      only:gamma
 use physcon,  only:gg 
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer, parameter :: neighcachesize = 1E5
 integer :: ip,ineigh,nneigh,ixyzcachesize,n,irad,ip_neigh
 real    :: xyzcache(3,neighcachesize)
 real    :: pmass,maxrad,dist2,rho,rad2_limit,mean_v,sigma_v
 real    :: sigma_v_allrad(nrad),virial_term_allrad(nrad)
 real,   allocatable :: dumxyzh(:,:)
 character(len=70) :: filename

 !- Particle mass
 pmass = massoftype(igas)

 !- Max radius 
 maxrad = rad_thresh(nrad)

 !- Build tree 
 allocate(dumxyzh(4,npart))
 dumxyzh = xyzh
 call set_linklist(npart,npart,dumxyzh,vxyzu)

 open(unit=2026,file='velocity_dispersion_'//TRIM(dumpfile)//'.dat',status='replace')
 open(unit=2027,file='virial_term_'//TRIM(dumpfile)//'.dat',status='replace')


 !- Loop over each length-scale 
 !$omp parallel do default(none) shared(npart,pmass,xyzh,vxyzu,rad_thresh,rholimit_cgs) &
 !$omp shared(node,hfact,maxrad,xyzcache,ifirstincell,unit_density,only_highdenpart,unit_velocity) &
 !$omp shared(sigma_v_allrad,virial_term_allrad,umass,udist) &
 !$omp private(ip,rho,nneigh,irad,rad2_limit,mean_v,sigma_v,n,ineigh,ip_neigh,dist2) &
 !$omp schedule(runtime)
 over_part: do ip = 1,npart
    print*,'ip/npart',ip,'/',npart
    rho = pmass*(hfact/abs(xyzh(4,ip)))**3

    if (only_highdenpart .and. rho < rholimit_cgs/unit_density) cycle over_part

    ! Get all neigh 
    call getneigh(node,xyzh(1:3,ip),0.,maxrad,3,listneigh,nneigh,xyzh,xyzcache,neighcachesize,ifirstincell,.false.)


    over_rad: do irad = 1,nrad
       rad2_limit = (rad_thresh(irad))**2 

       !- Find mean velocity 
       mean_v = 0. 
       n = 0
       over_neigh_mean: do ineigh = 1,nneigh
          ip_neigh = listneigh(ineigh)
          dist2 = mag2(xyzh(1:3,ip)-xyzh(1:3,ip_neigh))
          if (dist2 < rad2_limit) then 
             n = n + 1 
             mean_v = mean_v + sqrt(mag2(vxyzu(1:3,ip_neigh)))
          endif 
       enddo over_neigh_mean
       if (n == 0) cycle over_rad
       mean_v = mean_v/real(n) 

       !- Find dispersion (std)
       sigma_v = 0.
       over_neigh_sigma: do ineigh = 1,nneigh
          ip_neigh = listneigh(ineigh)
          dist2 = mag2(xyzh(1:3,ip)-xyzh(1:3,ip_neigh)) 
          if (dist2 < rad2_limit) then 
            sigma_v = sigma_v + (sqrt(mag2(vxyzu(1:3,ip_neigh))) - mean_v)**2 
          endif 
       enddo over_neigh_sigma
       sigma_v = sqrt(sigma_v/real(n)) 

       sigma_v_allrad(irad) = sigma_v*unit_velocity
       if (sigma_v > 0) then   
          virial_term_allrad(irad) = 2.d0*(n*pmass)/(sigma_v**2*rad_thresh(irad))
       else 
          virial_term_allrad(irad) = -1 ! to be removed during analysis 
       endif 
    enddo over_rad

    !$omp critical 
    write(2026,'(1i15,10f20.10)') ip, sigma_v_allrad(1:nrad)
    write(2027,'(1i15,10f20.10)') ip, virial_term_allrad(1:nrad)
    !$omp end critical 
 enddo over_part
 !$omp end parallel do

 close(2026)
 close(2027)
 deallocate(dumxyzh)

end subroutine do_analysis


real function mag2(vec)
 real,   intent(in) :: vec(3)

 mag2 = dot_product(vec,vec)

end function mag2

!--------------------------------------------------------------------------
end module analysis
