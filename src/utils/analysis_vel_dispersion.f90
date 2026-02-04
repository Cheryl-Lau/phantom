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

 integer, parameter :: nrad  = 50
 real    :: rad_min = 1.d-2
 real    :: rad_max = 3.d+1
 real    :: rholimit_cgs     = 1.d-22
 logical :: only_highdenpart = .false. 

 integer :: isink_centre   = 1
 real    :: box_radius     = 70.
 real    :: box_centre(3)  = (/ 0.,0.,0. /)
 logical :: use_sink       = .false. 
 logical :: only_inbox     = .true.

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use linklist, only:node,ifirstincell,listneigh,set_linklist
 use kdtree,   only:getneigh
 use dim,      only:maxneigh
 use units,    only:udist,utime,umass,unit_velocity,unit_density,unit_pressure,unit_ergg
 use io,       only:fatal,warning 
 use part,     only:hfact,rhoh,massoftype,igas
 use part,     only:xyzmh_ptmass,vxyz_ptmass,nptmass
 use eos,      only:gamma
 use physcon,  only:gg,pi
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer, parameter :: neighcachesize = 1E5
 integer :: ip,ineigh,nneigh,ixyzcachesize,n,irad,ip_neigh
 real    :: xyzcache(3,neighcachesize)
 real    :: pmass,dist2,rho,rad2_limit,mean_v,sigma_v,mean_rho
 real    :: logr_min,dlogr,rad 
 real    :: rad_thresh(nrad)
 real    :: sigma_v_allrad(nrad),virial_term_allrad(nrad),rho_avg_allrad(nrad)
 real    :: xmin_box,xmax_box,ymin_box,ymax_box,zmin_box,zmax_box,centre(3)
 real,   allocatable :: dumxyzh(:,:)
 character(len=70) :: filename

 !- Particle mass
 pmass = massoftype(igas)

 !- Set up rad_thresh(nrad) array
 logr_min = log(rad_min)
 dlogr    = log(rad_max/rad_min)/dble(nrad)
 do irad = 1,nrad 
    rad = exp(logr_min + dble(irad)*dlogr)
    rad_thresh(irad) = rad 
 enddo 
 print*,'measuring sigma at: ',rad_thresh(1:nrad)

 !- Build tree 
 allocate(dumxyzh(4,npart))
 dumxyzh = xyzh
 call set_linklist(npart,npart,dumxyzh,vxyzu)

 !- Set boundaries if requested 
 if (only_inbox) then 
    if (use_sink) then
       if (isink_centre > nptmass) call fatal('analysis_velfield','sink not found.') 
       centre = xyzmh_ptmass(1:3,isink_centre)
    else
       centre = box_centre
    endif 
    xmin_box = centre(1) - box_radius; xmax_box = centre(1) + box_radius 
    ymin_box = centre(2) - box_radius; ymax_box = centre(2) + box_radius 
    zmin_box = centre(3) - box_radius; zmax_box = centre(3) + box_radius 
 endif 

 open(unit=2026,file='velocity_dispersion_'//TRIM(dumpfile)//'.dat',status='replace')
 open(unit=2027,file='virial_term_'//TRIM(dumpfile)//'.dat',status='replace')
 open(unit=2028,file='density_sizescale_'//TRIM(dumpfile)//'.dat',status='replace')

 !- Write header line - the radii
 write(2026,'(15x,100f30.10)') rad_thresh(1:nrad)
 write(2027,'(15x,100f30.10)') rad_thresh(1:nrad)
 write(2028,'(15x,100f30.10)') rad_thresh(1:nrad)


 !- Loop over each length-scale 
 !$omp parallel do default(none) shared(npart,pmass,xyzh,vxyzu,rad_thresh,rholimit_cgs) &
 !$omp shared(node,hfact,rad_max,xyzcache,ifirstincell,only_highdenpart) &
 !$omp shared(only_inbox,xmin_box,xmax_box,ymin_box,ymax_box,zmin_box,zmax_box) &
 !$omp shared(umass,udist,unit_velocity,unit_density) &
 !$omp private(ip,rho,nneigh,irad,rad2_limit,mean_v,sigma_v,mean_rho,n,ineigh,ip_neigh,dist2) &
 !$omp private(sigma_v_allrad,virial_term_allrad,rho_avg_allrad) &
 !$omp schedule(runtime)
 over_part: do ip = 1,npart
    print*,'ip/npart',ip,'/',npart
    rho = pmass*(hfact/abs(xyzh(4,ip)))**3

    if (only_highdenpart .and. rho < rholimit_cgs/unit_density) cycle over_part

    if (only_inbox) then 
       if (xyzh(1,ip) < xmin_box .or. xyzh(1,ip) > xmax_box) cycle over_part
       if (xyzh(2,ip) < ymin_box .or. xyzh(2,ip) > ymax_box) cycle over_part
       if (xyzh(3,ip) < zmin_box .or. xyzh(3,ip) > zmax_box) cycle over_part
    endif 

    over_rad: do irad = 1,nrad
       rad2_limit = (rad_thresh(irad))**2 

       ! Get all neigh 
       call getneigh(node,xyzh(1:3,ip),0.,rad_thresh(irad),3,listneigh,nneigh,xyzh,xyzcache,neighcachesize,ifirstincell,.false.)

       !- Find mean velocity and density 
       mean_v = 0. 
       mean_rho = 0.
       n = 0
       over_neigh_mean: do ineigh = 1,nneigh
          ip_neigh = listneigh(ineigh)
          dist2 = mag2(xyzh(1:3,ip)-xyzh(1:3,ip_neigh))
          if (dist2 < rad2_limit) then 
             n = n + 1 
             mean_v = mean_v + sqrt(mag2(vxyzu(1:3,ip_neigh)))
             mean_rho = mean_rho + pmass*(hfact/abs(xyzh(4,ip_neigh)))**3
          endif 
       enddo over_neigh_mean
       if (n == 0) cycle over_rad
       mean_v = mean_v/real(n) 
       mean_rho = mean_rho/real(n) 

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

       !- Virial parameter 
       if (sigma_v > 0) then   
          virial_term_allrad(irad) = 2.d0*gg*(n*pmass*umass)/((sigma_v*unit_velocity)**2*rad_thresh(irad)*udist)
       else 
          virial_term_allrad(irad) = -1 ! to be removed during analysis 
       endif 

       !- Mean density 
       rho_avg_allrad(irad) = (n*pmass*umass)/(4./3.*pi*(rad_thresh(irad)*udist)**3)

    enddo over_rad

    !$omp critical 
    write(2026,'(1i15,100f30.10)') ip, sigma_v_allrad(1:nrad)
    write(2027,'(1i15,100f30.10)') ip, virial_term_allrad(1:nrad)
    write(2028,'(1i15,100es30.10)') ip, rho_avg_allrad(1:nrad)
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
